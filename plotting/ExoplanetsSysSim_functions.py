# To import required modules:
import numpy as np
import time
import sys
import matplotlib
import matplotlib.cm as cm #for color maps
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec #for specifying plot attributes
from matplotlib import ticker #for setting contour plots to log scale
import scipy.integrate #for numerical integration
import scipy.misc #for factorial function
from scipy.special import erf #error function, used in computing CDF of normal distribution
import corner #corner.py package for corner plots
#matplotlib.rc('text', usetex=True)





##### This module will be used to define almost ALL the functions required to load, analyze, and plot simulated data generated from the ExoplanetsSysSim models





# Useful fundamental constants:

AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Mearth = 5.972*10.**24 #Earth mass in kg
Rearth = 6.371*10.**8. #Earth radius in cm

# Fixed parameters:

N_Kep = 79935 #139232 #number of Kepler targets satisfying our cuts to give our observed catalog

res_ratios, res_width = [2.0, 1.5, 4/3., 5/4.], 0.05 #NOTE: in the model, the near-resonant planets have period ratios between X and (1+w)*X where X = [2/1, 3/2, 4/3, 5/4] and w = 0.05!

param_keys_all = [("num_targets_sim_pass_one", r'$N_{\rm stars,sim}$'),
                  ("max_incl_sys", r'$i_{\rm ref,max}$'),
                  ("log_rate_clusters", r'$\lambda_c$'),
                  ("max_clusters_in_sys", r'$N_{c,\rm max}$'),
                  ("log_rate_planets_per_cluster", r'$\lambda_p$'),
                  ("max_planets_in_clusters", r'$N_{p,\rm max}$'),
                  ("power_law_P", r'$\alpha_P$'),
                  ("min_period", r'$P_{\rm min}$'),
                  ("max_period", r'$P_{\rm max}$'),
                  ("power_law_r1", r'$\alpha_{R1}$'),
                  ("power_law_r2", r'$\alpha_{R2}$'),
                  ("break_radius (R_earth)", r'$R_{p,\rm break}$ $(R_\oplus)$'),
                  ("min_radius (R_earth)", r'$R_{p,\rm min}$ $(R_\oplus)$'),
                  ("max_radius (R_earth)", r'$R_{p,\rm max}$ $(R_\oplus)$'),
                  ("f_high_incl", r'$f_{\rm high,\sigma_i}$'),
                  ("sigma_incl", r'$\sigma_i$'),
                  ("sigma_incl_near_mmr", r'$\sigma_{i,\rm res}$'),
                  ("sigma_hk", r'$\sigma_e$'),
                  ("num_mutual_hill_radii", r'$\Delta_c$'),
                  ("sigma_log_radius_in_cluster", r'$\sigma_R$'),
                  ("sigma_logperiod_per_pl_in_cluster", r'$\sigma_N$')
                  ] #list of the symbols and names for all the model parameters; NOTE: although the params are named log rate of clusters and planets per cluster, we use the symbols and values for the rates





# Miscellaneous functions:

def a_from_P(P, Mstar):
    #This function converts period (days) to semi-major axis (AU) assuming mass of planet m << Mstar (Msun)
    y = (P/365.25)**(2./3.)*(Mstar/1.0)**(1./3.)
    return y

def P_from_a(a, Mstar):
    #This function converts semi-major axis (AU) to period (days) assuming mass of planet m << Mstar (Msun)
    y = 365.25*(a**(3./2.))*(Mstar/1.0)**(-1./2.)
    return y

def photoevap_boundary_Carrera2018(R, P):
    #R is the planet radius in Earth radii, P is the period in days
    #This function returns 1 if the planet is above the boundary, and 0 if the planet is below the boundary as defined by Eq. 5 in Carrera et al 2018
    Rtrans = 2.6*(P)**(-0.1467)
    if R > Rtrans:
        above_boundary = 1
    elif R < Rtrans:
        above_boundary = 0

    return above_boundary

def cdf_normal(x, mu=0., std=1.):
    # This function computes the CDF (i.e. the integral of the normal distribution between -inf and x) at x given mean 'mu' and standard deviation 'std'
    # Note: this function can deal with array inputs for x, mu and std, as long as the array inputs are the same shape
    return 0.5*(1. + erf((x - mu)/(std*np.sqrt(2))))





# Functions to load and analyze simulated physical catalogs:

def load_cat_phys(file_name):
    #This function loads the simulated physical catalog of planets:
    start = time.time()
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'physical_catalog_planets%s.txt' % run_number, 'r')
        lines = (line for line in file if not line.startswith('#'))
        #cat_phys = np.loadtxt(lines, skiprows=1, dtype={'names': ('target_id', 'star_id', 'planet_mass', 'planet_radius', 'period', 'ecc'), 'formats': ('i4', 'i4', 'f8', 'f8', 'f8', 'f8')})
        cat_phys = np.genfromtxt(lines, names=True, dtype=('i4', 'i4', 'f8', 'f8', 'f8', 'f8')) #faster than the above
    stop = time.time()
    print 'Time to load: %s s' % (stop - start)

    return cat_phys

def load_star_phys(file_name):
    #This function loads the catalog of stars with all simulated planets:
    start = time.time()
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'physical_catalog_stars%s.txt' % run_number, 'r')
        lines = (line for line in file if not line.startswith('#'))
        #star_phys = np.loadtxt(lines, skiprows=1, dtype={'names': ('target_id', 'star_id', 'star_mass', 'star_radius', 'num_planets'), 'formats': ('i4', 'i4', 'f8', 'f8', 'i4')})
        star_phys = np.genfromtxt(lines, names=True, dtype=('i4', 'i4', 'f8', 'f8', 'i4')) #faster than the above
    stop = time.time()
    print 'Time to load: %s s' % (stop - start)

    return star_phys

def load_planets_stars_phys_separate(file_name_path, run_number):
    #This function loads the simulated physical planets and stars from the individual files their properties were saved in
    
    start = time.time()

    clusterids_per_sys = [] #list to be filled with lists of all the cluster id's per system
    with open(file_name_path + 'clusterids_all%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2].split(', ')
                clusterids_sys = [int(i) for i in line]
                clusterids_per_sys.append(clusterids_sys)
    
    P_per_sys = [] #list to be filled with lists of all periods per system (days)
    N_sys_with_planets = 0 #counter for the number of simulated systems with planets
    with open(file_name_path + 'periods_all%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                N_sys_with_planets += 1
                line = line[1:-2].split(', ')
                P_sys = [float(i) for i in line]
                P_per_sys.append(P_sys)
    #print P_sys

    radii_per_sys = [] #list to be filled with lists of all planet radii per system (solar radii)
    with open(file_name_path + 'radii_all%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2].split(', ')
                radii_sys = [float(i) for i in line]
                radii_per_sys.append(radii_sys)
    #print radii_sys

    mass_per_sys = [] #list to be filled with lists of all planet radii per system (solar masses)
    with open(file_name_path + 'masses_all%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2].split(', ')
                mass_sys = [float(i) for i in line]
                mass_per_sys.append(mass_sys)
    #print mass_sys
    
    e_per_sys = [] #list to be filled with lists of all eccentricities per system
    with open(file_name_path + 'eccentricities_all%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2].split(', ')
                e_sys = [float(i) for i in line]
                e_per_sys.append(e_sys)
    #print e_sys

    Mstar_all = np.loadtxt(file_name_path + 'stellar_masses_with_planets%s.out' % run_number) #array of stellar masses of all the systems with a planetary system, in solar masses
    Rstar_all = np.loadtxt(file_name_path + 'stellar_radii_with_planets%s.out' % run_number) #array of stellar radii of all the systems with a planetary system, in solar radii

    stop = time.time()
    print 'Time to load: %s s' % (stop - start)

    return clusterids_per_sys, P_per_sys, radii_per_sys, mass_per_sys, e_per_sys, Mstar_all, Rstar_all

def compute_summary_stats_from_cat_phys(cat_phys=None, star_phys=None, file_name_path=None, run_number=''):
    #This function takes in a simulated observed catalog of planets 'cat_phys' in table format and returns many arrays (1D and 2D) of the summary stats

    if cat_phys is not None and star_phys is not None:
        i_sys = np.unique(cat_phys['target_id'])
        N_sys_with_planets = len(i_sys) #number of simulated systems with planets

        P_per_sys = [] #list to be filled with lists of all periods per system (days)
        e_per_sys = [] #list to be filled with lists of all eccentricities per system
        radii_per_sys = [] #list to be filled with lists of all planet radii per system (solar radii)
        mass_per_sys = [] #list to be filled with lists of all planet radii per system (solar masses)
        Mstar_all = [] #list to be filled with the stellar masses of all the systems with a planetary system (solar masses)
        Rstar_all = [] #list to be filled with the stellar radii of all the systems with a planetary system (solar radii)
        num_planets_cumu = np.concatenate((np.array([0]), np.cumsum(star_phys['num_planets'])))
        for i in range(len(num_planets_cumu) - 1):
            cat_phys_sys = cat_phys[num_planets_cumu[i]:num_planets_cumu[i+1]]
            
            P_per_sys.append(cat_phys_sys['period'])
            e_per_sys.append(cat_phys_sys['ecc'])
            radii_per_sys.append(cat_phys_sys['planet_radius'])
            mass_per_sys.append(cat_phys_sys['planet_mass'])
            Mstar_all.append(star_phys['star_mass'][i])
            Rstar_all.append(star_phys['star_radius'][i])
        Mstar_all = np.array(Mstar_all)
        Rstar_all = np.array(Rstar_all)

    if file_name_path != None:
        clusterids_per_sys, P_per_sys, radii_per_sys, mass_per_sys, e_per_sys, Mstar_all, Rstar_all = load_planets_stars_phys_separate(file_name_path, run_number)



    clusterids_all = [] #list to be zero-padded so each list of cluster id's is sorted and has the same length, and then converted to an array
    P_all = [] #list to be zero-padded so each list of periods is sorted and has the same length, and then converted to an array
    e_all = [] #list to be zero-padded so each list of eccentricities is sorted (by period) and has the same length, and then converted to an array
    radii_all = [] #list to be zero-padded so each list of radii is sorted (by period) and has the same length, and then converted to an array
    mass_all = [] #list to be zero-padded so each list of masses is sorted (by period) and has the same length, and then converted to an array

    Pmin = 0. #set a minimum period (days), discarding planets less than this period

    Mmax = np.max([len(x) for x in P_per_sys]) #maximum planet multiplicity generated by the clustering method
    clustertot_all = []
    pl_per_cluster_all = []
    for i in range(len(P_per_sys)):
        clusterids_unique = list(set(clusterids_per_sys[i]))
        clustertot_all.append(len(clusterids_unique))
        for c in clusterids_unique:
            pl_per_cluster_all.append(np.sum(np.array(clusterids_per_sys[i]) == c))
        
        i_sorted = np.argsort(P_per_sys[i]) #array of indices which would sort the system by period
        P_sorted = np.array(P_per_sys[i])[i_sorted]
        P_sorted_cut = P_sorted[P_sorted > Pmin]
        clusterids_sorted_cut = np.array(clusterids_per_sys[i])[i_sorted][P_sorted > Pmin]
        e_sorted_cut = np.array(e_per_sys[i])[i_sorted][P_sorted > Pmin]
        radii_sorted_cut = np.array(radii_per_sys[i])[i_sorted][P_sorted > Pmin]
        mass_sorted_cut = np.array(mass_per_sys[i])[i_sorted][P_sorted > Pmin]
        
        P_sys = list(P_sorted_cut) + [0]*(Mmax - len(P_sorted_cut)) #zero-pad the list up to Mmax elements
        clusterids_sys = list(clusterids_sorted_cut) + [0]*(Mmax - len(clusterids_sorted_cut)) #zero-pad the list up to Mmax elements
        e_sys = list(e_sorted_cut) + [0]*(Mmax - len(e_sorted_cut)) #zero-pad the list up to Mmax elements
        radii_sys = list(radii_sorted_cut) + [0]*(Mmax - len(radii_sorted_cut)) #zero-pad the list up to Mmax elements
        mass_sys = list(mass_sorted_cut) + [0]*(Mmax - len(mass_sorted_cut)) #zero-pad the list up to Mmax elements
        
        P_all.append(P_sys)
        clusterids_all.append(clusterids_sys)
        e_all.append(e_sys)
        radii_all.append(radii_sys)
        mass_all.append(mass_sys)
    P_all = np.array(P_all)
    clusterids_all = np.array(clusterids_all)
    e_all = np.array(e_all)
    radii_all = np.array(radii_all)
    mass_all = np.array(mass_all)

    clustertot_all = np.array(clustertot_all)
    pl_per_cluster_all = np.array(pl_per_cluster_all)

    #To convert the radii and masses to Earth units:
    radii_all = radii_all*(Rsun/Rearth) #radii in Earth radii
    mass_all = mass_all*(Msun/Mearth) #masses in Earth masses

    Mtot_all = np.sum(P_all > 0, axis=1) #array of all planet multiplicites
    #Mtot_all = Mtot_all[Mtot_all > 0]



    #To calculate the underlying period ratios, radii ratios, and separations in mutual Hill radii:
    Rm_all = [] #list to be filled with all the period ratios
    radii_ratio_all = [] #list to be filled with all the radii ratios
    N_mH_all = [] #list to be filled with all the separations between adjacent planet pairs in units of mutual Hill radii

    radii_above_all_flat = [] #list to be filled with the radii of planets above the photoevaporation boundary
    radii_below_all_flat = [] #list to be filled with the radii of planets below the photoevaporation boundary
    radii_ratio_above_all_flat = [] #list to be filled with the radii ratios of adjacent planet pairs above the photoevaporation boundary
    radii_ratio_below_all_flat = [] #list to be filled with the radii ratios of adjacent planet pairs below the photoevaporation boundary
    radii_ratio_across_all_flat = [] #list to be filled with the radii ratios of adjacent planet pairs across the photoevaporation boundary

    for i in range(len(P_all)):
        Mstar_system = Mstar_all[i] #mass of the star for this system, in solar masses
        P_all_system = P_all[i][P_all[i] > 0]
        e_all_system = e_all[i][P_all[i] > 0]
        radii_all_system = radii_all[i][P_all[i] > 0]
        mass_all_system = mass_all[i][P_all[i] > 0]
        
        #To calculate all the period ratios:
        Rm_all_system = list(P_all_system[1:]/P_all_system[0:-1]) #list of period ratios in this system
        Rm_all_system = np.array(Rm_all_system + [0]*(Mmax - 1 - len(Rm_all_system))) #to add filler 0's to Rm_all_system to pad it to Mmax - 1 elements
        Rm_all.append(Rm_all_system)
        
        #To calculate all the radii ratios:
        radii_ratio_all_system = list(radii_all_system[1:]/radii_all_system[0:-1]) #list of radii ratios in this system
        radii_ratio_all_system = np.array(radii_ratio_all_system + [0]*(Mmax - 1 - len(radii_ratio_all_system))) #to add filler 0's to radii_ratio_all_system to pad it to Mmax - 1 elements
        radii_ratio_all.append(radii_ratio_all_system)
        
        #To calculate all the separations in mutual Hill radii between adjacent planet pairs:
        a_all_system = a_from_P(P_all_system, Mstar_system)
        R_mH_all_system = ((a_all_system[0:-1] + a_all_system[1:])/2.)*(Mearth*(mass_all_system[0:-1] + mass_all_system[1:])/(3.*Mstar_system*Msun))**(1./3.) #mutual Hill radii between adjacent planet pairs in this system, in AU
        #R_sep_all_system = a_all_system[1:] - a_all_system[0:-1] #separations between adjacent planet pairs in this system, in AU, ignoring eccentricities
        R_sep_all_system = a_all_system[1:]*(1. - e_all_system[1:]) - a_all_system[0:-1]*(1. + e_all_system[0:-1]) #separations between adjacent planet pairs in this system, in AU, including eccentricities
        N_mH_all_system = list(R_sep_all_system/R_mH_all_system) #separations between adjacent planet pairs in this system, in mutual Hill radii
        N_mH_all_system = np.array(N_mH_all_system + [0]*(Mmax - 1 - len(N_mH_all_system))) #to add filler 0's to N_mH_all_system to pad it to Mmax - 1 elements
        N_mH_all.append(N_mH_all_system)
    
        #To separate the planets in the system as above and below the boundary:
        system_above_bools = np.array([photoevap_boundary_Carrera2018(radii_all_system[x], P_all_system[x]) for x in range(len(P_all_system))])
        #if len(system_above_bools) > 1:
        #print system_above_bools
        
        #To record the transit depths of the planets above and below the boundary:
        for j in range(len(radii_all_system)):
            radii_above_all_flat.append(radii_all_system[j]) if system_above_bools[j] == 1 else radii_below_all_flat.append(radii_all_system[j])

        #To record the transit depth ratios of the planets above, below, and across the boundary:
        radii_ratio_all_system = list(radii_all_system[1:]/radii_all_system[0:-1]) #list of radii ratios in this system
        for j in range(len(radii_ratio_all_system)):
            if system_above_bools[j] + system_above_bools[j+1] == 2: #both planets are above the boundary
                radii_ratio_above_all_flat.append(radii_ratio_all_system[j])
            elif system_above_bools[j] + system_above_bools[j+1] == 1: #one planet is above, the other is below the boundary
                radii_ratio_across_all_flat.append(radii_ratio_all_system[j])
            elif system_above_bools[j] + system_above_bools[j+1] == 0: #both planets are below the boundary
                radii_ratio_below_all_flat.append(radii_ratio_all_system[j])

    Rm_all = np.array(Rm_all)
    radii_ratio_all = np.array(radii_ratio_all)
    N_mH_all = np.array(N_mH_all)

    P_all_flat = P_all.flatten() #all the periods of all the planets
    P_all_flat = P_all_flat[P_all_flat > 0]

    Rm_all_flat = Rm_all.flatten() #all the period ratios of all the adjacent planets
    Rm_all_flat = Rm_all_flat[Rm_all_flat > 0]

    N_mH_all_flat = N_mH_all.flatten() #all the separations of all the adjacent planets in units of their mutual Hill radii
    N_mH_all_flat = N_mH_all_flat[N_mH_all_flat > 0]

    #e_all_flat = e_all.flatten() #all the eccentricities of all the planets
    #e_all_flat = e_all_flat[e_all_flat > 0]
    e_all_flat = e_all[P_all > 0] #all the eccentricities of all the planets (which can be zero)

    radii_all_flat = radii_all.flatten() #all the planet radii, in Earth radii
    radii_all_flat = radii_all_flat[radii_all_flat > 0]

    radii_ratio_all_flat = radii_ratio_all.flatten() #all the radii ratios of all the adjacent planets
    radii_ratio_all_flat = radii_ratio_all_flat[radii_ratio_all_flat > 0]

    mass_all_flat = mass_all.flatten() #all the planet masses, in Earth masses
    mass_all_flat = mass_all_flat[mass_all_flat > 0]

    radii_above_all_flat = np.array(radii_above_all_flat)
    radii_below_all_flat = np.array(radii_below_all_flat)
    radii_ratio_above_all_flat = np.array(radii_ratio_above_all_flat)
    radii_ratio_below_all_flat = np.array(radii_ratio_below_all_flat)
    radii_ratio_across_all_flat = np.array(radii_ratio_across_all_flat)

    #'sssp' stands for 'summary stats simulated physical'
    sssp_per_sys = {'clusterids_all': clusterids_all, 'P_all': P_all, 'radii_all': radii_all, 'mass_all': mass_all, 'e_all': e_all, 'Rm_all': Rm_all, 'radii_ratio_all': radii_ratio_all, 'N_mH_all': N_mH_all}
    sssp = {'Mstar_all': Mstar_all, 'Rstar_all': Rstar_all, 'Mtot_all': Mtot_all, 'clustertot_all': clustertot_all, 'pl_per_cluster_all': pl_per_cluster_all, 'P_all': P_all_flat, 'radii_all': radii_all_flat, 'mass_all': mass_all_flat, 'e_all': e_all_flat, 'Rm_all': Rm_all_flat, 'radii_ratio_all': radii_ratio_all_flat, 'N_mH_all': N_mH_all_flat, 'radii_above_all': radii_above_all_flat, 'radii_below_all': radii_below_all_flat, 'radii_ratio_above_all': radii_ratio_above_all_flat, 'radii_ratio_below_all': radii_ratio_below_all_flat, 'radii_ratio_across_all': radii_ratio_across_all_flat}
    return [sssp_per_sys, sssp]





# Functions to load and analyze simulated observed catalogs:

def read_targets_period_radius_bounds(file_name):
    #This function reads the number of simulated targets and bounds for the periods and radii
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'observed_catalog_planets%s.txt' % run_number, 'r')
        for line in file:
            if line[:26] == '# num_targets_sim_pass_one':
                N_sim = int(line[28:])
            elif line[:14] == '# max_incl_sys':
                max_incl_sys = float(line[16:])
                cos_factor = np.cos(max_incl_sys*np.pi/180.)
            elif line[:12] == '# min_period':
                P_min = float(line[14:])
            elif line[:12] == '# max_period':
                P_max = float(line[14:])
            elif line[:12] == '# min_radius':
                radii_min = float(line[24:])
            elif line[:12] == '# max_radius':
                radii_max = float(line[24:])

    return N_sim, cos_factor, P_min, P_max, radii_min, radii_max

def read_sim_params(file_name): #####
    #This function reads the simulation parameters from the file
    param_vals_all = [] #list to be filled with the values of all the model parameters
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'observed_catalog_planets%s.txt' % run_number, 'r')
        for line in file:
            for i in range(len(param_keys_all)):
                chars = len(param_keys_all[i][0])
                if line[:3+chars] == '# ' + param_keys_all[i][0] + ':':
                    if param_keys_all[i][0][:3] == 'log':
                        param_vals_all.append(np.round(np.exp(float(line[4+chars:])), 4))
                    elif (param_keys_all[i][0][:11] == 'num_targets') or (param_keys_all[i][0][:11] == 'mr_max_mass'):
                        param_vals_all.append(int(float(line[4+chars:])))
                    else:
                        param_vals_all.append(np.round(float(line[4+chars:]), 4))

    if len(param_vals_all) != len(param_keys_all):
        print 'Problem with reading parameter values...'

    return param_vals_all

def load_cat_obs(file_name):
    #This function loads the simulated observed catalog of planets
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'observed_catalog_planets%s.txt' % run_number, 'r')
        lines = (line for line in file if not line.startswith('#'))
        cat_obs = np.loadtxt(lines, skiprows=1, dtype={'names': ('target_id', 'star_id', 'period', 'depth', 'duration'), 'formats': ('i4', 'i4', 'f8', 'f8', 'f8')})

    return cat_obs

def load_star_obs(file_name):
    #This function loads the catalog of stars with simulated observed planets
    with open(file_name, 'r') as file: #open(loadfiles_directory + 'observed_catalog_stars%s.txt' % run_number, 'r')
        lines = (line for line in file if not line.startswith('#'))
        star_obs = np.loadtxt(lines, skiprows=1, dtype={'names': ('target_id', 'star_id', 'star_mass', 'star_radius', 'num_obs_planets'), 'formats': ('i4', 'i4', 'f8', 'f8', 'i4')})

    return star_obs

def load_planets_stars_obs_separate(file_name_path, run_number):
    #This function loads the simulated observed planets and stars from the individual files their properties were saved in
    
    P_per_sys = [] #list to be filled with lists of the observed periods per system (days)
    with open(file_name_path + 'periods%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2]
                line_per_sys = line.split('; ')
                #print len(line_per_sys)
                for x in line_per_sys:
                    P_sys = x.split()
                    P_sys = [float(i) for i in P_sys]
                    P_per_sys.append(P_sys)

    D_per_sys = [] #list to be filled with lists of the transit depths per system
    with open(file_name_path + 'depths%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2]
                line_per_sys = line.split('; ')
                #print len(line_per_sys)
                for x in line_per_sys:
                    D_sys = x.split()
                    D_sys = [float(i) for i in D_sys]
                    D_per_sys.append(D_sys)

    tdur_per_sys = [] #list to be filled with lists of the transit durations per system (days)
    with open(file_name_path + 'durations%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2]
                line_per_sys = line.split('; ')
                #print len(line_per_sys)
                for x in line_per_sys:
                    tdur_sys = x.split()
                    tdur_sys = [float(i) for i in tdur_sys]
                    tdur_per_sys.append(tdur_sys)

    Mstar_obs = [] #list to be filled with the stellar masses of the systems with observed planets (Msun)
    with open(file_name_path + 'stellar_masses_obs%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2]
                Mstars = line.split(', ')
                Mstars = [float(i) for i in Mstars]
                Mstar_obs = Mstar_obs + Mstars
    Mstar_obs = np.array(Mstar_obs)
    
    Rstar_obs = [] #list to be filled with the stellar radii of the systems with observed planets (Rsun)
    with open(file_name_path + 'stellar_radii_obs%s.out' % run_number, 'r') as file:
        for line in file:
            if line[0] != '#':
                line = line[1:-2]
                Rstars = line.split(', ')
                Rstars = [float(i) for i in Rstars]
                Rstar_obs = Rstar_obs + Rstars
    Rstar_obs = np.array(Rstar_obs)

    return P_per_sys, D_per_sys, tdur_per_sys, Mstar_obs, Rstar_obs

def compute_summary_stats_from_cat_obs(cat_obs=None, star_obs=None, file_name_path=None, run_number=''):
    #This function takes in a simulated observed catalog of planets 'cat_obs' in table format and returns many arrays (1D and 2D) of the summary stats
    
    if cat_obs is not None and star_obs is not None:
        i_sys = np.unique(cat_obs['target_id'])

        P_per_sys = [] #list to be filled with lists of the observed periods per system (days)
        D_per_sys = [] #list to be filled with lists of the transit depths per system
        tdur_per_sys = [] #list to be filled with lists of the transit durations per system (days)
        Mstar_obs = [] #list to be filled with the stellar masses of the systems with observed planets (Msun)
        Rstar_obs = [] #list to be filled with the stellar radii of the systems with observed planets (Rsun)
        for i in i_sys:
            P_sys = cat_obs['period'][cat_obs['target_id'] == i]
            D_sys = cat_obs['depth'][cat_obs['target_id'] == i]
            tdur_sys = cat_obs['duration'][cat_obs['target_id'] == i]
            Mstar_sys = star_obs['star_mass'][star_obs['target_id'] == i][0]
            Rstar_sys = star_obs['star_radius'][star_obs['target_id'] == i][0]
            
            P_per_sys.append(P_sys)
            D_per_sys.append(D_sys)
            tdur_per_sys.append(tdur_sys)
            Mstar_obs.append(Mstar_sys)
            Rstar_obs.append(Rstar_sys)
        Mstar_obs = np.array(Mstar_obs)
        Rstar_obs = np.array(Rstar_obs)

    if file_name_path != None:
        P_per_sys, D_per_sys, tdur_per_sys, Mstar_obs, Rstar_obs = load_planets_stars_obs_separate(file_name_path, run_number)



    P_obs = [] #list to be zero-padded so each list of periods is sorted and has the same length, and then converted to an array
    D_obs = [] #list to be zero-padded so each list of depths is sorted (by period) and has the same length, and then converted to an array
    tdur_obs = [] #list to be zero-padded so each list of transit durations is sorted (by period) and has the same length, and then converted to an array

    Pmin = 0. #set a minimum period (days), discarding planets less than this period

    Mmax = np.max([len(x) for x in P_per_sys]) #maximum planet multiplicity generated by the clustering method
    for i in range(len(P_per_sys)):
        i_sorted = np.argsort(P_per_sys[i]) #array of indices which would sort the system by period
        P_sorted = np.array(P_per_sys[i])[i_sorted]
        P_sorted_cut = P_sorted[P_sorted > Pmin]
        D_sorted_cut = np.array(D_per_sys[i])[i_sorted][P_sorted > Pmin]
        tdur_sorted_cut = np.array(tdur_per_sys[i])[i_sorted][P_sorted > Pmin]
        
        P_sys = list(P_sorted_cut) + [-1]*(Mmax - len(P_sorted_cut)) #zero-pad the list up to Mmax elements
        D_sys = list(D_sorted_cut) + [-1]*(Mmax - len(D_sorted_cut)) #zero-pad the list up to Mmax elements
        tdur_sys = list(tdur_sorted_cut) + [-1]*(Mmax - len(tdur_sorted_cut)) #zero-pad the list up to Mmax elements

        P_obs.append(P_sys)
        D_obs.append(D_sys)
        tdur_obs.append(tdur_sys)
    P_obs = np.array(P_obs)
    D_obs = np.array(D_obs)
    tdur_obs = np.array(tdur_obs)*24.*60. #tdur_obs converted to mins

    Mtot_obs = np.sum(P_obs > 0, axis=1) #array of observed planet multiplicites
    Nmult_obs = np.array([np.sum(Mtot_obs == x) for x in range(1,Mmax+1)]) #array of total numbers of systems with observed planet multiplicities of 1,2,3,...,Mmax planets
    radii_obs = np.sqrt(D_obs)*np.transpose([Rstar_obs])*(Rsun/Rearth) #array of planet radii, in Earth radii



    #To calculate the observed period ratios, period-normalized transit duration ratios, and transit depth ratios:
    Rm_obs = [] #list to be filled with the observed period ratios
    D_ratio_obs = [] #list to be filled with the observed transit depth ratios
    xi_obs = [] #list to be filled with the period-normalized transit duration ratios
    xi_res_obs = [] #list to be filled with the period-normalized transit duration ratios for planet pairs near resonance
    xi_res32_obs = []
    xi_res21_obs = []
    xi_nonres_obs = [] #list to be filled with the period_normalized transit duration ratios for planet pairs not in resonance

    D_above_obs_flat = [] #list to be filled with the transit depths of observed planets above the photoevaporation boundary
    D_below_obs_flat = [] #list to be filled with the transit depths of observed planets below the photoevaporation boundary
    D_ratio_above_obs_flat = [] #list to be filled with the transit depth ratios of adjacent observed planet pairs above the photoevaporation boundary
    D_ratio_below_obs_flat = [] #list to be filled with the transit depth ratios of adjacent observed planet pairs below the photoevaporation boundary
    D_ratio_across_obs_flat = [] #list to be filled with the transit depth ratios of adjacent observed planet pairs across the photoevaporation boundary

    for i in range(len(P_obs)):
        P_obs_system = P_obs[i][P_obs[i] > 0]
        radii_obs_system = radii_obs[i][P_obs[i] > 0]
        tdur_obs_system = tdur_obs[i][P_obs[i] > 0]
        D_obs_system = D_obs[i][P_obs[i] > 0]
        
        #To calculate all the observed period ratios:
        Rm_obs_system = list(P_obs_system[1:]/P_obs_system[0:-1]) #list of period ratios observed in this system
        Rm_obs_system = np.array(Rm_obs_system + [-1]*(Mmax - 1 - len(Rm_obs_system))) #to add filler 0's to Rm_obs_system to pad it to Mmax - 1 elements
        Rm_obs.append(Rm_obs_system)
        
        #To calculate all the observed transit depth ratios:
        D_ratio_obs_system = list(D_obs_system[1:]/D_obs_system[0:-1]) #list of transit depth ratios observed in this system
        D_ratio_obs_system = np.array(D_ratio_obs_system + [-1]*(Mmax - 1 - len(D_ratio_obs_system))) #to add filler 0's to D_ratio_obs_system to pad it to Mmax - 1 elements
        D_ratio_obs.append(D_ratio_obs_system)
        
        #To calculate all the period-normalized transit duration ratios:
        xi_obs_system = list((tdur_obs_system[0:-1]/tdur_obs_system[1:])*(P_obs_system[1:]/P_obs_system[0:-1])**(1./3.)) #list of period-normalized transit duration ratios in this system
        xi_obs_system = np.array(xi_obs_system + [-1]*(Mmax - 1 - len(xi_obs_system))) #to add filler 0's to xi_obs_system to pad it to Mmax - 1 elements
        xi_obs.append(xi_obs_system)
        
        #To separate the period-normalized transit duration ratios for planet pairs near vs. not in resonance:
        mask_res_system = np.zeros(len(Rm_obs_system), dtype=bool)
        mask_res32_system = np.zeros(len(Rm_obs_system), dtype=bool)
        mask_res21_system = np.zeros(len(Rm_obs_system), dtype=bool)
        
        for ratio in res_ratios:
            mask_res_system[(Rm_obs_system >= ratio) & (Rm_obs_system <= ratio*(1.+res_width))] = 1

        mask_res32_system[(Rm_obs_system >= 1.5) & (Rm_obs_system <= 1.5*(1.+res_width))] = 1
        mask_res21_system[(Rm_obs_system >= 2.) & (Rm_obs_system <= 2.*(1.+res_width))] = 1
        
        xi_res_obs_system = list(xi_obs_system[mask_res_system])
        xi_res32_obs_system = list(xi_obs_system[mask_res32_system])
        xi_res21_obs_system = list(xi_obs_system[mask_res21_system])
        xi_nonres_obs_system = list(xi_obs_system[~mask_res_system])
        xi_res_obs_system = np.array(xi_res_obs_system + [-1]*(10 - len(xi_res_obs_system)))
        xi_res32_obs_system = np.array(xi_res32_obs_system + [-1]*(10 - len(xi_res32_obs_system)))
        xi_res21_obs_system = np.array(xi_res21_obs_system + [-1]*(10 - len(xi_res21_obs_system)))
        xi_nonres_obs_system = np.array(xi_nonres_obs_system + [-1]*(10 - len(xi_nonres_obs_system)))
        xi_res_obs.append(xi_res_obs_system)
        xi_res32_obs.append(xi_res32_obs_system)
        xi_res21_obs.append(xi_res21_obs_system)
        xi_nonres_obs.append(xi_nonres_obs_system)

        #To separate the planets in the system as above and below the boundary:
        system_above_bools = np.array([photoevap_boundary_Carrera2018(radii_obs_system[x], P_obs_system[x]) for x in range(len(P_obs_system))])
        #if len(system_above_bools) > 1:
        #print system_above_bools
        
        #To record the transit depths of the planets above and below the boundary:
        for j in range(len(D_obs_system)):
            D_above_obs_flat.append(D_obs_system[j]) if system_above_bools[j] == 1 else D_below_obs_flat.append(D_obs_system[j])
            
        #To record the transit depth ratios of the planets above, below, and across the boundary:
        D_ratio_obs_system = list(D_obs_system[1:]/D_obs_system[0:-1]) #list of transit depth ratios observed in this system
        for j in range(len(D_ratio_obs_system)):
            if system_above_bools[j] + system_above_bools[j+1] == 2: #both planets are above the boundary
                D_ratio_above_obs_flat.append(D_ratio_obs_system[j])
            elif system_above_bools[j] + system_above_bools[j+1] == 1: #one planet is above, the other is below the boundary
                D_ratio_across_obs_flat.append(D_ratio_obs_system[j])
            elif system_above_bools[j] + system_above_bools[j+1] == 0: #both planets are below the boundary
                D_ratio_below_obs_flat.append(D_ratio_obs_system[j])

    Rm_obs = np.array(Rm_obs)
    D_ratio_obs = np.array(D_ratio_obs)
    xi_obs = np.array(xi_obs)
    xi_res_obs = np.array(xi_res_obs)
    xi_res32_obs = np.array(xi_res32_obs)
    xi_res21_obs = np.array(xi_res21_obs)
    xi_nonres_obs = np.array(xi_nonres_obs)

    P_obs_flat = P_obs.flatten() #all the observed periods of all the planets
    P_obs_flat = P_obs_flat[P_obs_flat > 0]

    Rm_obs_flat = Rm_obs.flatten() #all the observed period ratios of all the observed adjacent planets
    Rm_obs_flat = Rm_obs_flat[Rm_obs_flat > 0]

    D_obs_flat = D_obs.flatten() #all the transit depths
    D_obs_flat = D_obs_flat[D_obs_flat > 0]
    radii_obs_flat = radii_obs.flatten() #all the observed planet radii, in Earth radii
    radii_obs_flat = radii_obs_flat[radii_obs_flat > 0]

    D_ratio_obs_flat = D_ratio_obs.flatten() #all the transit depth ratios
    D_ratio_obs_flat = D_ratio_obs_flat[D_ratio_obs_flat > 0]

    tdur_obs_flat = tdur_obs.flatten() #all the observed transit durations, in mins
    tdur_obs_flat = tdur_obs_flat[tdur_obs_flat >= 0] #####

    xi_obs_flat = xi_obs.flatten() #all the observed period-normalized transit duration ratios
    xi_obs_flat = xi_obs_flat[xi_obs_flat > 0]
    xi_res_obs_flat = xi_res_obs.flatten() #the observed period-normalized transit duration ratios for planet pairs near resonance
    xi_res_obs_flat = xi_res_obs_flat[xi_res_obs_flat > 0]
    xi_res32_obs_flat = xi_res32_obs.flatten()
    xi_res32_obs_flat = xi_res32_obs_flat[xi_res32_obs_flat > 0]
    xi_res21_obs_flat = xi_res21_obs.flatten()
    xi_res21_obs_flat = xi_res21_obs_flat[xi_res21_obs_flat > 0]
    xi_nonres_obs_flat = xi_nonres_obs.flatten() #the observed period-normalized transit duration ratios for planet pairs not in resonance
    xi_nonres_obs_flat = xi_nonres_obs_flat[xi_nonres_obs_flat > 0]

    D_above_obs_flat = np.array(D_above_obs_flat)
    D_below_obs_flat = np.array(D_below_obs_flat)
    D_ratio_above_obs_flat = np.array(D_ratio_above_obs_flat)
    D_ratio_below_obs_flat = np.array(D_ratio_below_obs_flat)
    D_ratio_across_obs_flat = np.array(D_ratio_across_obs_flat)

    #'sss' stands for 'summary stats simulated'
    sss_per_sys = {'P_obs': P_obs, 'D_obs': D_obs, 'tdur_obs': tdur_obs, 'radii_obs': radii_obs, 'Rm_obs': Rm_obs, 'D_ratio_obs': D_ratio_obs, 'xi_obs': xi_obs, 'xi_res_obs': xi_res_obs, 'xi_res32_obs': xi_res32_obs, 'xi_res21_obs': xi_res21_obs, 'xi_nonres_obs': xi_nonres_obs}
    sss = {'Mstar_obs': Mstar_obs, 'Rstar_obs': Rstar_obs, 'Mtot_obs': Mtot_obs, 'Nmult_obs': Nmult_obs, 'P_obs': P_obs_flat, 'D_obs': D_obs_flat, 'tdur_obs': tdur_obs_flat, 'radii_obs': radii_obs_flat, 'Rm_obs': Rm_obs_flat, 'D_ratio_obs': D_ratio_obs_flat, 'xi_obs': xi_obs_flat, 'xi_res_obs': xi_res_obs_flat, 'xi_res32_obs': xi_res32_obs_flat, 'xi_res21_obs': xi_res21_obs_flat, 'xi_nonres_obs': xi_nonres_obs_flat, 'D_above_obs': D_above_obs_flat, 'D_below_obs': D_below_obs_flat, 'D_ratio_above_obs': D_ratio_above_obs_flat, 'D_ratio_below_obs': D_ratio_below_obs_flat, 'D_ratio_across_obs': D_ratio_across_obs_flat}
    return [sss_per_sys, sss]





# Functions to load and analyze the Kepler observed catalog:

def load_Kepler_planets_cleaned():
    
    planets_cleaned = np.genfromtxt('ExoplanetsSysSim_Clusters/SysSimExClusters/plotting/q1_q17_dr25_gaia_fgk_koi_cleaned.csv', dtype={'names': ('kepid', 'KOI', 'koi_disposition', 'koi_pdisposition', 'koi_score', 'P', 't_D', 'depth', 'Rp', 'Rstar'), 'formats': ('i8', 'S9', 'S15', 'S15', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',)}, delimiter=',') #orbit periods 'P' are in days; transit durations 't_D' are in hrs; transit depths 'depth' are in ppm; planetary radii 'Rp' are in Rearth; stellar radii 'Rstar' are in Rsolar
    planets_cleaned = planets_cleaned[1:]
    return planets_cleaned

def load_Kepler_stars_cleaned():

    stars_cleaned = np.genfromtxt('ExoplanetsSysSim_Clusters/SysSimExClusters/plotting/q1_q17_dr25_gaia_fgk_cleaned.csv', dtype={'names': ('kepid', 'mass', 'radius'), 'formats': ('i8', 'f8', 'f8')}, delimiter=',') #orbit periods 'P' are in days; transit durations 't_D' are in hrs; transit depths 'depth' are in ppm; planetary radii 'Rp' are in Rearth; stellar radii 'Rstar' are in Rsolar
    stars_cleaned = stars_cleaned[1:]
    return stars_cleaned

def compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max):

    planets_cleaned = load_Kepler_planets_cleaned()
    
    '''
    #To load and compute the exoplanet multiplicities, periods, and period ratios of the confirmed Kepler exoplanets:
    Q1Q17_DR25 = np.genfromtxt('q1_q17_dr25_koi.tab_selectcols_new.csv', dtype={'names': ('KepID', 'KOI', 'Archive_Disp', 'Kepler_Disp', 'Disp', 'P', 't_D', 'depth', 'Rp', 'Rstar'), 'formats': ('i8', 'S9', 'S15', 'S15', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',)}, delimiter=',', usecols=(1,2,3,4,5,10,11,12,13,14)) #orbit periods 'P' are in days; transit durations 't_D' are in hrs; transit depths 'depth' are in ppm; planetary radii 'Rp' are in Rearth; stellar radii 'Rstar' are in Rsolar
    Q1Q17_DR25 = Q1Q17_DR25[1:] #skip_header doesn't work so manually get rid of first row of NaNs

    #If using a stellar table with cuts already made and just matching kepids to get a koi catalog:
    in_stellar_catalog_kepois = np.loadtxt('ExoplanetsSysSim_Clusters/clusters_v0.7/kepoi_names.txt', delimiter=', ', dtype='S10')
    planets_cleaned = Q1Q17_DR25[(Q1Q17_DR25['Archive_Disp'] == 'CONFIRMED') + (Q1Q17_DR25['Archive_Disp'] == 'CANDIDATE')]

    in_stellar_catalog_indices = []
    for kepoi in in_stellar_catalog_kepois:
        if len(np.where(planets_cleaned['KOI'] == kepoi)[0]) == 1:
            in_stellar_catalog_indices.append(np.where(planets_cleaned['KOI'] == kepoi)[0][0])
        elif len(np.where(planets_cleaned['KOI'] == kepoi)[0]) == 0:
            continue
        else:
            print 'More than one match to the kepoi in the Exoplanets table!'
    print 'Number of CONFIRMED and CANDIDATE planets in sample: ', len(in_stellar_catalog_indices)

    planets_cleaned = planets_cleaned[in_stellar_catalog_indices]

    #The stellar properties (and thus also planet radii) in 'planets_cleaned' (i.e. the file loaded into 'Q1Q17_DR25') are not as reliable as the stellar properties in the Gaia catalog; will replace them in the following lines:
    kepid_Rstar_gaia_all = np.loadtxt('ExoplanetsSysSim_Clusters/clusters_v0.7/stellar_radii_q1_q17_dr25_gaia_fgk.txt', delimiter=',', skiprows=1, dtype={'names': ('KepID', 'Rstar'), 'formats': ('i8', 'f8',)})

    for i,kepid in enumerate(planets_cleaned['KepID']):
        stellar_radii_new = kepid_Rstar_gaia_all['Rstar'][np.where(kepid_Rstar_gaia_all['KepID'] == kepid)[0][0]]
        planets_cleaned['Rstar'][i] = stellar_radii_new
        planets_cleaned['Rp'][i] = (Rsun/Rearth)*stellar_radii_new*np.sqrt(planets_cleaned['depth'][i]/(1e6))
    '''
    
    #To make cuts in period and planetary radii:
    planets_cleaned = planets_cleaned[(planets_cleaned['P'] > P_min) & (planets_cleaned['P'] < P_max)]
    planets_cleaned = planets_cleaned[(planets_cleaned['Rp'] > radii_min) & (planets_cleaned['Rp'] < radii_max)]

    Rstar_confirmed = planets_cleaned['Rstar']
    


    #To compute the arrays of observables:
    KOI_systems = np.array([x[:6] for x in planets_cleaned['KOI']])
    checked_bools = np.zeros(len(planets_cleaned)) #0's denote KOI that were not checked yet; 1's denote already checked KOI
    M_confirmed = [] #list to be filled with the planet multiplicities of systems with confirmed planets
    R_confirmed = [] #list to be filled with period ratios of adjacent confirmed planet pairs

    P_sys_confirmed = [] #list to be filled with lists of planet periods per system
    radii_sys_confirmed = [] #list to be filled with lists of planet radii per system

    D_ratio_confirmed = [] #list to be filled with the transit depth ratios of adjacent confirmed planet pairs
    xi_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs
    xi_res_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near resonance
    xi_res32_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near 3:2 resonance
    xi_res21_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near 2:1 resonance
    xi_nonres_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs not in resonance
    t_D_confirmed = planets_cleaned['t_D'] #array of the transit durations (hrs) of all the confirmed planets
    D_confirmed = planets_cleaned['depth']/(1e6) #array of the transit depths (fraction) of all the confirmed planets
    radii_confirmed = planets_cleaned['Rp'] #array of the planetary radii (Rearth) of all the confirmed planets

    D_above_confirmed = [] #list to be filled with the transit depths of confirmed planets above the photoevaporation boundary
    D_below_confirmed = [] #list to be filled with the transit depths of confirmed planets below the photoevaporation boundary
    D_ratio_above_confirmed = [] #list to be filled with the transit depth ratios of adjacent confirmed planet pairs above the photoevaporation boundary
    D_ratio_below_confirmed = [] #list to be filled with the transit depth ratios of adjacent confirmed planet pairs below the photoevaporation boundary
    D_ratio_across_confirmed = [] #list to be filled with the transit depth ratios of adjacent confirmed planet pairs across the photoevaporation boundary

    for i in range(len(KOI_systems)):
        if checked_bools[i] == 0: #if the KOI has not been checked (included while looking at another planet in the same system)
            system_i = np.where(KOI_systems == KOI_systems[i])[0]
            checked_bools[system_i] = 1
            
            #To get the periods and transit durations in this system:
            system_P = planets_cleaned['P'][system_i] #periods of all the planets in this system
            system_t_D = planets_cleaned['t_D'][system_i] #transit durations of all the planets in this system
            system_D = planets_cleaned['depth'][system_i] #transit depths of all the planets in this system
            system_radii = planets_cleaned['Rp'][system_i] #radii of all the planets in this system
            system_sort_i = np.argsort(system_P) #indices that would sort the periods of the planets in this system
            system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
            system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
            system_D = system_D[system_sort_i] #transit depths of all the planets in this system, sorted by period
            system_radii = system_radii[system_sort_i] #radii of all the planets in this system, sorted by period
            
            #To count the total number of planets in this system:
            M_confirmed.append(len(system_P))
            P_sys_confirmed.append(list(system_P) + [0]*(6 - len(system_P)))
            radii_sys_confirmed.append(list(system_radii) + [0]*(6 - len(system_radii)))
            
            #To compute the period ratios and period-normalized transit duration ratios in this system (and separate into planet pairs near vs. not in resonance):
            system_R = system_P[1:]/system_P[0:-1] #period ratios of all the adjacent planet pairs in this system
            system_D_ratio = system_D[1:]/system_D[0:-1] #transit depth ratios of all the adjacent planet pairs in this system
            system_xi = (system_t_D[0:-1]/system_t_D[1:])*(system_P[1:]/system_P[0:-1])**(1./3.) #period-normalized transit duration ratios of all the adjacent planet pairs in this system
            mask_res_system = np.zeros(len(system_R), dtype=bool)
            mask_res32_system = np.zeros(len(system_R), dtype=bool)
            mask_res21_system = np.zeros(len(system_R), dtype=bool)
            
            for ratio in res_ratios:
                mask_res_system[(system_R >= ratio) & (system_R <= ratio*(1.+res_width))] = 1

            mask_res32_system[(system_R >= 1.5) & (system_R <= 1.5*(1.+res_width))] = 1
            mask_res21_system[(system_R >= 2.) & (system_R <= 2.*(1.+res_width))] = 1
            system_xi_res = system_xi[mask_res_system]
            system_xi_res32 = system_xi[mask_res32_system]
            system_xi_res21 = system_xi[mask_res21_system]
            system_xi_nonres = system_xi[~mask_res_system]
            #if sum(mask_res_system) > 0:
            #print system_R[mask_res_system], system_xi_res
            
            for R in system_R:
                R_confirmed.append(R)
            for D_ratio in system_D_ratio:
                D_ratio_confirmed.append(D_ratio)
            for xi in system_xi:
                xi_confirmed.append(xi)
            for xi in system_xi_res:
                xi_res_confirmed.append(xi)
            for xi in system_xi_res32:
                xi_res32_confirmed.append(xi)
            for xi in system_xi_res21:
                xi_res21_confirmed.append(xi)
            for xi in system_xi_nonres:
                xi_nonres_confirmed.append(xi)

            #To separate the planets in the system as above and below the boundary:
            system_above_bools = np.array([photoevap_boundary_Carrera2018(system_radii[x], system_P[x]) for x in range(len(system_P))])
            #if len(system_above_bools) > 1:
            #print system_above_bools
            
            #To record the transit depths of the planets above and below the boundary:
            for j in range(len(system_D)):
                D_above_confirmed.append(system_D[j]) if system_above_bools[j] == 1 else D_below_confirmed.append(system_D[j])
                
            #To record the transit depth ratios of the planets above, below, and across the boundary:
            for j in range(len(system_D_ratio)):
                if system_above_bools[j] + system_above_bools[j+1] == 2: #both planets are above the boundary
                    D_ratio_above_confirmed.append(system_D_ratio[j])
                elif system_above_bools[j] + system_above_bools[j+1] == 1: #one planet is above, the other is below the boundary
                    D_ratio_across_confirmed.append(system_D_ratio[j])
                elif system_above_bools[j] + system_above_bools[j+1] == 0: #both planets are below the boundary
                    D_ratio_below_confirmed.append(system_D_ratio[j])

    P_sys_confirmed = np.array(P_sys_confirmed)
    radii_sys_confirmed = np.array(radii_sys_confirmed)

    P_confirmed = planets_cleaned['P']
    M_confirmed = np.array(M_confirmed)
    Nmult_confirmed = np.array([np.sum(M_confirmed == x) for x in range(1,np.max(M_confirmed)+1)])
    R_confirmed = np.array(R_confirmed)
    D_ratio_confirmed = np.array(D_ratio_confirmed)
    xi_confirmed = np.array(xi_confirmed)
    xi_res_confirmed = np.array(xi_res_confirmed)
    xi_res32_confirmed = np.array(xi_res32_confirmed)
    xi_res21_confirmed = np.array(xi_res21_confirmed)
    xi_nonres_confirmed = np.array(xi_nonres_confirmed)

    D_above_confirmed = np.array(D_above_confirmed)/(1e6)
    D_below_confirmed = np.array(D_below_confirmed)/(1e6)
    D_ratio_above_confirmed = np.array(D_ratio_above_confirmed)
    D_ratio_below_confirmed = np.array(D_ratio_below_confirmed)
    D_ratio_across_confirmed = np.array(D_ratio_across_confirmed)

    #'ssk' stands for 'summary stats Kepler'
    ssk_per_sys = {'P_obs': P_sys_confirmed, 'radii_obs': radii_sys_confirmed}
    ssk = {'Rstar_obs': Rstar_confirmed, 'Mtot_obs': M_confirmed, 'Nmult_obs': Nmult_confirmed, 'P_obs': P_confirmed, 'D_obs': D_confirmed, 'tdur_obs': t_D_confirmed, 'radii_obs': radii_confirmed, 'Rm_obs': R_confirmed, 'D_ratio_obs': D_ratio_confirmed, 'xi_obs': xi_confirmed, 'xi_res_obs': xi_res_confirmed, 'xi_res32_obs': xi_res32_confirmed, 'xi_res21_obs': xi_res21_confirmed, 'xi_nonres_obs': xi_nonres_confirmed, 'D_above_obs': D_above_confirmed, 'D_below_obs': D_below_confirmed, 'D_ratio_above_obs': D_ratio_above_confirmed, 'D_ratio_below_obs': D_ratio_below_confirmed, 'D_ratio_across_obs': D_ratio_across_confirmed}
    return [ssk_per_sys, ssk]





# Functions for computing distances between the simulated and Kepler observed catalogs:

def CRPD_dist(En, On): #NOTE: this distance can potentially give negative values?!
    #This function computes the Cressie Read Power Divergence statistic for observed planet multiplicities
    #En and On must be arrays of the total numbers of systems with 1,2,3,... observed planets, in the simulated (i.e. expected) and the actual (i.e. observed Kepler) data, respectively
    n_max = max(len(En), len(On))
    En = np.array(list(En) + [0]*(n_max - len(En)))
    On = np.array(list(On) + [0]*(n_max - len(On)))
    
    E_array = En/float(np.sum(En)) #normalized numbers (fractions) of simulated systems with 1,2,3,... observed planets
    O_array = On/float(np.sum(On)) #normalized numbers (fractions) of actual Kepler systems with 1,2,3,... observed planets
    rho = 0.
    for i,E_i in enumerate(E_array):
        if En[i] != 0:
            rho += O_array[i]*((O_array[i]/E_array[i])**(2./3.) - 1.)
    rho = (9./5.)*rho

    return rho

def KS_dist_mult(x1, x2):
    #This function computes the K-S distance for two discrete, integer distributions (for multiplicities)
    #This function returns two values: the K-S distance and the x value corresponding to that distance
    x12_max = np.max((np.max(x1), np.max(x2))) #maximum value of x1 and x2
    x1_counts, x1_bins = np.histogram(x1, bins=x12_max, range=(0.5, x12_max+0.5))
    x2_counts, x2_bins = np.histogram(x2, bins=x12_max, range=(0.5, x12_max+0.5))
    pdf_diffs = x1_counts/np.float(len(x1)) - x2_counts/np.float(len(x2))
    cdf_diffs = np.cumsum(pdf_diffs)
    KS_dist = np.max(np.abs(cdf_diffs)) #K-S distance
    KS_x = np.arange(1, x12_max+1)[np.where(np.abs(cdf_diffs) == KS_dist)[0][0]] #x value where the distance is the largest
    
    return KS_dist, KS_x

def KS_dist(x1, x2):
    #This function computes the K-S distance for two continuous distributions (no repeated values)
    #This function returns two values: the K-S distance and the x value corresponding to that distance
    x_all = np.concatenate((x1, x2)) #combined array
    i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
    pdf_diffs = np.concatenate((np.ones(len(x1))/np.float(len(x1)), -np.ones(len(x2))/np.float(len(x2))))[i_all_sorted]
    cdf_diffs = np.cumsum(pdf_diffs)
    KS_dist = np.max(np.abs(cdf_diffs)) #K-S distance
    KS_x = x_all[i_all_sorted][np.where(np.abs(cdf_diffs) == KS_dist)[0][0]] #x value (a value in either x1 or x2) where the distance is the largest
    
    return KS_dist, KS_x

def AD_dist(x1, x2):
    #This function computes and returns the AD distance for two continuous distributions (no repeated values), according to A. N. Pettitt (1976) Eq. (1.2)
    n, m = len(x1), len(x2)
    if n > 1 and m > 1:
        N = n + m
        x_all = np.concatenate((x1, x2)) #combined array
        i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
        M_i_diffs = np.concatenate((np.ones(n), np.zeros(m)))[i_all_sorted]
        M_i_array = np.cumsum(M_i_diffs)[:-1] #array of M_i except for last element, i.e. from i=1 to i=N-1
        i_array = 1. + np.arange(N-1) #array of i from i=1 to i=N-1
        AD_dist = (1./(n*m))*np.sum(((M_i_array*N - n*i_array)**2.)/(i_array*(N - i_array))) #AD distance
    else:
        print 'Not enough points to compute AD distance; returning inf.'
        AD_dist = np.inf
    
    return AD_dist

def AD_dist2(x1, x2): #I tested this and it returns the same results as AD_dist()
    #This function computes and returns the AD distance for two continuous distributions (no repeated values), according to Scholz & Stephens (1987) Eq. (3)
    n1, n2 = len(x1), len(x2)
    if n > 1 and m > 1:
        N = n1 + n2
        x_all = np.concatenate((x1, x2)) #combined array
        i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
        
        M_1j_diffs = np.concatenate((np.ones(n1), np.zeros(n2)))[i_all_sorted]
        M_1j_array = np.cumsum(M_1j_diffs)[:-1] #array of M_1j except for last element, i.e. from j=1 to j=N-1
        M_2j_diffs = np.concatenate((np.zeros(n1), np.ones(n2)))[i_all_sorted]
        M_2j_array = np.cumsum(M_2j_diffs)[:-1] #array of M_2j except for last element, i.e. from j=1 to j=N-1
        j_array = 1. + np.arange(N-1) #array of j from j=1 to j=N-1
        
        AD_dist = (1./N)*((1./n1)*np.sum(((N*M_1j_array - n1*j_array)**2.)/(j_array*(N - j_array))) + (1./n2)*np.sum(((N*M_2j_array - n2*j_array)**2.)/(j_array*(N - j_array)))) #AD distance
    else:
        print 'Not enough points to compute AD distance; returning inf.'
        AD_dist = np.inf
    
    return AD_dist

def AD_mod_dist(x1, x2):
    #This function is the same as 'AD_dist' (or 'AD_dist2') except without the factor of 'nm/N' before the integral
    n, m = len(x1), len(x2)
    if n > 1 and m > 1:
        N = n + m
        x_all = np.concatenate((x1, x2)) #combined array
        i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
        M_i_diffs = np.concatenate((np.ones(n), np.zeros(m)))[i_all_sorted]
        M_i_array = np.cumsum(M_i_diffs)[:-1] #array of M_i except for last element, i.e. from i=1 to i=N-1
        i_array = 1. + np.arange(N-1) #array of i from i=1 to i=N-1
        AD_dist = (N/((n*m)**2.))*np.sum(((M_i_array*N - n*i_array)**2.)/(i_array*(N - i_array))) #AD distance
    else:
        print 'Not enough points to compute AD distance; returning inf.'
        AD_dist = np.inf
    
    return AD_dist

def load_weights(file_name, KS_or_AD):
    #This function loads the weights for the K-S or A-D (or A-D mod) distances
    
    with open(file_name, 'r') as file:
        for line in file:
            if line[0:10] == 'Weights_'+KS_or_AD:
                weights = np.array([float(x) for x in line[13:-2].split(', ')])
    if len(weights) == 17:
        return weights
    else:
        print 'Problem with loading weights from file...'
        return weights

def load_weights_and_compute_weighted_dist(file_name, KS_or_AD, dists,  dists_exclude=[], print_dists=True):
    #This function loads the weights, computes the weighted distances, and prints them
    
    weights = load_weights(file_name, KS_or_AD)
    if len(dists_exclude) > 0:
        weights[dists_exclude] = 0
    dists_weighted = dists*weights
    
    if print_dists == True:
        print '#####'
        print 'Distances for (delta_f, CRPD, CRPD_switched, Multiplicity, P, P ratio, t_dur, xi, xi_nonres, xi_res, depth, depth_above, depth_below, depth_ratio, depth_ratio_above, depth_ratio_below, depth_ratio_across):'
        print 'Distances (%s): ' % KS_or_AD, [float(format(x, '.5f')) for x in dists]
        print 'Weights (%s): ' % KS_or_AD, [float(format(x, '.5f')) for x in weights]
        print 'Weighted distances (%s): ' % KS_or_AD, [float(format(x, '.5f')) for x in dists_weighted]
        print 'Total weighted distance (%s) = %s' % (KS_or_AD, sum(dists_weighted))
        print '#####'
    return sum(dists_weighted)

def compute_distances_sim_Kepler(file_name, cat_obs, star_obs, P_min, P_max, radii_min, radii_max, AD_mod='true', dists_exclude=[], print_dists=True, run_number=''):
    #This function computes the K-S (and their positions), A-D, and other distances as well as additional statistics:
    
    if cat_obs is None and star_obs is None:
        N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(file_name + 'periods%s.out' % run_number)
        sss_per_sys, sss = compute_summary_stats_from_cat_obs(file_name_path=file_name, run_number=run_number)
    else:
        N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(file_name)
        sss_per_sys, sss = compute_summary_stats_from_cat_obs(cat_obs=cat_obs, star_obs=star_obs)
    ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)

    print '(Planets Kepler obs, Planet pairs Kepler obs) = (%s, %s)' % (len(ssk['P_obs']), len(ssk['Rm_obs']))
    print '(Planets obs, Planet pairs obs) = (%s, %s)' % (len(sss['P_obs']), len(sss['Rm_obs']))



    # Misc distances/statistics:
    delta_f = np.abs(len(sss['P_obs'])/(float(N_sim)/cos_factor) - len(ssk['P_obs'])/float(N_Kep)) #absolute difference in the rates of observed planets per star
    
    Nmult_obs_sim_5plus = np.array(list(sss['Nmult_obs'][:4]) + [sum(sss['Nmult_obs'][4:])])
    Nmult_obs_Kep_5plus = np.array(list(ssk['Nmult_obs'][:4]) + [sum(ssk['Nmult_obs'][4:])])
    d_mult_CRPD = CRPD_dist(Nmult_obs_sim_5plus, Nmult_obs_Kep_5plus)
    d_mult_CRPD_switched = CRPD_dist(Nmult_obs_Kep_5plus, Nmult_obs_sim_5plus)
    
    R_res32_sim, R_res32_confirmed = np.float(sum((sss['Rm_obs'] >= 1.5) & (sss['Rm_obs'] <= 1.5*(1.+res_width))))/np.float(len(sss['Rm_obs'])), np.float(sum((ssk['Rm_obs'] >= 1.5) & (ssk['Rm_obs'] <= 1.5*(1.+res_width))))/np.float(len(ssk['Rm_obs'])) #fractions of planet pairs within 5% of 3:2 MMR, for simulated and Kepler data
    R_res21_sim, R_res21_confirmed = np.float(sum((sss['Rm_obs'] >= 2.) & (sss['Rm_obs'] <= 2.*(1.+res_width))))/np.float(len(sss['Rm_obs'])), np.float(sum((ssk['Rm_obs'] >= 2.) & (ssk['Rm_obs'] <= 2.*(1.+res_width))))/np.float(len(ssk['Rm_obs'])) #fractions of planet pairs within 5% of 2:1 MMR, for simulated and Kepler data
    R_res32_diff = np.abs(R_res32_sim - R_res32_confirmed) #difference in fractions of planet pairs close to 3:2 MMR between simulated and Kepler data
    R_res21_diff = np.abs(R_res21_sim - R_res21_confirmed) #difference in fractions of planet pairs close to 2:1 MMR between simulated and Kepler data

    # KS distances and their positions:
    M_KS, M_KS_pos = KS_dist_mult(sss['Mtot_obs'][sss['Mtot_obs'] > 0], ssk['Mtot_obs'])
    P_KS, P_KS_pos = KS_dist(sss['P_obs'], ssk['P_obs'])
    R_KS, R_KS_pos = KS_dist(sss['Rm_obs'], ssk['Rm_obs'])
    tdur_KS, tdur_KS_pos = KS_dist(sss['tdur_obs'], ssk['tdur_obs']*60.)
    logxi_KS, logxi_KS_pos = KS_dist(np.log10(sss['xi_obs']), np.log10(ssk['xi_obs']))
    logxi_res_KS, logxi_res_KS_pos = KS_dist(np.log10(sss['xi_res_obs']), np.log10(ssk['xi_res_obs']))
    logxi_res32_KS, logxi_res32_KS_pos = KS_dist(np.log10(sss['xi_res32_obs']), np.log10(ssk['xi_res32_obs']))
    logxi_res21_KS, logxi_res21_KS_pos = KS_dist(np.log10(sss['xi_res21_obs']), np.log10(ssk['xi_res21_obs']))
    logxi_nonres_KS, logxi_nonres_KS_pos = KS_dist(np.log10(sss['xi_nonres_obs']), np.log10(ssk['xi_nonres_obs']))
    D_KS, D_KS_pos = KS_dist(sss['D_obs'], ssk['D_obs'])
    radii_KS, radii_KS_pos = KS_dist(sss['radii_obs'], ssk['radii_obs'])
    D_above_KS, D_above_KS_pos = KS_dist(sss['D_above_obs'], ssk['D_above_obs'])
    D_below_KS, D_below_KS_pos = KS_dist(sss['D_below_obs'], ssk['D_below_obs'])
    D_ratio_KS, D_ratio_KS_pos = KS_dist(sss['D_ratio_obs'], ssk['D_ratio_obs'])
    D_ratio_above_KS, D_ratio_above_KS_pos = KS_dist(sss['D_ratio_above_obs'], ssk['D_ratio_above_obs'])
    D_ratio_below_KS, D_ratio_below_KS_pos = KS_dist(sss['D_ratio_below_obs'], ssk['D_ratio_below_obs'])
    D_ratio_across_KS, D_ratio_across_KS_pos = KS_dist(sss['D_ratio_across_obs'], ssk['D_ratio_across_obs'])

    # AD distances:
    if AD_mod == 'true':
        AD_stat = AD_mod_dist
    elif AD_mod == 'false':
        AD_stat = AD_dist
    else:
        print 'Invalid input for AD_mod; must be a string true or false.'

    P_AD = AD_stat(sss['P_obs'], ssk['P_obs'])
    R_AD = AD_stat(sss['Rm_obs'], ssk['Rm_obs'])
    tdur_AD = AD_stat(sss['tdur_obs'], ssk['tdur_obs']*60.)
    logxi_AD = AD_stat(np.log10(sss['xi_obs']), np.log10(ssk['xi_obs']))
    logxi_res_AD = AD_stat(np.log10(sss['xi_res_obs']), np.log10(ssk['xi_res_obs']))
    logxi_nonres_AD = AD_stat(np.log10(sss['xi_nonres_obs']), np.log10(ssk['xi_nonres_obs']))
    D_AD = AD_stat(sss['D_obs'], ssk['D_obs'])
    radii_AD = AD_stat(sss['radii_obs'], ssk['radii_obs'])
    D_above_AD = AD_stat(sss['D_above_obs'], ssk['D_above_obs'])
    D_below_AD = AD_stat(sss['D_below_obs'], ssk['D_below_obs'])
    D_ratio_AD = AD_stat(sss['D_ratio_obs'], ssk['D_ratio_obs'])
    D_ratio_above_AD = AD_stat(sss['D_ratio_above_obs'], ssk['D_ratio_above_obs'])
    D_ratio_below_AD = AD_stat(sss['D_ratio_below_obs'], ssk['D_ratio_below_obs'])
    D_ratio_across_AD = AD_stat(sss['D_ratio_across_obs'], ssk['D_ratio_across_obs'])

    distances_KS = [delta_f, M_KS, d_mult_CRPD, d_mult_CRPD_switched, P_KS, R_KS, tdur_KS, logxi_KS, logxi_nonres_KS, logxi_res_KS, D_KS, D_above_KS, D_below_KS, D_ratio_KS, D_ratio_above_KS, D_ratio_below_KS, D_ratio_across_KS]
    distances_AD = [delta_f, M_KS, d_mult_CRPD, d_mult_CRPD_switched, P_AD, R_AD, tdur_AD, logxi_AD, logxi_nonres_AD, logxi_res_AD, D_AD, D_above_AD, D_below_AD, D_ratio_AD, D_ratio_above_AD, D_ratio_below_AD, D_ratio_across_AD]

    weights_file_name = 'Clustered_P_R_broken_R_weights_ADmod_%s_targs399675_evals1000.txt' % AD_mod
    weights_KS = load_weights('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name, 'KS')
    weights_AD = load_weights('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name, 'AD')

    distances_weighted_KS = distances_KS*weights_KS
    distances_weighted_AD = distances_AD*weights_AD

    delta_f_w, M_KS_w, d_mult_CRPD_w, d_mult_CRPD_switched_w, P_KS_w, R_KS_w, tdur_KS_w, logxi_KS_w, logxi_nonres_KS_w, logxi_res_KS_w, D_KS_w, D_above_KS_w, D_below_KS_w, D_ratio_KS_w, D_ratio_above_KS_w, D_ratio_below_KS_w, D_ratio_across_KS_w = distances_weighted_KS
    delta_f_w, M_KS_w, d_mult_CRPD_w, d_mult_CRPD_switched_w, P_AD_w, R_AD_w, tdur_AD_w, logxi_AD_w, logxi_nonres_AD_w, logxi_res_AD_w, D_AD_w, D_above_AD_w, D_below_AD_w, D_ratio_AD_w, D_ratio_above_AD_w, D_ratio_below_AD_w, D_ratio_across_AD_w = distances_weighted_AD

    weights_file_name_old = 'Clustered_P_R_broken_R_weights_ADmod_%s_targs696160_evals1000.txt' % AD_mod
    
    dist_weighted_KS = load_weights_and_compute_weighted_dist('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name, 'KS', distances_KS,  dists_exclude=dists_exclude, print_dists=print_dists)
    #dist_weighted_KS_old = load_weights_and_compute_weighted_dist('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name_old, 'KS', distances_KS,  dists_exclude=dists_exclude, print_dists=print_dists)

    dist_weighted_AD = load_weights_and_compute_weighted_dist('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name, 'AD', distances_AD,  dists_exclude=dists_exclude, print_dists=print_dists)
    #dist_weighted_AD_old = load_weights_and_compute_weighted_dist('ExoplanetsSysSim_Clusters/SysSimExClusters/src/' + weights_file_name_old, 'AD', distances_AD,  dists_exclude=dists_exclude, print_dists=print_dists)

    dists_misc = {'delta_f': delta_f, 'd_mult_CRPD': d_mult_CRPD, 'd_mult_CRPD_switched': d_mult_CRPD_switched, 'R_res32_diff': R_res32_diff, 'R_res21_diff': R_res21_diff}
    dists_misc_w = {'delta_f': delta_f_w, 'd_mult_CRPD': d_mult_CRPD_w, 'd_mult_CRPD_switched': d_mult_CRPD_switched_w}
    dists_KS = {'tot_weighted': dist_weighted_KS, 'M': M_KS, 'M_pos': M_KS_pos, 'P': P_KS, 'P_pos': P_KS_pos, 'Rm': R_KS, 'Rm_pos': R_KS_pos, 'tdur': tdur_KS, 'tdur_pos': tdur_KS_pos, 'logxi': logxi_KS, 'logxi_pos': logxi_KS_pos, 'logxi_res': logxi_res_KS, 'logxi_res_pos': logxi_res_KS_pos, 'logxi_res32': logxi_res32_KS, 'logxi_res32_pos': logxi_res32_KS_pos, 'logxi_res21': logxi_res21_KS, 'logxi_res21_pos': logxi_res21_KS_pos, 'logxi_nonres': logxi_nonres_KS, 'logxi_nonres_pos': logxi_nonres_KS_pos, 'D': D_KS, 'D_pos': D_KS_pos, 'radii': radii_KS, 'radii_pos': radii_KS_pos, 'D_above': D_above_KS, 'D_above_pos': D_above_KS_pos, 'D_below': D_below_KS, 'D_below_pos': D_below_KS_pos, 'D_ratio': D_ratio_KS, 'D_ratio_pos': D_ratio_KS_pos, 'D_ratio_above': D_ratio_above_KS, 'D_ratio_above_pos': D_ratio_above_KS_pos, 'D_ratio_below': D_ratio_below_KS, 'D_ratio_below_pos': D_ratio_below_KS_pos, 'D_ratio_across': D_ratio_across_KS, 'D_ratio_across_pos': D_ratio_across_KS_pos}
    dists_KS_w = {'M': M_KS_w, 'P': P_KS_w, 'Rm': R_KS_w, 'tdur': tdur_KS_w, 'logxi': logxi_KS_w, 'logxi_res': logxi_res_KS_w, 'logxi_nonres': logxi_nonres_KS_w, 'D': D_KS_w, 'D_above': D_above_KS_w, 'D_below': D_below_KS_w, 'D_ratio': D_ratio_KS_w, 'D_ratio_above': D_ratio_above_KS_w, 'D_ratio_below': D_ratio_below_KS_w, 'D_ratio_across': D_ratio_across_KS_w}
    dists_AD = {'tot_weighted': dist_weighted_AD, 'P': P_AD, 'Rm': R_AD, 'tdur': tdur_AD, 'logxi': logxi_AD, 'logxi_res': logxi_res_AD, 'logxi_nonres': logxi_nonres_AD, 'D': D_AD, 'radii': radii_AD, 'D_above': D_above_AD, 'D_below': D_below_AD, 'D_ratio': D_ratio_AD, 'D_ratio_above': D_ratio_above_AD, 'D_ratio_below': D_ratio_below_AD, 'D_ratio_across': D_ratio_across_AD}
    dists_AD_w = {'P': P_AD_w, 'Rm': R_AD_w, 'tdur': tdur_AD_w, 'logxi': logxi_AD_w, 'logxi_res': logxi_res_AD_w, 'logxi_nonres': logxi_nonres_AD_w, 'D': D_AD_w, 'D_above': D_above_AD_w, 'D_below': D_below_AD_w, 'D_ratio': D_ratio_AD_w, 'D_ratio_above': D_ratio_above_AD_w, 'D_ratio_below': D_ratio_below_AD_w, 'D_ratio_across': D_ratio_across_AD_w}
    return [dists_misc, dists_misc_w, dists_KS, dists_KS_w, dists_AD, dists_AD_w]





# Functions to make plots comparing the simulated and Kepler populations:

def setup_fig_single(fig_size, left, bottom, right, top, wspace=0, hspace=0):
    fig = plt.figure(figsize=fig_size)
    plot = GridSpec(1,1, left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    ax = plt.subplot(plot[0,0])
    return ax

def plot_panel_pdf_simple(ax, x_sim, x_Kep, x_min=None, x_max=None, y_min=None, y_max=None, n_bins=100, normalize=True, log_x=False, log_y=False, c_sim=['k'], c_Kep=['k'], ls_sim=['-'], ls_Kep=['-'], lw=1, alpha=0.2, labels_sim=['Simulated'], labels_Kep=['Kepler'], xticks_custom=None, xlabel_text='x', ylabel_text='Fraction', afs=20, tfs=20, lfs=16, legend=False):
    if x_min == None:
        x_min = np.nanmin([np.min(x) if len(x) > 0 else np.nan for x in x_sim+x_Kep])
    if x_max == None:
        x_max = np.nanmax([np.max(x) if len(x) > 0 else np.nan for x in x_sim+x_Kep])
    
    if log_x:
        bins = np.logspace(np.log10(x_min), np.log10(x_max), n_bins+1)
    else:
        bins = np.linspace(x_min, x_max, n_bins+1)

    for i,x in enumerate(x_sim):
        if normalize:
            weights = np.ones(len(x))/len(x)
        else:
            weights = np.ones(len(x))
        plt.hist(x, bins=bins, histtype='step', weights=weights, log=log_y, color=c_sim[i], ls=ls_sim[i], lw=lw, label=labels_sim[i])
    for i,x in enumerate(x_Kep):
        if normalize:
            weights = np.ones(len(x))/len(x)
        else:
            weights = np.ones(len(x))
        plt.hist(x, bins=bins, histtype='stepfilled', weights=weights, log=log_y, color=c_Kep[i], ls=ls_Kep[i], alpha=alpha, label=labels_Kep[i])

    if log_x:
        plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=afs)
    if xticks_custom != None:
        ax.set_xticks(xticks_custom)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])
    plt.xlabel(xlabel_text, fontsize=tfs)
    plt.ylabel(ylabel_text, fontsize=tfs)
    if legend:
        plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=lfs) #show the legend

def plot_fig_pdf_simple(fig_size, x_sim, x_Kep, x_min=None, x_max=None, y_min=None, y_max=None, n_bins=100, normalize=True, log_x=False, log_y=False, c_sim=['k'], c_Kep=['k'], ls_sim=['-'], ls_Kep=['-'], lw=1, alpha=0.2, labels_sim=['Simulated'], labels_Kep=['Kepler'], xticks_custom=None, xlabel_text='x', ylabel_text='Fraction', afs=20, tfs=20, lfs=16, legend=False, fig_lbrt=[0.15, 0.2, 0.95, 0.925], save_name='no_name_fig.pdf', save_fig=False):
    
    left, bottom, right, top = fig_lbrt
    ax = setup_fig_single(fig_size, left=left, bottom=bottom, right=right, top=top)
    
    plot_panel_pdf_simple(ax, x_sim, x_Kep, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, n_bins=n_bins, normalize=normalize, log_x=log_x, log_y=log_y, c_sim=c_sim, c_Kep=c_Kep, ls_sim=ls_sim, ls_Kep=ls_Kep, lw=lw, alpha=alpha, labels_sim=labels_sim, labels_Kep=labels_Kep, xticks_custom=xticks_custom, xlabel_text=xlabel_text, ylabel_text=ylabel_text, afs=afs, tfs=tfs, lfs=lfs, legend=legend)
    
    if save_fig:
        plt.savefig(save_name)
        plt.close()

def plot_fig_cdf_simple(fig_size, x_sim, x_Kep, x_min=None, x_max=None, log_x=False, c_sim=['k'], c_Kep=['k'], ls_sim=['-'], ls_Kep=['--'], lw=1, labels_sim=['Simulated'], labels_Kep=['Kepler'], xticks_custom=None, xlabel_text='x', ylabel_text='CDF', afs=20, tfs=20, lfs=16, legend=False, label_dist=False, fig_lbrt=[0.15, 0.2, 0.95, 0.925], save_name='no_name_fig.pdf', save_fig=False):
    
    left, bottom, right, top = fig_lbrt
    ax = setup_fig_single(fig_size, left=left, bottom=bottom, right=right, top=top)
    
    if x_min == None:
        x_min = np.nanmin([np.min(x) if len(x) > 0 else np.nan for x in x_sim+x_Kep])
    if x_max == None:
        x_max = np.nanmax([np.max(x) if len(x) > 0 else np.nan for x in x_sim+x_Kep])
    
    for i,x in enumerate(x_sim):
        plt.plot(np.sort(x), (np.arange(len(x))+1.)/np.float(len(x)), drawstyle='steps-post', color=c_sim[i], ls=ls_sim[i], lw=lw, label=labels_sim[i])
    for i,x in enumerate(x_Kep):
        plt.plot(np.sort(x), (np.arange(len(x))+1.)/np.float(len(x)), drawstyle='steps-post', color=c_Kep[i], ls=ls_Kep[i], lw=lw, label=labels_Kep[i])
    if label_dist:
        if len(x_sim) == len(x_Kep):
            for i in range(len(x_sim)):
                dist_KS, dist_KS_pos = KS_dist(x_sim[i], x_Kep[i])
                dist_AD = AD_mod_dist(x_sim[i], x_Kep[i])
                plt.annotate(s='', xy=(dist_KS_pos, sum(x_sim[i] <= dist_KS_pos)/np.float(len(x_sim[i]))), xytext=(dist_KS_pos, sum(x_Kep[i] <= dist_KS_pos)/np.float(len(x_Kep[i]))), arrowprops=dict(arrowstyle='<->'))
                plt.figtext(right-0.025, bottom+0.15+(len(x_sim)-(i+1.))*0.2, r'$\mathcal{D}_{\rm KS} = %s$' % np.round(dist_KS, 3), color=c_Kep[i], ha='right', fontsize=lfs)
                plt.figtext(right-0.025, bottom+0.05+(len(x_sim)-(i+1.))*0.2, r'$\mathcal{D}_{\rm AD^\prime} = %s$' % np.round(dist_AD, 3), color=c_Kep[i], ha='right', fontsize=lfs)
        else:
            print 'Error: x_sim != x_Kep'

    if log_x:
        plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=afs)
    if xticks_custom != None:
        ax.set_xticks(xticks_custom)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim([x_min, x_max])
    plt.xlabel(xlabel_text, fontsize=tfs)
    plt.ylabel(ylabel_text, fontsize=tfs)
    if legend:
        plt.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, frameon=False, fontsize=lfs) #show the legend

    if save_fig:
        plt.savefig(save_name)
        plt.close()

def load_cat_obs_and_plot_fig_pdf_composite(loadfiles_directory, run_number='', label_dist=True, AD_mod='true', dists_exclude=np.array([], dtype='i8'), n_bins=100, lw=1, alpha=0.2, afs=12, tfs=12, lfs=12, save_name='no_name_fig.pdf', save_fig=False):

    #To load and analyze the simulated and Kepler observed catalogs:
    
    N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)
    param_vals_all = read_sim_params(loadfiles_directory + 'periods%s.out' % run_number)
    
    sss_per_sys, sss = compute_summary_stats_from_cat_obs(file_name_path=loadfiles_directory, run_number=run_number)
    
    ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)
    
    dists_misc, dists_misc_w, dists_KS, dists_KS_w, dists_AD, dists_AD_w = compute_distances_sim_Kepler(loadfiles_directory, None, None, P_min, P_max, radii_min, radii_max, AD_mod=AD_mod, dists_exclude=dists_exclude, run_number=run_number)
    
    
    
    #To plot the 'observed' distributions with the actual observed Kepler distributions:
    
    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(4,3,left=0.075,bottom=0.075,right=0.975,top=0.975,wspace=0.15,hspace=0.5)
    
    #To print the parameter values:
    nrows = 7
    for i in range(len(param_keys_all)): #range(len(param_keys_all))
        plt.figtext(x=0.02+0.12*int(i/float(nrows)), y=0.95-0.025*(i%nrows), s=r'%s = %s' % (param_keys_all[i][1], param_vals_all[i]), fontsize=tfs)
    
    ax = plt.subplot(plot[1,0])
    plt.title(r'$\mathcal{D}_W({\rm KS}) = %1.2f$; $\mathcal{D}_W({\rm AD}) = %1.2f$' % (dists_KS['tot_weighted'], dists_AD['tot_weighted']), fontsize=lfs)
    x = sss['Mtot_obs'][sss['Mtot_obs'] > 0]
    max_M = np.max((np.max(sss['Mtot_obs']), np.max(ssk['Mtot_obs'])))
    counts, bins = np.histogram(x, bins=max_M+1, range=(-0.5, max_M+0.5))
    bins_mid = (bins[:-1] + bins[1:])/2.
    plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color='k', lw=lw, label='%s x5 pl (Sim)' % (sum(x)/5.))
    counts, bins = np.histogram(ssk['Mtot_obs'], bins=bins)
    plt.plot(bins_mid, counts/float(np.sum(counts)), 'o--', color='k', alpha=alpha, label='%s pl (Kep)' % sum(ssk['Mtot_obs']))
    plt.gca().set_yscale("log")
    ax.tick_params(axis='both', labelsize=afs)
    plt.xlim([1., max_M])
    #plt.xlim([1., 8.])
    plt.xlabel('Observed planets per system', fontsize=tfs)
    plt.ylabel('Fraction', fontsize=tfs)
    plt.legend(loc='lower left', bbox_to_anchor=(0.01,0.01), frameon=False, ncol=1, fontsize=lfs) #show the legend
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$|f_{\rm sim} - f_{\rm Kep}| = %1.4f$ ($%1.2f$)' % (dists_misc['delta_f'], dists_misc_w['delta_f']), ha='right', fontsize=lfs, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['M'], dists_KS_w['M']), ha='right', fontsize=lfs, transform = ax.transAxes)
        plt.text(x=0.98, y=0.4, s=r'${\rm CRPD} = %1.4f$ ($%1.2f$)' % (dists_misc['d_mult_CRPD'], dists_misc_w['d_mult_CRPD']), ha='right', fontsize=lfs, transform = ax.transAxes)
        plt.text(x=0.98, y=0.2, s=r'${\rm CRPD_r} = %1.4f$ ($%1.2f$)' % (dists_misc['d_mult_CRPD_switched'], dists_misc_w['d_mult_CRPD_switched']), ha='right', fontsize=lfs, transform = ax.transAxes)
    
    ax = plt.subplot(plot[2,0])
    plot_panel_pdf_simple(ax, [sss['P_obs']], [ssk['P_obs']], x_min=P_min, x_max=P_max, y_min=1e-3, y_max=0.1, n_bins=n_bins, log_x=True, log_y=True, lw=lw, alpha=alpha, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['P'], dists_KS_w['P']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['P'], dists_AD_w['P']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[3,0])
    R_max_cut = 30. #upper cut-off for plotting period ratios; np.max(sss['Rm_obs'])
    plot_panel_pdf_simple(ax, [sss['Rm_obs'][sss['Rm_obs'] < R_max_cut]], [ssk['Rm_obs'][ssk['Rm_obs'] < R_max_cut]], x_min=1., x_max=R_max_cut, n_bins=n_bins, log_x=True, lw=lw, alpha=alpha, xticks_custom=[1,2,3,4,5,10,20], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['Rm'], dists_KS_w['Rm']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['Rm'], dists_AD_w['Rm']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[0,1])
    plot_panel_pdf_simple(ax, [sss['tdur_obs']], [ssk['tdur_obs']*60.], x_max=1000., n_bins=n_bins, lw=lw, alpha=alpha, xlabel_text=r'$t_{\rm dur}$ (mins)', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['tdur'], dists_KS_w['tdur']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['tdur'], dists_AD_w['tdur']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[1,1])
    plot_panel_pdf_simple(ax, [np.log10(sss['xi_obs'])], [np.log10(ssk['xi_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, lw=lw, alpha=alpha, xlabel_text=r'$\log{\xi}$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['logxi'], dists_KS_w['logxi']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['logxi'], dists_AD_w['logxi']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[2,1])
    plot_panel_pdf_simple(ax, [sss['D_obs']], [ssk['D_obs']], x_min=1e-5, x_max=1e-2, n_bins=n_bins, log_x=True, lw=lw, alpha=alpha, xlabel_text=r'$\delta$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D'], dists_KS_w['D']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D'], dists_AD_w['D']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[3,1])
    plot_panel_pdf_simple(ax, [sss['D_above_obs'], sss['D_below_obs']], [ssk['D_above_obs'], ssk['D_below_obs']], x_min=1e-5, x_max=1e-2, n_bins=n_bins, log_x=True, c_sim=['b','r'], c_Kep=['b','r'], ls_sim=['-','-'], ls_Kep=['-','-'], lw=lw, alpha=alpha, labels_sim=['Above', 'Below'], labels_Kep=[None, None], xlabel_text=r'$\delta$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    plt.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, frameon=False, fontsize=12) #show the legend
    if label_dist:
        plt.text(x=0.98, y=0.85, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_above'], dists_KS_w['D_above']), ha='right', color='b', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.7, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_above'], dists_AD_w['D_above']), ha='right', color='b', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.55, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_below'], dists_KS_w['D_below']), ha='right', color='r', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.4, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_below'], dists_AD_w['D_below']), ha='right', color='r', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[0,2])
    plot_panel_pdf_simple(ax, [sss['D_ratio_obs']], [ssk['D_ratio_obs']], x_min=0.1, x_max=10., n_bins=n_bins, log_x=True, lw=lw, alpha=alpha, xlabel_text=r'$\delta_{i+1}/\delta_i$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    if label_dist:
        plt.text(x=0.98, y=0.8, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_ratio'], dists_KS_w['D_ratio']), ha='right', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.6, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_ratio'], dists_AD_w['D_ratio']), ha='right', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[1,2])
    plot_panel_pdf_simple(ax, [sss['D_ratio_above_obs'], sss['D_ratio_below_obs'], sss['D_ratio_across_obs']], [ssk['D_ratio_above_obs'], ssk['D_ratio_below_obs'], ssk['D_ratio_across_obs']], x_min=0.1, x_max=10., n_bins=n_bins, log_x=True, c_sim=['b','r','k'], c_Kep=['b','r','k'], ls_sim=['-','-','-'], ls_Kep=['-','-','-'], lw=lw, alpha=alpha, labels_sim=['Above', 'Below', 'Across'], labels_Kep=[None, None, None], xlabel_text=r'$\delta_{i+1}/\delta_i$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    plt.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, frameon=False, fontsize=12) #show the legend
    if label_dist:
        plt.text(x=0.98, y=0.85, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_ratio_above'], dists_KS_w['D_ratio_above']), ha='right', color='b', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.7, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_ratio_above'], dists_AD_w['D_ratio_above']), ha='right', color='b', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.55, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_ratio_below'], dists_KS_w['D_ratio_below']), ha='right', color='r', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.4, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_ratio_below'], dists_AD_w['D_ratio_below']), ha='right', color='r', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.25, s=r'$\mathcal{D}_{\rm KS} = %1.4f$ ($%1.2f$)' % (dists_KS['D_ratio_across'], dists_KS_w['D_ratio_across']), ha='right', color='k', fontsize=12, transform = ax.transAxes)
        plt.text(x=0.98, y=0.1, s=r'$\mathcal{D}_{\rm AD} = %1.4f$ ($%1.2f$)' % (dists_AD['D_ratio_across'], dists_AD_w['D_ratio_across']), ha='right', color='k', fontsize=12, transform = ax.transAxes)
    
    ax = plt.subplot(plot[2,2])
    plot_panel_pdf_simple(ax, [sss['radii_obs']], [ssk['radii_obs']], x_min=radii_min, x_max=radii_max, n_bins=n_bins, lw=lw, alpha=alpha, xlabel_text=r'$R_p (R_\oplus)$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    
    ax = plt.subplot(plot[3,2])
    plot_panel_pdf_simple(ax, [sss['Rstar_obs']], [ssk['Rstar_obs']], x_max=3., n_bins=n_bins, lw=lw, alpha=alpha, xlabel_text=r'$R_\star (R_\odot)$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)
    
    if save_fig:
        plt.savefig(save_name)
        plt.close()

def load_cat_obs_and_plot_figs_multis_gallery(loadfiles_directory, run_number='', save_name_base='no_name_fig', save_fig=False):

    #To load and analyze the simulated and Kepler observed catalogs:
    
    N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)
    param_vals_all = read_sim_params(loadfiles_directory + 'periods%s.out' % run_number)
    
    sss_per_sys, sss = compute_summary_stats_from_cat_obs(file_name_path=loadfiles_directory, run_number=run_number)
    
    ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)



    # To plot the observed multi-systems by period to visualize the systems (similar to Fig 1 in Fabrycky et al. 2014):
    N_multi = sum(sss['Mtot_obs'] >= 3) #number of simulated multi-systems with 3 or more planets
    N_multi_confirmed = sum(ssk['Mtot_obs'] >= 3)

    i_sorted_P0 = np.argsort(sss_per_sys['P_obs'][sss['Mtot_obs'] >= 3,0]) #array of indices that would sort the arrays of multi-systems by the innermost period of each system
    i_sorted_P0 = i_sorted_P0[np.sort(np.random.choice(np.arange(len(i_sorted_P0)), int(round(N_multi/(N_sim/N_Kep))), replace=False))]
    P_obs_multi = sss_per_sys['P_obs'][sss['Mtot_obs'] >= 3][i_sorted_P0]
    radii_obs_multi = sss_per_sys['radii_obs'][sss['Mtot_obs'] >= 3][i_sorted_P0]

    i_sorted_P0_confirmed = np.argsort(ssk_per_sys['P_obs'][ssk['Mtot_obs'] >= 3,0]) #array of indices that would sort the arrays of multi-systems by the innermost period of each system
    P_obs_multi_confirmed = ssk_per_sys['P_obs'][ssk['Mtot_obs'] >= 3][i_sorted_P0_confirmed]
    radii_obs_multi_confirmed = ssk_per_sys['radii_obs'][ssk['Mtot_obs'] >= 3][i_sorted_P0_confirmed]

    N_sys_per_plot = 150 #number of systems to plot per figure
    for i in range(int(np.ceil(float(len(i_sorted_P0))/N_sys_per_plot))):
        fig = plt.figure(figsize=(10,10))
        plot = GridSpec(1,2,left=0.025,bottom=0.1,right=0.975,top=0.95,wspace=0,hspace=0.1)
        
        ax = plt.subplot(plot[0,0])
        plt.title('Kepler multi-planet systems', fontsize=12)
        for j in range(len(P_obs_multi_confirmed[i*N_sys_per_plot:(i+1)*N_sys_per_plot])):
            P_sys = P_obs_multi_confirmed[i*N_sys_per_plot + j]
            radii_sys = radii_obs_multi_confirmed[i*N_sys_per_plot + j]
            P_sys = P_sys[P_sys > 0]
            radii_sys = radii_sys[radii_sys > 0]
            plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=2.*radii_sys**2.)
            if (j+1)%10 == 0:
                plt.axhline(y=j+1, lw=0.05, color='k')
        plt.gca().set_xscale("log")
        ax.set_xticks([3,10,30,100,300])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_yticks([])
        plt.xlim([2., 500.])
        plt.ylim([0., N_sys_per_plot])
        plt.xlabel(r'$P$ (days)', fontsize=12)

        ax = plt.subplot(plot[0,1])
        plt.title('Simulated multi-planet systems', fontsize=12)
        for j in range(len(P_obs_multi[i*N_sys_per_plot:(i+1)*N_sys_per_plot])):
            P_sys = P_obs_multi[i*N_sys_per_plot + j]
            radii_sys = radii_obs_multi[i*N_sys_per_plot + j]
            P_sys = P_sys[P_sys > 0]
            radii_sys = radii_sys[radii_sys > 0]
            plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=2.*radii_sys**2.)
            if (j+1)%10 == 0:
                plt.axhline(y=j+1, lw=0.05, color='k')
        plt.gca().set_xscale("log")
        ax.set_xticks([3,10,30,100,300])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_yticks([])
        plt.xlim([2., 500.])
        plt.ylim([0., N_sys_per_plot])
        plt.xlabel(r'$P$ (days)', fontsize=12)

        save_name = save_name_base + '_%s.pdf' % i
        if save_fig:
            plt.savefig(save_name)
            plt.close()

def load_cat_obs_and_plot_fig_period_radius(loadfiles_directory, run_number='', lw=1, save_name='no_name_fig.pdf', save_fig=False):

    #To load and analyze the simulated and Kepler observed catalogs:
    
    N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)
    param_vals_all = read_sim_params(loadfiles_directory + 'periods%s.out' % run_number)
    
    sss_per_sys, sss = compute_summary_stats_from_cat_obs(file_name_path=loadfiles_directory, run_number=run_number)
    
    ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)
    
    
    
    #To plot a period vs radius scatter plot with binned statistics to compare the simulated and Kepler catalogs:
    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(1,2,left=0.075,bottom=0.1,right=0.975,top=0.75,wspace=0.2,hspace=0)
    
    #To print the parameter values:
    nrows = 7
    for i in range(len(param_keys_all)): #range(len(param_keys_all))
        plt.figtext(x=0.02+0.12*int(i/float(nrows)), y=0.95-0.025*(i%nrows), s=r'%s = %s' % (param_keys_all[i][1], param_vals_all[i]), fontsize=12)
    
    P_bins = 5
    P_lines, radii_lines = np.logspace(np.log10(P_min), np.log10(P_max), P_bins+1), np.array([0.5, 1., 2., 4., 6., 8., 10.])
    radii_bins = len(radii_lines)-1
    
    ax = plt.subplot(plot[0,0])
    N_sample = int(np.round(len(sss_per_sys['P_obs'])*cos_factor)) #number of simulated planets we would expect if we assigned orbits isotropically
    i_sample = np.random.choice(np.arange(len(sss_per_sys['P_obs'])), N_sample, replace=False)
    plt.scatter(sss['P_obs'][i_sample], sss['radii_obs'][i_sample], c='k', marker='o')
    plt.scatter(ssk['P_obs'], ssk['radii_obs'], c='r', marker='o')
    for x in P_lines:
        plt.axvline(x, lw=lw, color='g')
    for y in radii_lines:
        plt.axhline(y, lw=lw, color='g')
    plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=20)
    ax.set_xticks([3,10,30,100,300])
    ax.set_yticks([0.5,2,4,6,8,10])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim([P_min, P_max])
    plt.ylim([radii_min, radii_max])
    plt.xlabel(r'$P$ (days)', fontsize=20)
    plt.ylabel(r'$R_p (R_\oplus)$', fontsize=20)

    ax = plt.subplot(plot[0,1])
    N_obs_grid = np.zeros((radii_bins, P_bins))
    N_confirmed_grid = np.zeros((radii_bins, P_bins))
    for j in range(radii_bins):
        for i in range(P_bins):
            N_obs_cell = np.sum((sss['P_obs'] > P_lines[i]) & (sss['P_obs'] < P_lines[i+1]) & (sss['radii_obs'] > radii_lines[j]) & (sss['radii_obs'] < radii_lines[j+1]))
            N_confirmed_cell = np.sum((ssk['P_obs'] > P_lines[i]) & (ssk['P_obs'] < P_lines[i+1]) & (ssk['radii_obs'] > radii_lines[j]) & (ssk['radii_obs'] < radii_lines[j+1]))
            N_obs_grid[j,i] = N_obs_cell
            N_confirmed_grid[j,i] = N_confirmed_cell
            plt.text(x=0.02+i*(1./P_bins), y=(j+1)*(1./radii_bins)-0.025, s='%s' % np.round(N_obs_cell*cos_factor, 1), ha='left', va='top', color='k', fontsize=16, transform = ax.transAxes)
            plt.text(x=0.02+i*(1./P_bins), y=(j+1)*(1./radii_bins)-0.075, s='%s' % N_confirmed_cell, ha='left', va='top', color='r', fontsize=16, transform = ax.transAxes)
            plt.text(x=0.02+i*(1./P_bins), y=(j+1)*(1./radii_bins)-0.125, s='%s' % np.round((N_obs_cell*cos_factor)/float(N_confirmed_cell), 2), ha='left', va='top', color='b', fontsize=16, fontweight='bold', transform = ax.transAxes)
    N_obs_normed_grid = N_obs_grid*cos_factor
    plt.imshow(N_obs_normed_grid/N_confirmed_grid, cmap='coolwarm', aspect='auto', interpolation="nearest", origin='lower') #extent=(3, 300, 0.5, 10)
    cbar = plt.colorbar()
    cbar.set_label(r'$N_{\rm Sim}/N_{\rm Kep}$', rotation=270, va='bottom', fontsize=20)
    #plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=20)
    plt.xticks(np.linspace(-0.5, P_bins-0.5, P_bins), [3,10,30,100,300])
    plt.yticks(np.linspace(-0.5, radii_bins-0.5, radii_bins+1), radii_lines)
    plt.xlabel(r'$P$ (days)', fontsize=20)

    plot = GridSpec(1,1,left=0.83,bottom=0.8,right=0.895,top=0.93,wspace=0,hspace=0) #just for the 'legend'
    ax = plt.subplot(plot[0,0])
    plt.text(x=0.05, y=0.9, s=r'$N_{\rm Sim}$', ha='left', va='top', color='k', fontsize=14, transform = ax.transAxes)
    plt.text(x=0.05, y=0.7, s=r'$N_{\rm Kep}$', ha='left', va='top', color='r', fontsize=14, transform = ax.transAxes)
    plt.text(x=0.05, y=0.5, s=r'$N_{\rm Sim}/N_{\rm Kep}$', ha='left', va='top', color='b', fontsize=14, transform = ax.transAxes)
    plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')

    if save_fig:
        plt.savefig(save_name)
        plt.close()





# Functions to analyze GP models/outputs and to make corner plots for visualizing the parameter space:

def load_training_points(dims, file_name_path='', best=100000, every=10):
    # 'dims' is the number of dimensions (model parameters); 'best' is the number of total best-ranked points from which this file of points was generated; 'every' is how many per of the 'best' ranked points are actually saved in the file

    data_table_recomputed = np.genfromtxt(file_name_path + 'Active_params_distances_recomputed_table_best%s_every%s.txt' % (best, every), delimiter=' ', names=True, dtype='f8')
    cols = len(data_table_recomputed.dtype.names)
    
    active_params_names = data_table_recomputed.dtype.names[:dims]
    print 'Active param names: ', active_params_names
    
    #xtrain = data_table_recomputed[np.array(active_params_names)]
    #xtrain = xtrain.view((float, dims))
    xtrain = data_table_recomputed.view(np.float64).reshape(data_table_recomputed.shape + (cols,))
    xtrain = xtrain[:,:dims]

    ytrain = data_table_recomputed['dist_tot_weighted']

    data_train = {'active_params_names': active_params_names, 'xtrain': xtrain, 'ytrain': ytrain}
    return data_train

def load_GP_table_prior_draws(file_name, file_name_path=''):
    #Example file_name: 'GP_emulator_points%s_meanf%s_small_hparams_prior_draws%s.csv' % (n_train, mean_f, n_draws)
    xprior_accepted_table = np.genfromtxt(file_name_path + file_name, skip_header=4, delimiter=',', names=True, dtype='f8')
    return xprior_accepted_table

def load_table_points_min_GP(file_name, file_name_path=''):
    #Example file_name: 'GP_train2000_meanf50.0_minimize_mean_iterations10.csv'
    xmin_table = np.genfromtxt(file_name_path + file_name, skip_header=2, delimiter=',', names=True, dtype='f8')
    return xmin_table

def load_GP_2d_grids(dims, n_train, mean_f, sigma_f, lscales, file_name_path='', grid_dims=50):
    file_name_mean = 'GP_train%s_meanf%s_sigmaf%s_lscales%s_grids2d_%sx%s_mean.csv' % (n_train, mean_f, sigma_f, lscales, grid_dims, grid_dims)
    file_name_std = 'GP_train%s_meanf%s_sigmaf%s_lscales%s_grids2d_%sx%s_std.csv' % (n_train, mean_f, sigma_f, lscales, grid_dims, grid_dims)
    GP_mean_2d_grids = np.genfromtxt(file_name_path + file_name_mean, delimiter=',', skip_header=5, dtype='f8')
    GP_std_2d_grids = np.genfromtxt(file_name_path + file_name_std, delimiter=',', skip_header=5, dtype='f8')
    with open(file_name_path + file_name_mean, 'r') as f:
        for line in f:
            if line[0:8] == '# xlower':
                xlower = [float(x) for x in line[12:-2].split(' ')]
            if line[0:8] == '# xupper':
                xupper = [float(x) for x in line[12:-2].split(' ')]
            if line[0:12] == '# xmin_guess':
                xmin_guess = [float(x) for x in line[16:-2].split(' ')]

    if grid_dims != np.shape(GP_mean_2d_grids)[1]:
        print 'PROBLEM: mismatch with grid_dims!'
    GP_mean_2d_grids = np.reshape(GP_mean_2d_grids, (sum(range(dims)), grid_dims, grid_dims))
    GP_std_2d_grids = np.reshape(GP_std_2d_grids, (sum(range(dims)), grid_dims, grid_dims))

    GP_grids = {'xlower': xlower, 'xupper': xupper, 'xmin_guess': xmin_guess, 'mean_grids': GP_mean_2d_grids, 'std_grids': GP_std_2d_grids}
    return GP_grids

def transform_sum_diff_params(xpoints, i, j):
    print 'Transforming columns (i,j) to (i+j,j-i).'
    xpoints_transformed = np.copy(xpoints)
    xpoints_transformed[:,i], xpoints_transformed[:,j] = xpoints_transformed[:,i] + xpoints_transformed[:,j], xpoints_transformed[:,j] - xpoints_transformed[:,i]
    return xpoints_transformed

def transform_sum_diff_params_inverse(xpoints, i, j):
    print 'Transforming columns (i,j) to ((i-j)/2,(i+j)/2).'
    xpoints_transformed = np.copy(xpoints)
    xpoints_transformed[:,i], xpoints_transformed[:,j] = (xpoints_transformed[:,i] - xpoints_transformed[:,j])/2., (xpoints_transformed[:,i] + xpoints_transformed[:,j])/2.
    return xpoints_transformed

def make_cuts_GP_mean_std_post(x_names, xprior_table, max_mean=np.inf, max_std=np.inf, max_post=np.inf):
    dims = len(x_names)
    
    xpoints_all = xprior_table[np.array(x_names)]
    xpoints_cut = xpoints_all[(xprior_table['GP_mean'] < max_mean) & (xprior_table['GP_std'] < max_std) & (xprior_table['GP_posterior_draw'] < max_post)]
    #xpoints_cut = xpoints_cut.view((float, dims))
    xpoints_cut = xpoints_cut.view(np.float64).reshape(xpoints_cut.shape + (dims+3,))
    xpoints_cut = xpoints_cut[:,:dims]
    print 'Total points: %s; points left: %s' % (len(xpoints_all), len(xpoints_cut))

    return xpoints_cut

def plot_fig_hists_GP_draws(fig_size, xprior_table, bins=100, save_name='no_name_fig.pdf', save_fig=False):
    fig = plt.figure(figsize=fig_size)
    plot = GridSpec(2,2,left=0.1,bottom=0.1,right=0.95,top=0.95,wspace=0.2,hspace=0.2)
    
    ax = plt.subplot(plot[0,0])
    plt.hist(xprior_table['GP_mean'], bins=bins)
    plt.xlabel('GP mean prediction')
    plt.ylabel('Points')
    
    ax = plt.subplot(plot[0,1])
    plt.hist(xprior_table['GP_posterior_draw'], bins=bins)
    plt.xlabel('GP prediction')
    plt.ylabel('Points')
    
    ax = plt.subplot(plot[1,0])
    plt.hist(xprior_table['GP_std'], bins=bins)
    plt.xlabel('GP std of prediction')
    plt.ylabel('Points')
    
    ax = plt.subplot(plot[1,1])
    plt.scatter(xprior_table['GP_mean'], xprior_table['GP_std'], marker='.')
    plt.xlabel('GP mean prediction')
    plt.ylabel('GP std of prediction')
    
    if save_fig:
        plt.savefig(save_name)
        plt.close()

def plot_cornerpy_wrapper(x_symbols, xpoints, xpoints_extra=None, c_extra='r', s_extra=1, quantiles=[0.16, 0.5, 0.84], show_titles=True, label_kwargs={'fontsize': 20}, title_kwargs={'fontsize':20}, save_name='no_name_fig.pdf', save_fig=False):

    dims = len(x_symbols)
    
    fig = corner.corner(xpoints, labels=x_symbols, quantiles=quantiles, show_titles=show_titles, label_kwargs=label_kwargs, title_kwargs=title_kwargs)

    # If want to plot an additional set of points:
    if xpoints_extra != None:
        axes = np.array(fig.axes).reshape((dims, dims))
        for i in range(dims):
            for j in range(i):
                ax = axes[i,j]
                ax.scatter(xpoints_extra[:,j], xpoints_extra[:,i], color=c_extra, s=s_extra)

    if save_fig:
        plt.savefig(save_name)
        plt.close()

def plot_contours_and_points_corner(x_symbols, xlower, xupper, contour_2d_grids, xpoints=None, points_size=1., points_alpha=1., fig_size=(16,16), fig_lbrtwh=[0.05, 0.05, 0.98, 0.98, 0.05, 0.05], afs=10, tfs=12, lfs=10, save_name='no_name_fig.pdf', save_fig=False):
    
    dims = len(x_symbols)
    grid_dims = np.shape(contour_2d_grids)[2]
    
    fig = plt.figure(figsize=fig_size)
    left, bottom, right, top, wspace, hspace = fig_lbrtwh
    plot = GridSpec(dims, dims, left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    grid_index = 0
    for i in range(dims): # indexes the row or y-axis variable
        for j in range(i): # indexes the column or x-axis variable
            xaxis = np.linspace(xlower[j], xupper[j], grid_dims)
            yaxis = np.linspace(xlower[i], xupper[i], grid_dims)
            xgrid, ygrid = np.meshgrid(xaxis, yaxis)
            
            ax = plt.subplot(plot[i,j])
            cplot = plt.contour(xgrid, ygrid, contour_2d_grids[grid_index])
            plt.clabel(cplot, inline=1, fontsize=lfs)
            if xpoints != None:
                plt.scatter(xpoints[:,j], xpoints[:,i], s=points_size, alpha=points_alpha, c='k')
            ax.tick_params(labelleft=False, labelbottom=False)
            if i == dims-1:
                ax.tick_params(labelbottom=True)
                plt.xlabel(x_symbols[j], fontsize=tfs)
            if j == 0:
                ax.tick_params(labelleft=True)
                plt.ylabel(x_symbols[i], fontsize=tfs)
            plt.xticks(rotation=45, fontsize=afs)
            plt.yticks(rotation=45, fontsize=afs)

            grid_index += 1

    for i in range(dims):
        ax = plt.subplot(plot[i,i])
        ax.tick_params(labelleft=False, labelbottom=False)
        if i == dims-1:
            #ax.tick_params(labelbottom=True)
            plt.xlabel(x_symbols[i], fontsize=tfs)
        if i == 0:
            #ax.tick_params(labelleft=True)
            plt.ylabel(x_symbols[i], fontsize=tfs)
        plt.xticks(rotation=45, fontsize=afs)
        plt.yticks(rotation=45, fontsize=afs)

    if save_fig:
        plt.savefig(save_name)
        plt.close()

def plot_function_heatmap_contours_given_irregular_points_corner(x_symbols, xpoints, fpoints, xlower=None, xupper=None, show_points=True, points_size=1., points_alpha=1., fig_size=(16,16), fig_lbrtwh=[0.05, 0.05, 0.98, 0.98, 0.05, 0.05], afs=10, tfs=12, lfs=10, save_name='no_name_fig.pdf', save_fig=False):

    dims = len(x_symbols)
    if xlower == None:
        xlower = np.min(xpoints, axis=0)
    if xupper == None:
        xupper = np.max(xpoints, axis=0)
    
    fig = plt.figure(figsize=fig_size)
    left, bottom, right, top, wspace, hspace = fig_lbrtwh
    plot = GridSpec(dims, dims, left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    grid_index = 0
    for i in range(dims): # indexes the row or y-axis variable
        for j in range(i): # indexes the column or x-axis variable
            ax = plt.subplot(plot[i,j])
            #plt.tricontour(xpoints[:,j], xpoints[:,i], fpoints, linewidths=0.5, colors='k')
            contours = plt.tricontourf(xpoints[:,j], xpoints[:,i], fpoints, cmap="RdBu_r")
            #plt.colorbar(contours)
            if show_points:
                plt.scatter(xpoints[:,j], xpoints[:,i], s=points_size, alpha=points_alpha, c='k')
            ax.axis((xlower[j], xupper[j], xlower[i], xupper[i]))
            ax.tick_params(labelleft=False, labelbottom=False)
            if i == dims-1:
                ax.tick_params(labelbottom=True)
                plt.xlabel(x_symbols[j], fontsize=tfs)
            if j == 0:
                ax.tick_params(labelleft=True)
                plt.ylabel(x_symbols[i], fontsize=tfs)
            plt.xticks(rotation=45, fontsize=afs)
            plt.yticks(rotation=45, fontsize=afs)
    
        grid_index += 1

    for i in range(dims):
        ax = plt.subplot(plot[i,i])
        plt.colorbar(contours)
        ax.tick_params(labelleft=False, labelbottom=False)
        if i == dims-1:
            #ax.tick_params(labelbottom=True)
            plt.xlabel(x_symbols[i], fontsize=tfs)
        if i == 0:
            #ax.tick_params(labelleft=True)
            plt.ylabel(x_symbols[i], fontsize=tfs)
            plt.xticks(rotation=45, fontsize=afs)
            plt.yticks(rotation=45, fontsize=afs)

    if save_fig:
        plt.savefig(save_name)
        plt.close()

