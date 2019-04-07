#using ExoplanetsSysSim
#using StatsFuns
#if !@isdefined JLD using JLD end
if !@isdefined DataFrames
    using DataFrames
end
if !@isdefined CSV
    using CSV
end
#import DataFrames.DataFrame, DataFrames.isna
import DataFrames.DataFrame, DataFrames.ismissing
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable

## Old code for generating stellar properties # TODO: WARNING: Should eventually use version in main branch to make sure have recent improvements

## stellar_table
function setup_star_table_christiansen(sim_param::SimParam; force_reread::Bool = false)
  #global df
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)), ".." , "data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv"))))
  df = setup_star_table_christiansen(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  ExoplanetsSysSim.StellarTable.set_star_table(df)
  return df  
end

function setup_star_table_christiansen(filename::String; force_reread::Bool = false)
  #global df, usable
  #println("Testing")
  wf = WindowFunction.setup_window_function(sim_param)
  #println("Testing 2")
  df = ExoplanetsSysSim.StellarTable.df
  #usable = ExoplanetsSysSim.StellarTable.usable
#if ismatch(r".jld$",filename)
  if occursin(r".jld$", filename)
  try
    data = JLD.load(filename)
    df::DataFrame = data["stellar_catalog"]
    usable::Array{Int64,1} = data["stellar_catalog_usable"]
    ExoplanetsSysSim.StellarTable.set_star_table(df, usable)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end

  else
  try
    #println(filename)
    df = CSV.read(filename,allowmissing=:all)
    #println("Testing 4")
    #flush(STDOUT)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end
  end # if ismatch

  has_mass = .! (ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
  has_radius = .! (ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))
  has_dens = .! (ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))
  has_rest = .! (ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))

#=
  in_Q1Q12 = []
  for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
=#

#Comment out following block for 'is_FGK' if using stellar catalogs with cuts already made:
#=
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
#println("logg?: ", df[x,:logg])
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
=#

  is_usable = has_radius .& has_mass .& has_rest #.& in_Q1Q12 #.& has_dens #.& is_FGK
  println("Total number of stars in our stellar catalog: ", sum(is_usable))
  #if contains(filename,"q1_q12_christiansen.jld")
#if contains(filename,"q1_q12_christiansen")   # TODO: Ask Danely what he's trying to do here.
#is_usable = is_usable #& in_Q1Q12
#end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :dataspan, :dutycycle ]
  deletecols!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = findall(is_usable)
  df = df[usable, symbols_to_keep]
  #ExoplanetsSysSim.StellarTable.set_star_table(df, usable)
  ExoplanetsSysSim.StellarTable.set_star_table(df)
#end
  return df
end

function test_stellar_table() # TODO: Write test for stellar catalog functions
end

