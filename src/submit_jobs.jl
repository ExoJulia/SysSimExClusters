"""
Given a file "f", writes our default PBS script settings to the file.
"""
function write_pbs_settings(f)
    println(f, "#!/bin/tcsh")
    println(f, "#PBS -A ebf11_a_g_sc_default")
    println(f, "#PBS -l nodes=1:ppn=1")
    println(f, "#PBS -l walltime=48:00:00")
    println(f, "#PBS -l pmem=8gb")
    println(f, "#PBS -j oe")
    #println(f, "#PBS -m abe")
    #println(f, "#PBS -M myh7@psu.edu")
    println(f, "")
    println(f, "cd \$PBS_O_WORKDIR")
    println(f, "")
    println(f, "setenv JULIA_DEPOT_PATH /gpfs/group/ebf11/default/myh7/julia_pkgdir/") # Point to where we installed things for Julia and are developing ExoplanetsSysSim
end

"""
Generates a PBS script for running "optimize.jl".
"""
function generate_pbs_optimize(run_number)
    f_name = "optimize_job_"*string(run_number)*".pbs"
    f = open(f_name, "w")
    write_pbs_settings(f)
    println(f, "/gpfs/group/ebf11/default/sw/julia-0.7.0/bin/julia optimize.jl "*string(run_number))
    close(f)

    return f_name
end

"""
Generates and submits "n_jobs" of PBS scripts to run "optimize.jl".
"""
function submit_jobs_optimize(n_jobs::Int64)
    for i in 1:n_jobs
        f_name = generate_pbs_optimize(i)
        run(`qsub $f_name`) #this line submits the job by running 'qsub' in the command line!
        println("Job ", f_name, " submitted.")
    end
end

"""
Generates a PBS script for running "generate_pbs_compute_distances_given_params_random.jl".
"""
function generate_pbs_compute_distances_given_params_random(run_number)
    f_name = "compute_distances_random_job_"*string(run_number)*".pbs"
    f = open(f_name, "w")
    write_pbs_settings(f)
    println(f, "/gpfs/group/ebf11/default/sw/julia-0.7.0/bin/julia compute_distances_given_params_random.jl "*string(run_number))
    close(f)

    return f_name
end

"""
Generates and submits "n_jobs" of PBS scripts to run "generate_pbs_compute_distances_given_params_random.jl".
"""
function submit_jobs_compute_distances_given_params_random(n_jobs::Int64)
    for i in 1:n_jobs
        f_name = generate_pbs_compute_distances_given_params_random(i)
        run(`qsub $f_name`) #this line submits the job by running 'qsub' in the command line!
        println("Job ", f_name, " submitted.")
    end
end

"""
Generates a PBS script for running "generate_pbs_compute_distances_given_params_parallel.jl".
"""
function generate_pbs_compute_distances_given_params_parallel(run_number)
    f_name = "compute_distances_parallel_job_"*string(run_number)*".pbs"
    f = open(f_name, "w")
    write_pbs_settings(f)
    println(f, "/gpfs/group/ebf11/default/sw/julia-0.7.0/bin/julia compute_distances_given_params_parallel.jl "*string(run_number))
    close(f)

    return f_name
end

"""
Generates and submits "n_jobs" of PBS scripts to run "generate_pbs_compute_distances_given_params_parallel.jl".
"""
function submit_jobs_compute_distances_given_params_parallel(n_jobs::Int64)
    for i in 1:n_jobs
        f_name = generate_pbs_compute_distances_given_params_parallel(i)
        run(`qsub $f_name`) #this line submits the job by running 'qsub' in the command line!
        println("Job ", f_name, " submitted.")
    end
end

"""
Generates a PBS script for running "GP_draw_points_parallel.jl".
"""
function generate_pbs_GP_draw_points_parallel(run_number)
    f_name = "GP_draw_points_parallel_job_"*string(run_number)*".pbs"
    f = open(f_name, "w")
    write_pbs_settings(f)
    println(f, "/gpfs/group/ebf11/default/sw/julia-0.7.0/bin/julia GP_draw_points_parallel.jl "*string(run_number))
    close(f)

    return f_name
end

"""
Generates and submits "n_jobs" of PBS scripts to run "GP_draw_points_parallel.jl".
"""
function submit_jobs_GP_draw_points_parallel(n_jobs::Int64)
    for i in 1:n_jobs
        f_name = generate_pbs_GP_draw_points_parallel(i)
        run(`qsub $f_name`) #this line submits the job by running 'qsub' in the command line!
        println("Job ", f_name, " submitted.")
    end
end





##### NOTE: must run this script from the same directory in which you want to submit the jobs from!

##### To actually write a number of PBS scripts and submit them:

n_jobs = 20 # total number of jobs to submit

#submit_jobs_optimize(n_jobs)
#submit_jobs_compute_distances_given_params_random(n_jobs)
#submit_jobs_compute_distances_given_params_parallel(1)
#submit_jobs_GP_draw_points_parallel(1)
