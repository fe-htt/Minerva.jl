"""
Simulation orchestration functions for ensemble parameter studies.
Provides core functionality for running multiple simulations with different parameter combinations.
"""

"""
    save_sim_results(all_comb, plot_quantities, sol, sens_name, output_directory)

Save ensemble simulation results to CSV files.

# Arguments
- `all_comb`: all parameter combinations tested
- `plot_quantities`: variables to extract from solution
- `sol`: ensemble solution object
- `sens_name`: name of sensitivity study (creates subdirectory)
- `output_directory`: base output directory path
"""
function save_sim_results(all_comb, plot_quantities, sol, sim_name, output_directory)
  for (i,val) in enumerate(all_comb)
    local S = DataFrame()
    # Populate DataFrame with the simulation's results
    for mydata in plot_quantities
      insertcols!(S, Symbol(mydata) => sol[i][mydata])
    end
    # Rename the df column names to simpler names
    for n in names(S)
      rename!(S, n => replace(n, "(t)"=>"", "₊" => "."))
    end
    # Write to .csv file
    CSV.write(output_directory*"/"*sim_name*"/"*join(val,"-")*".csv", S)
  end
end


"""
    process_sim_settings(sim_values::Vector{Vector}, sim_params::Vector)

Processes simulation settings read from a Julia file defining sensitivity parameters.

The input file should define:
- `sim_values::Vector{Vector}`: values for each parameter
- `sim_params::Vector`: parameters to vary (ModelingToolkit parameters)

# Returns
Tuple of:
- `length(sim_params)`: number of parameters
- `lsr`: length of each parameter value vector
- `N_traj`: total number of simulations (product of all value counts)
- `all_comb`: iterator over all parameter combinations

# Example
```julia
# In sim_settings.jl:
#sim_name = "eta_study"
#sim_params = [sys.DSS.eta, sys.PI.etaT]
#sim_values = [[0.995, 0.999], [0.90, 0.95]]
```
"""
function process_sim_settings(sim_values, sim_params)
    lsr = [length(i) for i in sim_values]
    return length(sim_params), 
           lsr, 
           prod(lsr), 
           collect(IterTools.product(sim_values...))
end





"""
    prob_func(prob, i)

Create combination of parameters for ensemble of simulations.

# Arguments
- `prob`: base ODEProblem
- `i`: current parameter

# Returns
nothing
"""
function prob_func(prob, i, repeat)
  advance_msg(i)
  remake(prob, p = [sens_params[j] => all_comb[i][j] for j in eachindex(sens_params)]) # assigns to all sens parameters their value in the N-upla all_comb[i]
end



"""
    setup_ensemble_problem(prob, sens_params, all_comb, N_traj)

Create ensemble problem with progress reporting.

# Arguments
- `prob`: base ODEProblem
- `sens_params`: arameters to vary
- `all_comb`: all parameter combinations
- `N_traj`: total number of trajectories

# Returns
`EnsembleProblem` ready to solve
"""
function setup_ensemble_problem(prob, sim_params, all_comb, N_traj)
    # Progress tracking closure
    function advance_msg(i, col_widths, all_comb, N_traj)
        percent = (i / N_traj) * 100
        values = vcat([@sprintf("%d", i), @sprintf("%.2f", percent)], collect(all_comb[i]))
        row_str = join([@sprintf("%-*s", w, v) for (w, v) in zip(col_widths, values)], " | ")
        println(row_str)
    end
    
    # Problem function for ensemble
    function prob_func(prob, i, repeat)
        col_widths = [8, 12]
        append!(col_widths, fill(12, length(sim_params)))
        advance_msg(i, col_widths, all_comb, N_traj)
        remake(prob, p = [sim_params[j] => all_comb[i][j] for j in eachindex(sim_params)])
    end
    
    return EnsembleProblem(prob, prob_func = prob_func)
end


"""
    run_ensemble_sim(prob, sim_params, all_comb, N_traj; 
                    solver = Rodas5P(),
                    abstol = 1e-8,
                    maxiters = 1e6,
                    callback = nothing,
                    ensemble_alg = EnsembleSerial())

Run ensemble simulation with specified parameters.

# Arguments
- `prob`: base ODEProblem
- `sim_params`: parameters to vary
- `all_comb`: all parameter combinations
- `N_traj`: number of trajectories

# Keyword Arguments
- `solver`: ODE solver (default: Rodas5P())
- `abstol`: absolute tolerance
- `maxiters`: maximum iterations
- `callback`: optional callback (e.g., TerminateSteadyState)
- `ensemble_alg`: ensemble algorithm (default: EnsembleSerial())

# Returns
Ensemble solution object

# Example
```julia
cb_ss = TerminateSteadyState(1e-10, 1e-8)
sol = run_ensemble_sim(prob, sim_params, all_comb, N_traj;
                        callback = cb_ss,
                        abstol = 1e-10)
```
"""
function run_ensemble_sim(prob, sim_params, all_comb, N_traj; 
                        solver = Rodas5P(),
                        abstol = 1e-8,
                        maxiters = 1e6,
                        callback = nothing,
                        ensemble_alg = EnsembleSerial())
    
    ensemble_prob = setup_ensemble_problem(prob, sim_params, all_comb, N_traj)
    
    if callback !== nothing
        return solve(ensemble_prob, solver, ensemble_alg; 
                    trajectories = N_traj, 
                    abstol = abstol, 
                    maxiters = maxiters,
                    callback = callback)
    else
        return solve(ensemble_prob, solver, ensemble_alg; 
                    trajectories = N_traj, 
                    abstol = abstol, 
                    maxiters = maxiters)
    end
end


"""
    print_sim_header(sim_params)

Print formatted header for simulation progress tracking.

# Arguments
- `sim_params`: vector of parameter names/symbols

# Returns
Formatted header string with separators
"""
function print_sim_header(sim_params)
    col_widths = [8, 12]
    append!(col_widths, fill(12, length(sim_params)))
    header = vcat(["Sim #", "Progress %"], sim_params)
    separator = join(["="^w for w in col_widths], "===")
    header_str = join([@sprintf("%-*s", w, h) for (w, h) in zip(col_widths, header)], " | ")
    return separator * "\n" * header_str * "\n" * separator
end