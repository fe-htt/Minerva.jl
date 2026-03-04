using Printf
using DifferentialEquations
using Minerva

# Get current file name and current directory
current_file = @__FILE__
current_directory = @__DIR__

# Define the base for the output files name
file_name = splitext(basename(current_directory))[1]

# Input file
include(joinpath(current_directory, "input.jl"))

# Read simulation settings and process input settings
include(joinpath(current_directory, "sim_settings.jl"))
ls, lsr, N_traj, all_comb = process_sim_settings(sim_values, sim_params)
println("# N. OF RUNS: ", N_traj)

# Run_information
println(print_sim_header(sim_params))


function prob_func(prob, i, repeat)
  advance_msg(i)
  remake(prob, p = [sim_params[j] => all_comb[i][j] for j in eachindex(sim_params)]) # assigns to all sens parameters their value in the N-upla all_comb[i]
end

# Print progress message for each completed run of the ensemble of simulations
col_widths = [8, 12];
append!(col_widths, fill(12, length(sim_params)));  # Adjust for parameters
function advance_msg(i)
  percent = (i / N_traj) * 100
  values = vcat([@sprintf("%d", i), @sprintf("%.2f", percent)], collect(all_comb[i]))
  row_str = join([@sprintf("%-*s", w, v) for (w, v) in zip(col_widths, values)], " | ")
  println(row_str)
  println(values)
end

# ============================================================================ #
#tag RUN ensemble
# ============================================================================ #
cb_ss = TerminateSteadyState(1e-10, 1e-8; min_t = nothing) # abstol, reltol
global prob = ODEProblem(sys, u0, tspan, jac=true) # jac=true, autodiff=true ---- ODEProblem(sys, u0, tspan, [], jac=true)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
global sol = solve(ensemble_prob, Rodas5P(), EnsembleSerial(); trajectories=N_traj, abstol=1e-8, maxiters=1e6) # progress=true, Rodas5(), Rosenbrock23(), saveat=5.0, maxiters=1e7, abstol=1e-15, reltol=1e-20
#global sol = solve(ensemble_prob, Rodas5P(), EnsembleSerial(); trajectories=N_traj, abstol=1e-8, maxiters=1e6, callback=cb_ss) # progress=true, Rodas5(), Rosenbrock23(), saveat=5.0, maxiters=1e7, abstol=1e-15, reltol=1e-20
#summ = EnsembleSummary(sol)#, 0:1000.0:tspan[end])





# Write output file for each run
print("# Saving output to file ...")
output_path = joinpath(current_directory, sim_name)
mkpath(output_path)
save_sim_results(all_comb, plot_quantities, sol, sim_name, current_directory)
# Save ensamble summary as JLD2
# @save current_directory*"/"*sim_name*"/ensemble_summary_1.jld2" summ  # Save
println(" done.")



