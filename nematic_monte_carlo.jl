using Plots
# %%
mutable struct LiquidCrystalSystem
    x::Vector{Float64}
    y::Vector{Float64}
    theta::Vector{Float64}
end
# %%
function total_energy(system::LiquidCrystalSystem, H::Float64, molecule_length::Float64)
    """
        Compute the total energy of the liquid crystal system

        INPUT
        - system: current state of the liquid crystal system
        - H: external field strength
        - molecule_length: length of each liquid crystal molecule

        OUTPUT
        - E: total energy of the system
    """
    E = 0.0
    num = length(system.x)
    for i in 1:num
        for j in 1:num
            if system.theta[i] != system.theta[j]
                factor1 = ((system.x[j] - system.x[i]) * sin(system.theta[j]) - (system.y[j] - system.y[i]) * cos(system.theta[j])) / sin(system.theta[j] - system.theta[i]) / molecule_length
                factor2 = ((system.x[i] - system.x[j]) * sin(system.theta[i]) - (system.y[i] - system.y[j]) * cos(system.theta[i])) / sin(system.theta[i] - system.theta[j]) / molecule_length
                if 0 <= factor1 <= 1 && 0 <= factor2 <= 1
                    E += 1
                end
            end
        end
        E -= (H * sin(system.theta[i]))^2
    end
    return E
end
# %%
function boltzmann_factor(E::Float64, T::Float64)
    """
        Compute the Boltzmann factor for given energy and temperature

        INPUT
        - E: energy difference
        - T: temperature

        OUTPUT
        - factor: Boltzmann factor exp(-E / kT), with k ≈ 1e-23
    """
    return exp(-1.0e23 * E / T)
end
# %%
function propose_move(system::LiquidCrystalSystem, H::Float64, T::Float64, N::Int, system_size::Float64, molecule_length::Float64)
    """
        Propose a random move for a randomly selected molecule

        INPUT
        - system: current state of the system
        - H: external field strength
        - T: temperature
        - N: number of molecules
        - system_size: size of the square system (length of one side)
        - molecule_length: length of each molecule

        OUTPUT
        - new_system: new state of the system after proposed move (may be same as before if move is invalid)
    """
    i = rand(1:N)
    old_x, old_y, old_theta = system.x[i], system.y[i], system.theta[i]
    new_system = LiquidCrystalSystem(copy(system.x), copy(system.y), copy(system.theta))
    new_system.x[i] += 2 * rand() - 1
    new_system.y[i] += 2 * rand() - 1
    new_system.theta[i] += (π/2) * (rand() - 0.5)

    x_end = new_system.x[i] + molecule_length * cos(new_system.theta[i])
    y_end = new_system.y[i] + molecule_length * sin(new_system.theta[i])

    if new_system.x[i] < 0 || new_system.x[i] > system_size || new_system.y[i] < 0 || new_system.y[i] > system_size || x_end < 0 || x_end > system_size || y_end < 0 || y_end > system_size
        new_system.x[i], new_system.y[i], new_system.theta[i] = old_x, old_y, old_theta
    end

    return new_system
end
# %%
function metropolis_step(system::LiquidCrystalSystem, H::Float64, T::Float64, N::Int, system_size::Float64, molecule_length::Float64)
    """
        Perform a single Metropolis update step

        INPUT
        - system: current state of the system
        - H: external field strength
        - T: temperature
        - N: number of molecules
        - system_size: size of the system
        - molecule_length: length of each molecule

        OUTPUT
        - new_system: updated system state after Metropolis step
    """
    E_old = total_energy(system, H, molecule_length)
    new_system = propose_move(system, H, T, N, system_size, molecule_length)
    E_new = total_energy(new_system, H, molecule_length)

    if E_new <= E_old || rand() <= boltzmann_factor(E_new - E_old, T)
        return new_system
    else
        return system
    end
end
# %%
function monte_carlo_sweep(system::LiquidCrystalSystem, H::Float64, T::Float64, N::Int, system_size::Float64, molecule_length::Float64)
    """
        Perform one full Monte Carlo sweep (N Metropolis steps)

        INPUT
        - system: initial system state
        - H: external field strength
        - T: temperature
        - N: number of molecules
        - system_size: size of the system
        - molecule_length: length of each molecule

        OUTPUT
        - result: system state after one sweep
    """
    result = LiquidCrystalSystem(copy(system.x), copy(system.y), copy(system.theta))
    for _ in 1:N
        result = metropolis_step(result, H, T, N, system_size, molecule_length)
    end
    return result
end
# %%
function run_simulation(H::Float64, T::Float64, N::Int, num_steps::Int, system_size::Float64, molecule_length::Float64)
    """
        Run the full Monte Carlo simulation for the liquid crystal system

        INPUT
        - H: external field strength
        - T: temperature
        - N: number of molecules
        - num_steps: number of Monte Carlo sweeps
        - system_size: size of the system
        - molecule_length: length of each molecule

        OUTPUT
        - results: array of LiquidCrystalSystem objects representing the system at each time step
    """
    results = Vector{LiquidCrystalSystem}(undef, num_steps + 1)
    results[1] = LiquidCrystalSystem([system_size * rand() for _ in 1:N], [system_size * rand() for _ in 1:N], [π/2 for _ in 1:N])

    for i in 1:N
        while results[1].x[i] + molecule_length * cos(results[1].theta[i]) < 0 ||
              results[1].x[i] + molecule_length * cos(results[1].theta[i]) > system_size ||
              results[1].y[i] + molecule_length * sin(results[1].theta[i]) < 0 ||
              results[1].y[i] + molecule_length * sin(results[1].theta[i]) > system_size
            results[1].x[i] = system_size * rand()
            results[1].y[i] = system_size * rand()
        end
    end

    for i in 1:num_steps
        results[i+1] = monte_carlo_sweep(results[i], H, T, N, system_size, molecule_length)
    end
    return results
end
# %%
magnetic_field = 0.0
temperature = 300.
num_particles = 100
num_steps = 200
system_size = 100.0
molecule_length = 15.0
# %%
@time simulation_results = run_simulation(0.0, 300.0, num_particles, num_steps, system_size, molecule_length);
open("./LC.data.H$(magnetic_field).T$(temperature).num$(num_particles).time$(num_steps).txt", "w") do io
    println(io, simulation_results)
end
# %%
anim = @animate for i in 1:(num_steps + 1)
    t = i - 1
    p = quiver(
        simulation_results[i].x,
        simulation_results[i].y,
        quiver = (molecule_length * cos.(simulation_results[i].theta), molecule_length * sin.(simulation_results[i].theta)),
        xlim = (0, system_size), ylim = (0, system_size),
        aspect_ratio = 1, color = 1,
        title = "Liquid Crystal (N=$(num_particles), H=0) at t=$t"
    )
    plot(p)
end
# %%
gif(anim, "./anim.H$(magnetic_field).T$(temperature).num$(num_particles).time$(num_steps).gif", fps = 10)
# %%
function order_parameter(data, N)
    """
        Compute the nematic order parameter over time

        INPUT
        - data: array of system states (LiquidCrystalSystem)
        - N: number of molecules

        OUTPUT
        - result: vector of order parameters at each time step
    """
    result = zeros(length(data))
    for i in 1:length(data)
        sum_cos = sum(x -> cos(2*x), data[i].theta)
        sum_sin = sum(x -> sin(2*x), data[i].theta)
        mean_cos = sum_cos / N
        mean_sin = sum_sin / N
        result[i] = sqrt(mean_cos^2 + mean_sin^2)
    end
    return result
end
# %%
Ss = order_parameter(simulation_results, num_particles)
# %%
figS = plot(
    0:num_steps, Ss, 
    xlabel="Time",
    ylabel="Order parameter",
    title="N=$(num_particles), H=$(magnetic_field)",
    legend=nothing,
    ylims=(0.,1.)
)
savefig(figS, "./S.H$(magnetic_field).T$(temperature).num$(num_particles).time$(num_steps).png")
# %%
# %%
Ns = [10,50,100,150,200,250,300,400]
order_parameters_dict = Dict{Int, Vector{Float64}}()
for N in Ns
    println("Running simulation for N = $N")
    result = run_simulation(magnetic_field, temperature, N, num_steps, system_size, molecule_length)
    S_series = order_parameter(result, N)
    order_parameters_dict[N] = S_series
end
# %%
fig_comparison = plot(
    xlabel="Time", ylabel="Order Parameter", title="Comparison for Various N"
)
for N in Ns
    plot!(fig_comparison, 0:num_steps, order_parameters_dict[N], label="N=$N")
end
# %%
savefig(fig_comparison, "./S_comparison.H$(magnetic_field).T$(temperature).time$(num_steps).png")
