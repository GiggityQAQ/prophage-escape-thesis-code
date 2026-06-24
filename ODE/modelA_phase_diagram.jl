using DifferentialEquations
using Plots
using Statistics
using Colors
using Plots.PlotMeasures

# =========================================================
# Model A phase diagram
#
# This script simulates Model A under periodic antibiotic
# pressure and generates one SPC dominant-attractor phase
# diagram.
#
# Model A assumption:
#   All prophages have the same baseline induction rate.
#   Antibiotic pressure affects sensitive cells through
#   antibiotic killing, but it does not increase prophage
#   induction.
#
# For each (T_on, T_off) parameter point, four initial-condition
# tests are performed:
#   1) Sresident
#   2) Cresident
#   3) President
#   4) DirectPC
#
# Each test is classified into one dominant mechanism:
#   S: sensitive mechanism
#   P: prophage-carried resistance mechanism
#   C: chromosome-carried resistance mechanism
#   N: extinction
#
# The final phase-diagram label for one parameter point is the
# set of unique outcomes observed across the four tests.
# =========================================================


# =========================================================
# User settings
# =========================================================

T_on_values  = 10.0:80.0:1200.0
T_off_values = 10.0:80.0:1200.0

# Use a smaller step size for a higher-resolution phase diagram.
# This will greatly increase the total simulation time.
# T_on_values  = 10.0:10.0:1200.0
# T_off_values = 10.0:10.0:1200.0

test_kinds = (
    :Sresident,
    :Cresident,
    :President,
    :DirectPC
)

rare_density = 1e-4

# Convergence is checked by comparing endpoint states across
# the last ncheck_cycles environmental cycles.
ncheck_cycles = 5
convergence_tol = 1e-6

# Resident establishment runs until convergence (early stop allowed).
max_resident_cycles = 300

# Post-invasion analysis runs a fixed number of complete antibiotic
# cycles instead of stopping early at convergence. Early stopping can
# freeze the system on a slow pre-extinction transient and over-report
# coexistence of the S/P/C mechanisms.
post_invasion_cycles = 200

# If total host abundance in the final analysis window is below
# this threshold, the outcome is classified as extinction.
host_extinction_threshold = 1e-10

# ODE solver output and tolerance settings.
saveat_value = 1.0
reltol_value = 1e-7
abstol_value = 1e-9

# Output file for the phase diagram.
output_filename = "modelA_phase_SPC_dominant_square.png"


# =========================================================
# Parameters
# =========================================================

Base.@kwdef mutable struct ParamsA
    lambda0::Float64 = 1.0
    alpha::Float64   = 1.0
    a_drug::Float64  = 1.0

    CR::Float64      = 0.05
    CRR::Float64     = 1.01 * 0.05

    phi::Float64     = 1e-2
    q::Float64       = 0.50

    eta::Float64     = 0.02
    B::Float64       = 50.0
    m::Float64       = 0.20
end


# =========================================================
# ODE system
#
# State variables:
#   1 S_S   : non-lysogenic ARG-negative host
#   2 S_R   : non-lysogenic ARG-positive host
#   3 L_SS  : ARG-negative host carrying ARG-negative prophage
#   4 L_SR  : ARG-negative host carrying ARG-positive prophage
#   5 L_RS  : ARG-positive host carrying ARG-negative prophage
#   6 L_RR  : ARG-positive host carrying ARG-positive prophage
#   7 V_S   : ARG-negative free phage
#   8 V_R   : ARG-positive free phage
#   9 A     : antibiotic state, 0 = absent, 1 = present
#
# Model A:
#   The same induction rate p.eta is used for all lysogens.
#   There is no antibiotic-induced increase in induction.
# =========================================================

function modelA!(du, u, p::ParamsA, t)
    S_S, S_R, L_SS, L_SR, L_RS, L_RR, V_S, V_R, A = u

    # Total bacterial host density controls density-dependent competition.
    N = S_S + S_R + L_SS + L_SR + L_RS + L_RR

    # Only non-lysogenic hosts are susceptible to infection in this model.
    S_total = S_S + S_R

    du[1] = S_S * (p.lambda0 - p.alpha * N - p.a_drug * A) -
            p.phi * V_S * S_S -
            p.phi * V_R * S_S

    du[2] = S_R * ((1.0 - p.CR) * p.lambda0 - p.alpha * N) -
            p.phi * V_S * S_R -
            p.phi * V_R * S_R

    du[3] = L_SS * (p.lambda0 - p.alpha * N - p.a_drug * A - p.eta) +
            p.q * p.phi * V_S * S_S

    du[4] = L_SR * ((1.0 - p.CR) * p.lambda0 - p.alpha * N - p.eta) +
            p.q * p.phi * V_R * S_S

    du[5] = L_RS * ((1.0 - p.CR) * p.lambda0 - p.alpha * N - p.eta) +
            p.q * p.phi * V_S * S_R

    du[6] = L_RR * ((1.0 - p.CRR) * p.lambda0 - p.alpha * N - p.eta) +
            p.q * p.phi * V_R * S_R

    du[7] = p.B * p.eta * (L_SS + L_RS) +
            p.B * (1.0 - p.q) * p.phi * V_S * S_total -
            p.phi * V_S * S_total -
            p.m * V_S

    du[8] = p.B * p.eta * (L_SR + L_RR) +
            p.B * (1.0 - p.q) * p.phi * V_R * S_total -
            p.phi * V_R * S_total -
            p.m * V_R

    du[9] = 0.0

    return nothing
end


# =========================================================
# Solver utilities
# =========================================================

function solve_segment(u0, t0, t1, Avalue, p::ParamsA;
                       saveat=1.0,
                       reltol=1e-7,
                       abstol=1e-9)

    # A is fixed within each segment because the antibiotic
    # environment is piecewise constant.
    u = copy(u0)
    u[9] = Avalue

    prob = ODEProblem(modelA!, u, (t0, t1), p)

    sol = solve(
        prob,
        Tsit5();
        saveat=saveat,
        reltol=reltol,
        abstol=abstol
    )

    return sol
end


function run_one_complete_cycle(u0, t0, T_off, T_on, p::ParamsA;
                                saveat=1.0,
                                reltol=1e-7,
                                abstol=1e-9)

    current_u = copy(u0)
    current_t = t0

    all_t = Float64[]
    all_X = Matrix{Float64}(undef, 0, 9)
    first_segment = true

    # First segment: antibiotic-free recovery phase.
    if T_off > 0
        sol_off = solve_segment(
            current_u,
            current_t,
            current_t + T_off,
            0.0,
            p;
            saveat=saveat,
            reltol=reltol,
            abstol=abstol
        )

        X_off = reduce(hcat, sol_off.u)'

        append!(all_t, sol_off.t)
        all_X = vcat(all_X, X_off)

        current_u = copy(sol_off.u[end])
        current_t = sol_off.t[end]
        first_segment = false
    end

    # Second segment: antibiotic exposure phase.
    if T_on > 0
        sol_on = solve_segment(
            current_u,
            current_t,
            current_t + T_on,
            1.0,
            p;
            saveat=saveat,
            reltol=reltol,
            abstol=abstol
        )

        X_on = reduce(hcat, sol_on.u)'

        # Avoid duplicating the boundary time point between OFF and ON segments.
        if first_segment
            append!(all_t, sol_on.t)
            all_X = vcat(all_X, X_on)
        else
            append!(all_t, sol_on.t[2:end])
            all_X = vcat(all_X, X_on[2:end, :])
        end

        current_u = copy(sol_on.u[end])
        current_t = sol_on.t[end]
    end

    return current_u, current_t, all_t, all_X
end


# =========================================================
# Convergence checking
# =========================================================

function max_host_endpoint_difference(u_prev::Vector{Float64},
                                      u_curr::Vector{Float64})

    # Convergence is evaluated using host-related variables only.
    return maximum(abs.(u_curr[1:6] .- u_prev[1:6]))
end


function has_converged_last_cycles(cycle_end_states::Vector{Vector{Float64}};
                                   ncheck=5,
                                   tol=1e-6)

    if length(cycle_end_states) < ncheck
        return false
    end

    last_states = cycle_end_states[end-ncheck+1:end]

    for k in 2:length(last_states)
        diff = max_host_endpoint_difference(last_states[k - 1], last_states[k])

        if diff >= tol
            return false
        end
    end

    return true
end


function run_to_convergence_with_history(u0, T_off, T_on, p::ParamsA;
                                         saveat=1.0,
                                         reltol=1e-7,
                                         abstol=1e-9,
                                         max_cycles=300,
                                         ncheck=5,
                                         tol=1e-6)

    current_u = copy(u0)
    current_t = 0.0

    all_t = Float64[]
    all_X = Matrix{Float64}(undef, 0, 9)
    cycle_end_states = Vector{Vector{Float64}}()

    first_cycle = true

    for cyc in 1:max_cycles
        current_u, current_t, cyc_t, cyc_X = run_one_complete_cycle(
            current_u,
            current_t,
            T_off,
            T_on,
            p;
            saveat=saveat,
            reltol=reltol,
            abstol=abstol
        )

        # Avoid duplicating the first time point of each new cycle.
        if first_cycle
            append!(all_t, cyc_t)
            all_X = vcat(all_X, cyc_X)
            first_cycle = false
        else
            append!(all_t, cyc_t[2:end])
            all_X = vcat(all_X, cyc_X[2:end, :])
        end

        push!(cycle_end_states, copy(current_u))

        if has_converged_last_cycles(
            cycle_end_states;
            ncheck=ncheck,
            tol=tol
        )
            return all_t, all_X, cycle_end_states, true, cyc
        end
    end

    return all_t, all_X, cycle_end_states, false, max_cycles
end


# =========================================================
# Fixed post-invasion simulation
#
# After invaders are introduced, the simulation always runs for
# exactly ncycles complete antibiotic cycles. Unlike the
# convergence-based routine, it never stops early, so slowly
# declining lineages are given time to reach true extinction.
# =========================================================

function run_fixed_post_invasion_cycles(u0, T_off, T_on, p::ParamsA;
                                        ncycles=200,
                                        saveat=1.0,
                                        reltol=1e-7,
                                        abstol=1e-9)

    current_u = copy(u0)
    current_t = 0.0

    all_t = Float64[]
    all_X = Matrix{Float64}(undef, 0, 9)
    cycle_end_states = Vector{Vector{Float64}}()

    first_cycle = true

    for cyc in 1:ncycles
        current_u, current_t, cyc_t, cyc_X = run_one_complete_cycle(
            current_u,
            current_t,
            T_off,
            T_on,
            p;
            saveat=saveat,
            reltol=reltol,
            abstol=abstol
        )

        # Avoid duplicating the first time point of each new cycle.
        if first_cycle
            append!(all_t, cyc_t)
            all_X = vcat(all_X, cyc_X)
            first_cycle = false
        else
            append!(all_t, cyc_t[2:end])
            all_X = vcat(all_X, cyc_X[2:end, :])
        end

        push!(cycle_end_states, copy(current_u))
    end

    return all_t, all_X, cycle_end_states
end


# =========================================================
# Initial conditions
#
# The four seed functions define the starting resident
# communities or direct-competition community before the
# invasion step.
# =========================================================

function seed_S_resident()
    return [
        0.1,
        0.0,
        0.1,
        0.0,
        0.0,
        0.0,
        1e-8,
        0.0,
        0.0
    ]
end


function seed_C_resident()
    return [
        0.0,
        0.1,
        0.0,
        0.0,
        0.1,
        0.0,
        1e-8,
        0.0,
        0.0
    ]
end


function seed_P_resident()
    return [
        0.1,
        0.0,
        0.0,
        0.1,
        0.0,
        0.0,
        0.0,
        1e-8,
        0.0
    ]
end


function seed_direct_PC()
    return [
        0.0,
        0.05,
        0.0,
        0.1,
        0.05,
        0.0,
        1e-8,
        1e-8,
        0.0
    ]
end


# =========================================================
# Build initial condition for each test
#
# Resident-based tests first run a resident community to a
# long-term state. Then rare invaders are added.
#
# DirectPC does not include a resident-establishment phase.
# =========================================================

function build_initial_condition(kind::Symbol, T_off, T_on, p::ParamsA;
                                 rare_density=1e-4,
                                 saveat=1.0,
                                 reltol=1e-7,
                                 abstol=1e-9,
                                 max_resident_cycles=300,
                                 ncheck=5,
                                 tol=1e-6)

    if kind == :Sresident
        _, _, cycle_states, converged, _ =
            run_to_convergence_with_history(
                seed_S_resident(),
                T_off,
                T_on,
                p;
                saveat=saveat,
                reltol=reltol,
                abstol=abstol,
                max_cycles=max_resident_cycles,
                ncheck=ncheck,
                tol=tol
            )

        if !converged
            return nothing, false
        end

        u0 = copy(cycle_states[end])

        # Rare P invader: prophage-carried ARG strategy.
        u0[4] += rare_density
        u0[8] += rare_density

        # Rare C invader: chromosome-carried ARG strategy.
        u0[2] += rare_density
        u0[5] += rare_density

        u0[9] = 0.0
        return u0, true

    elseif kind == :Cresident
        _, _, cycle_states, converged, _ =
            run_to_convergence_with_history(
                seed_C_resident(),
                T_off,
                T_on,
                p;
                saveat=saveat,
                reltol=reltol,
                abstol=abstol,
                max_cycles=max_resident_cycles,
                ncheck=ncheck,
                tol=tol
            )

        if !converged
            return nothing, false
        end

        u0 = copy(cycle_states[end])

        # Rare P invader: prophage-carried ARG strategy.
        u0[4] += rare_density
        u0[8] += rare_density

        u0[9] = 0.0
        return u0, true

    elseif kind == :President
        _, _, cycle_states, converged, _ =
            run_to_convergence_with_history(
                seed_P_resident(),
                T_off,
                T_on,
                p;
                saveat=saveat,
                reltol=reltol,
                abstol=abstol,
                max_cycles=max_resident_cycles,
                ncheck=ncheck,
                tol=tol
            )

        if !converged
            return nothing, false
        end

        u0 = copy(cycle_states[end])

        # Rare C invader: chromosome-carried ARG strategy.
        u0[2] += rare_density
        u0[5] += rare_density

        # A tiny V_S background allows ARG-negative phages to be present
        # after introducing chromosome-carried resistance.
        u0[7] += 1e-8

        u0[9] = 0.0
        return u0, true

    elseif kind == :DirectPC
        u0 = seed_direct_PC()
        u0[9] = 0.0
        return u0, true

    else
        error("Unknown test kind: $kind")
    end
end


# =========================================================
# Classification
#
# Classification is based on mean abundance over the final
# analysis window rather than a single endpoint.
#
# S mechanism:
#   S_S + L_SS
#
# P mechanism:
#   L_SR
#
# C mechanism:
#   S_R + L_RS + L_RR
#
# The dominant mechanism is the one with the largest mean
# abundance in the final window.
# =========================================================

function final_window_indices(t, T_off, T_on; n_cycles=5)
    cycle_len = T_off + T_on
    tmin = t[end] - n_cycles * cycle_len
    return findall(t .>= tmin)
end


function classify_one_run_modelA(t, X, T_off, T_on;
                                 host_extinction_threshold=1e-10,
                                 final_cycles=5)

    idx = final_window_indices(t, T_off, T_on; n_cycles=final_cycles)

    S_mech = mean(X[idx, 1] .+ X[idx, 3])
    P_mech = mean(X[idx, 4])
    C_mech = mean(X[idx, 2] .+ X[idx, 5] .+ X[idx, 6])

    total_host = S_mech + P_mech + C_mech

    if total_host < host_extinction_threshold
        return "N"
    end

    vals = [S_mech, P_mech, C_mech]
    labs = ["S", "P", "C"]

    return labs[argmax(vals)]
end


function combine_four_test_outcomes(outcomes::Vector{String})
    uniq = unique(outcomes)

    # If at least one non-extinct outcome exists, remove N from
    # the combined phase label.
    if length(uniq) > 1
        uniq = filter(x -> x != "N", uniq)
    end

    order = Dict(
        "N" => 1,
        "S" => 2,
        "P" => 3,
        "C" => 4
    )

    uniq = sort(uniq, by=x -> order[x])

    return join(uniq, ",")
end


# =========================================================
# One parameter point
#
# For one (T_on, T_off) pair:
#   1. Build initial condition for each test.
#   2. Run post-invasion dynamics.
#   3. Classify the final dominant mechanism.
#   4. Combine unique outcomes across the four tests.
# =========================================================

function classify_point_modelA(T_on, T_off, p::ParamsA;
                               rare_density=1e-4,
                               saveat=1.0,
                               reltol=1e-7,
                               abstol=1e-9,
                               max_resident_cycles=300,
                               post_invasion_cycles=200,
                               ncheck=5,
                               convergence_tol=1e-6,
                               host_extinction_threshold=1e-10)

    outcomes = String[]

    for kind in test_kinds
        u0, ok = build_initial_condition(
            kind,
            T_off,
            T_on,
            p;
            rare_density=rare_density,
            saveat=saveat,
            reltol=reltol,
            abstol=abstol,
            max_resident_cycles=max_resident_cycles,
            ncheck=ncheck,
            tol=convergence_tol
        )

        if !ok
            push!(outcomes, "N")
            continue
        end

        # Fixed-length post-invasion run (no early stop). This lets
        # slowly declining lineages reach true extinction instead of
        # being frozen on a transient.
        t, X, _ =
            run_fixed_post_invasion_cycles(
                u0,
                T_off,
                T_on,
                p;
                ncycles=post_invasion_cycles,
                saveat=saveat,
                reltol=reltol,
                abstol=abstol
            )

        outcome = classify_one_run_modelA(
            t,
            X,
            T_off,
            T_on;
            host_extinction_threshold=host_extinction_threshold,
            final_cycles=ncheck
        )

        push!(outcomes, outcome)
    end

    return combine_four_test_outcomes(outcomes)
end


# =========================================================
# Phase scan
#
# This function scans all combinations of T_on and T_off.
# Rows correspond to T_off values, and columns correspond to
# T_on values.
# =========================================================

function scan_phase_modelA(p::ParamsA;
                           T_on_values=10.0:80.0:1200.0,
                           T_off_values=10.0:80.0:1200.0,
                           rare_density=1e-4,
                           saveat=1.0,
                           reltol=1e-7,
                           abstol=1e-9,
                           max_resident_cycles=300,
                           post_invasion_cycles=200,
                           ncheck=5,
                           convergence_tol=1e-6,
                           host_extinction_threshold=1e-10)

    labels = Matrix{String}(undef, length(T_off_values), length(T_on_values))

    for i in eachindex(T_off_values)
        println("Row ", i, "/", length(T_off_values), " | T_off = ", T_off_values[i])

        for j in eachindex(T_on_values)
            labels[i, j] = classify_point_modelA(
                T_on_values[j],
                T_off_values[i],
                p;
                rare_density=rare_density,
                saveat=saveat,
                reltol=reltol,
                abstol=abstol,
                max_resident_cycles=max_resident_cycles,
                post_invasion_cycles=post_invasion_cycles,
                ncheck=ncheck,
                convergence_tol=convergence_tol,
                host_extinction_threshold=host_extinction_threshold
            )
        end
    end

    return collect(T_on_values), collect(T_off_values), labels
end


# =========================================================
# Plot helper
# =========================================================

function phase_axis_limits(T_on_values, T_off_values)
    axis_upper = maximum([
        1200.0,
        maximum(T_on_values),
        maximum(T_off_values)
    ])

    return (0.0, axis_upper)
end


# =========================================================
# Plotting
#
# The phase diagram uses categorical colors for combined
# SPC outcomes. The plot is forced to be square so that
# T_on and T_off have the same visual scale.
# =========================================================

function plot_phase_modelA(T_on_values, T_off_values, labels;
                           filename="modelA_phase_SPC_dominant_square.png")

    color_dict = Dict(
        "N"     => RGB(0.85, 0.85, 0.85),
        "S"     => RGB(0.62, 0.78, 0.38),
        "P"     => RGB(0.20, 0.47, 0.78),
        "C"     => RGB(0.86, 0.55, 0.12),
        "S,P"   => RGB(0.00, 0.75, 0.85),
        "S,C"   => RGB(0.78, 0.72, 0.18),
        "P,C"   => RGB(0.00, 0.48, 0.22),
        "S,P,C" => RGB(0.48, 0.18, 0.65)
    )

    order = [
        "N",
        "S",
        "P",
        "C",
        "S,P",
        "S,C",
        "P,C",
        "S,P,C"
    ]

    present = unique(vec(labels))
    cats = [c for c in order if c in present]

    cat_to_int = Dict(c => i for (i, c) in enumerate(cats))

    Z = [
        cat_to_int[labels[i, j]]
        for i in axes(labels, 1), j in axes(labels, 2)
    ]

    cols = [color_dict[c] for c in cats]
    cmap = cgrad(cols, categorical=true)

    lims = phase_axis_limits(T_on_values, T_off_values)

    plt = heatmap(
        T_on_values,
        T_off_values,
        Z;
        color=cmap,
        clims=(0.5, length(cats) + 0.5),
        colorbar=false,
        xlabel="Antibiotic Presence Period (T_on)",
        ylabel="Antibiotic Absence Period (T_off)",
        xticks=0:300:lims[2],
        yticks=0:300:lims[2],
        xlims=lims,
        ylims=lims,
        aspect_ratio=:equal,
        framestyle=:box,
        grid=false,
        size=(1800, 1800),
        dpi=300,
        title="Model A Phase Diagram\nSPC Dominant Attractor Set",
        titlefontsize=14,
        guidefontsize=14,
        tickfontsize=11,
        legendfontsize=9,
        legend=:outerright,
        left_margin=20mm,
        bottom_margin=24mm,
        right_margin=42mm,
        top_margin=16mm
    )

    # Add invisible dummy points to create a categorical legend.
    for cat in cats
        scatter!(
            plt,
            [NaN],
            [NaN];
            label=cat,
            marker=:square,
            markersize=10,
            markercolor=color_dict[cat],
            markerstrokecolor=color_dict[cat]
        )
    end

    # Draw thin vertical grid lines to show phase-grid boundaries.
    for x in T_on_values
        vline!(
            plt,
            [x];
            color=:white,
            linewidth=0.25,
            label=false,
            alpha=0.30
        )
    end

    # Draw thin horizontal grid lines to show phase-grid boundaries.
    for y in T_off_values
        hline!(
            plt,
            [y];
            color=:white,
            linewidth=0.25,
            label=false,
            alpha=0.30
        )
    end

    # Add text labels only for regions large enough to avoid clutter.
    min_cells_for_label = 4

    for cat in cats
        idxs = findall(x -> x == cat, labels)

        if length(idxs) >= min_cells_for_label
            xs = [T_on_values[Tuple(I)[2]] for I in idxs]
            ys = [T_off_values[Tuple(I)[1]] for I in idxs]

            annotate!(
                plt,
                mean(xs),
                mean(ys),
                text(cat, 13, :black)
            )
        end
    end

    savefig(plt, filename)
    println("Saved phase diagram: $filename")

    return plt
end


# =========================================================
# Main
#
# Run the full phase scan, generate the phase diagram, display
# the plot, and return the plot object.
# =========================================================

function main()
    p = ParamsA()

    T_on_vals, T_off_vals, labels = scan_phase_modelA(
        p;
        T_on_values=T_on_values,
        T_off_values=T_off_values,
        rare_density=rare_density,
        saveat=saveat_value,
        reltol=reltol_value,
        abstol=abstol_value,
        max_resident_cycles=max_resident_cycles,
        post_invasion_cycles=post_invasion_cycles,
        ncheck=ncheck_cycles,
        convergence_tol=convergence_tol,
        host_extinction_threshold=host_extinction_threshold
    )

    plt = plot_phase_modelA(
        T_on_vals,
        T_off_vals,
        labels;
        filename=output_filename
    )

    display(plt)

    return plt
end


plt = main()