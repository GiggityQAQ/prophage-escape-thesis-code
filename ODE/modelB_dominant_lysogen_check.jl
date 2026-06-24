using DifferentialEquations
using Plots
using Colors
using Plots.PlotMeasures
using Statistics

# =========================================================
# Model B quick-test phase diagram
#
# This script generates a simplified Model B phase diagram.
#
# Model B assumption:
#   ARG-negative prophages have antibiotic-induced extra
#   induction during antibiotic exposure:
#
#       eta_S = eta + epsilon * A
#
#   ARG-positive prophages keep the baseline induction rate:
#
#       eta_R = eta
#
# Quick-test classification rule:
#   For each single simulation run, only one lysogen type is
#   retained as the final outcome.
#
#   The dominant lysogen is defined as the lysogen population
#   with the largest maximum density over the final analysis
#   window.
#
# Compared with the lysogen-presence phase diagram:
#   - This script does not report within-run coexistence.
#   - Comma "," labels are not used.
#   - Different outcomes across initial-condition tests are
#     still combined using " | ".
#
# For each (T_on, T_off) parameter point, four initial-condition
# tests are performed:
#   1) Sresident
#   2) Cresident
#   3) President
#   4) DirectPC
#
# Example:
#   Sresident  -> L_SR
#   Cresident  -> L_SR
#   President  -> L_SR
#   DirectPC   -> L_SS
#
# Final grid label:
#   L_SS | L_SR
#
# Only one phase diagram image is generated.
# =========================================================


# =========================================================
# User settings
# =========================================================

# Parameter grid for the quick test.
T_on_values  = 10.0:80.0:1200.0
T_off_values = 10.0:80.0:1200.0

# Use this smaller grid for faster testing.
# T_on_values  = [10.0, 300.0, 600.0, 900.0, 1200.0]
# T_off_values = [10.0, 300.0, 600.0, 900.0, 1200.0]

# Use this finer grid for a higher-resolution final figure.
# This will substantially increase simulation time.
# T_on_values  = 10.0:10.0:1200.0
# T_off_values = 10.0:10.0:1200.0

test_kinds = (
    :Sresident,
    :Cresident,
    :President,
    :DirectPC
)

rare_density = 1e-4

# Convergence is evaluated by comparing endpoint states across
# the last ncheck_cycles environmental cycles.
ncheck_cycles = 5
convergence_tol = 1e-6

# Resident establishment runs until convergence (early stop allowed).
max_resident_cycles = 300

# Post-invasion analysis runs a fixed number of complete antibiotic
# cycles (matching the time-series script) instead of stopping early at
# convergence. Early stopping can freeze the system on a slow
# pre-extinction transient and over-report lysogen coexistence.
post_invasion_cycles = 200

# If the largest lysogen density in the final analysis window
# is below this threshold, the outcome is classified as N.
lysogen_extinction_threshold = 1e-10

# ODE solver output interval and tolerance settings.
saveat_value = 1.0
reltol_value = 1e-7
abstol_value = 1e-9

# Output file for the phase diagram.
output_filename = "modelB_phase_dominant_lysogen_type_by_test_square.png"


# =========================================================
# Parameters
# =========================================================

Base.@kwdef mutable struct ParamsB
    lambda0::Float64 = 1.0
    alpha::Float64   = 1.0
    a_drug::Float64  = 1.0

    CR::Float64      = 0.05
    CRR::Float64     = 1.01 * 0.05

    phi::Float64     = 1e-2
    q::Float64       = 0.50

    eta::Float64     = 0.02
    epsilon::Float64 = 0.02

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
# In Model B, antibiotic-induced extra induction applies only
# to ARG-negative prophages.
# =========================================================

function modelB!(du, u, p::ParamsB, t)
    S_S, S_R, L_SS, L_SR, L_RS, L_RR, V_S, V_R, A = u

    # Total bacterial host density controls density-dependent competition.
    N = S_S + S_R + L_SS + L_SR + L_RS + L_RR

    # Only non-lysogenic hosts are available for new phage infection.
    S_total = S_S + S_R

    # Model B induction rates.
    # ARG-negative prophages receive extra induction when A = 1.
    eta_S = p.eta + p.epsilon * A
    eta_R = p.eta

    du[1] = S_S * (p.lambda0 - p.alpha * N - p.a_drug * A) -
            p.phi * V_S * S_S -
            p.phi * V_R * S_S

    du[2] = S_R * ((1.0 - p.CR) * p.lambda0 - p.alpha * N) -
            p.phi * V_S * S_R -
            p.phi * V_R * S_R

    du[3] = L_SS * (p.lambda0 - p.alpha * N - p.a_drug * A - eta_S) +
            p.q * p.phi * V_S * S_S

    du[4] = L_SR * ((1.0 - p.CR) * p.lambda0 - p.alpha * N - eta_R) +
            p.q * p.phi * V_R * S_S

    du[5] = L_RS * ((1.0 - p.CR) * p.lambda0 - p.alpha * N - eta_S) +
            p.q * p.phi * V_S * S_R

    du[6] = L_RR * ((1.0 - p.CRR) * p.lambda0 - p.alpha * N - eta_R) +
            p.q * p.phi * V_R * S_R

    du[7] = p.B * eta_S * (L_SS + L_RS) +
            p.B * (1.0 - p.q) * p.phi * V_S * S_total -
            p.phi * V_S * S_total -
            p.m * V_S

    du[8] = p.B * eta_R * (L_SR + L_RR) +
            p.B * (1.0 - p.q) * p.phi * V_R * S_total -
            p.phi * V_R * S_total -
            p.m * V_R

    du[9] = 0.0

    return nothing
end


# =========================================================
# Solver utilities
# =========================================================

function solve_segment(u0, t0, t1, Avalue, p::ParamsB;
                       saveat=1.0,
                       reltol=1e-7,
                       abstol=1e-9)

    # The antibiotic state is fixed within one segment.
    # Avalue = 0.0 for the antibiotic-free phase.
    # Avalue = 1.0 for the antibiotic-present phase.
    u = copy(u0)
    u[9] = Avalue

    prob = ODEProblem(modelB!, u, (t0, t1), p)

    sol = solve(
        prob,
        Tsit5();
        saveat=saveat,
        reltol=reltol,
        abstol=abstol
    )

    return sol
end


function run_one_complete_cycle(u0, t0, T_off, T_on, p::ParamsB;
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

    # Convergence is evaluated using the six host-related variables.
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


function run_to_convergence_with_history(u0, T_off, T_on, p::ParamsB;
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
        current_u, current_t, cyc_t, cyc_X =
            run_one_complete_cycle(
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
# This matches the protocol used by the time-series script.
# =========================================================

function run_fixed_post_invasion_cycles(u0, T_off, T_on, p::ParamsB;
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
# communities or the direct-competition community.
#
# The order of entries follows the state-variable order used
# by modelB!.
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
# Resident-based tests:
#   1. Run the resident community to convergence.
#   2. Add rare invaders.
#
# DirectPC:
#   Starts directly from the mixed initial condition.
# =========================================================

function build_initial_condition(kind::Symbol, T_off, T_on, p::ParamsB;
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

        # Tiny V_S background after introducing chromosome-carried resistance.
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
# Dominant lysogen classification within one run
#
# The final analysis window covers the last final_cycles
# complete antibiotic cycles.
#
# For each lysogen type:
#   L_SS, L_SR, L_RS, L_RR
#
# the maximum density in the final window is calculated.
#
# The single lysogen type with the largest final-window maximum
# is returned as the outcome of the run.
#
# If all lysogens remain below the extinction threshold, the
# run is classified as N.
# =========================================================

function final_window_indices(t, T_off, T_on; n_cycles=5)
    cycle_len = T_off + T_on
    tmin = t[end] - n_cycles * cycle_len
    return findall(t .>= tmin)
end


function classify_one_run_dominant_lysogen(t, X, T_off, T_on;
                                           lysogen_extinction_threshold=1e-10,
                                           final_cycles=5)

    idx = final_window_indices(t, T_off, T_on; n_cycles=final_cycles)

    labels = ["L_SS", "L_SR", "L_RS", "L_RR"]

    scores = [
        maximum(X[idx, 3]),  # L_SS
        maximum(X[idx, 4]),  # L_SR
        maximum(X[idx, 5]),  # L_RS
        maximum(X[idx, 6])   # L_RR
    ]

    best_idx = argmax(scores)
    best_score = scores[best_idx]

    if best_score <= lysogen_extinction_threshold
        return "N"
    end

    return labels[best_idx]
end


# =========================================================
# Across-test outcome combination
#
# Each single run returns only one dominant lysogen type.
#
# Across the four initial-condition tests, different dominant
# outcomes are combined using " | ".
#
# Example:
#   Sresident  -> L_SR
#   Cresident  -> L_SR
#   President  -> L_SR
#   DirectPC   -> L_SS
#
# Combined outcome:
#   L_SS | L_SR
# =========================================================

function single_label_rank(label::String)
    order = Dict(
        "N"    => 0,
        "L_SS" => 1,
        "L_SR" => 2,
        "L_RS" => 3,
        "L_RR" => 4
    )

    return get(order, label, 99)
end


function combined_label_sort_key(label::String)

    # Combined labels are sorted first by the number of distinct
    # across-test outcomes, then by lysogen rank.
    parts = split(label, " | ")
    ranks = [single_label_rank(String(p)) for p in parts]

    return (length(parts), sum(ranks), label)
end


function combine_four_test_dominant_outcomes(outcomes::Vector{String})
    uniq = unique(outcomes)

    # If at least one non-extinct outcome exists, remove N from
    # the combined phase label.
    if length(uniq) > 1
        uniq = filter(x -> x != "N", uniq)
    end

    uniq = sort(uniq, by=single_label_rank)

    return join(uniq, " | ")
end


# =========================================================
# One parameter point
#
# For one (T_on, T_off) pair:
#   1. Run all four initial-condition tests.
#   2. Classify the dominant lysogen type for each test.
#   3. Combine distinct outcomes across the four tests.
# =========================================================

function classify_point_modelB_dominant_by_test(T_on, T_off, p::ParamsB;
                                                rare_density=1e-4,
                                                saveat=1.0,
                                                reltol=1e-7,
                                                abstol=1e-9,
                                                max_resident_cycles=300,
                                                post_invasion_cycles=200,
                                                ncheck=5,
                                                convergence_tol=1e-6,
                                                lysogen_extinction_threshold=1e-10)

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
        # being frozen on a transient, matching the time-series script.
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

        outcome = classify_one_run_dominant_lysogen(
            t,
            X,
            T_off,
            T_on;
            lysogen_extinction_threshold=lysogen_extinction_threshold,
            final_cycles=ncheck
        )

        push!(outcomes, outcome)
    end

    return combine_four_test_dominant_outcomes(outcomes)
end


# =========================================================
# Phase scan
#
# This function scans all combinations of T_on and T_off.
# Rows correspond to T_off values, and columns correspond to
# T_on values.
# =========================================================

function scan_phase_modelB_dominant_by_test(p::ParamsB;
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
                                            lysogen_extinction_threshold=1e-10)

    labels = Matrix{String}(undef, length(T_off_values), length(T_on_values))

    for i in eachindex(T_off_values)
        println("Row ", i, "/", length(T_off_values), " | T_off = ", T_off_values[i])

        for j in eachindex(T_on_values)
            labels[i, j] = classify_point_modelB_dominant_by_test(
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
                lysogen_extinction_threshold=lysogen_extinction_threshold
            )
        end
    end

    return collect(T_on_values), collect(T_off_values), labels
end


# =========================================================
# Plot helpers
# =========================================================

function half_grid_step(vals)

    # Estimate half of the grid spacing so heatmap cells are
    # centered correctly on the provided T_on and T_off values.
    if length(vals) <= 1
        return 0.5
    end

    sorted_vals = sort(collect(vals))
    return minimum(diff(sorted_vals)) / 2
end


function square_axis_limits(xvals, yvals)

    # Use identical x- and y-axis limits to keep the phase diagram square.
    x_half = half_grid_step(xvals)
    y_half = half_grid_step(yvals)

    x_min = minimum(xvals) - x_half
    x_max = maximum(xvals) + x_half

    y_min = minimum(yvals) - y_half
    y_max = maximum(yvals) + y_half

    lim_min = min(x_min, y_min)
    lim_max = max(x_max, y_max)

    return (lim_min, lim_max)
end


# =========================================================
# Plotting
#
# The phase diagram uses dynamic categorical colors because
# the observed combined labels depend on the simulation result.
#
# The plot is exported as a square figure:
#   - size=(1800, 1800) creates a square canvas.
#   - aspect_ratio=:equal keeps x and y units visually equal.
#   - xlims and ylims share the same limits.
#   - margins reserve space for labels, title, and legend.
# =========================================================

function plot_phase_modelB_dominant_by_test(T_on_values, T_off_values, labels;
                                            filename="modelB_phase_dominant_lysogen_type_by_test_square.png")

    present = unique(vec(labels))

    # Sort labels by the number of across-test outcomes and lysogen rank.
    cats = sort(present, by=combined_label_sort_key)

    cat_to_int = Dict(c => i for (i, c) in enumerate(cats))

    Z = [
        cat_to_int[labels[i, j]]
        for i in axes(labels, 1), j in axes(labels, 2)
    ]

    # Generate categorical colors automatically from observed labels.
    ncat = length(cats)

    if ncat <= 1
        cols = [RGB(0.30, 0.30, 0.30)]
    else
        cols = distinguishable_colors(
            ncat,
            [RGB(1, 1, 1), RGB(0, 0, 0)],
            dropseed=true
        )
    end

    color_dict = Dict(cats[i] => cols[i] for i in eachindex(cats))
    cmap = cgrad(cols, categorical=true)

    lims = square_axis_limits(T_on_values, T_off_values)

    plt = heatmap(
        T_on_values,
        T_off_values,
        Z;
        color=cmap,
        clims=(0.5, ncat + 0.5),
        colorbar=false,
        xlabel="Antibiotic Presence Period (T_on)",
        ylabel="Antibiotic Absence Period (T_off)",
        xticks=0:300:1200,
        yticks=0:300:1200,
        xlims=lims,
        ylims=lims,
        aspect_ratio=:equal,
        framestyle=:box,
        grid=false,
        size=(1800, 1800),
        dpi=300,
        title="Model B Phase Diagram\nDominant Lysogen Type by Initial Condition",
        titlefontsize=13,
        guidefontsize=13,
        tickfontsize=10,
        legendfontsize=7,
        legend=:outerright,
        left_margin=18mm,
        bottom_margin=22mm,
        right_margin=80mm,
        top_margin=14mm
    )

    # Add invisible dummy points to create a categorical legend.
    for cat in cats
        scatter!(
            plt,
            [NaN],
            [NaN];
            label=cat,
            marker=:square,
            markersize=8,
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
    min_cells_for_label = 25

    for cat in cats
        idxs = findall(x -> x == cat, labels)

        if length(idxs) >= min_cells_for_label
            xs = [T_on_values[Tuple(I)[2]] for I in idxs]
            ys = [T_off_values[Tuple(I)[1]] for I in idxs]

            annotate!(
                plt,
                mean(xs),
                mean(ys),
                text(cat, 7, :black)
            )
        end
    end

    savefig(plt, filename)
    println("Saved figure: ", filename)

    return plt
end


# =========================================================
# Main
#
# Run the full quick-test phase scan, generate the phase
# diagram, display the plot, and return the label matrix.
# =========================================================

function main()
    p = ParamsB()

    T_on_vals, T_off_vals, labels =
        scan_phase_modelB_dominant_by_test(
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
            lysogen_extinction_threshold=lysogen_extinction_threshold
        )

    plt = plot_phase_modelB_dominant_by_test(
        T_on_vals,
        T_off_vals,
        labels;
        filename=output_filename
    )

    display(plt)

    return T_on_vals, T_off_vals, labels
end


T_on_vals, T_off_vals, labels = main()