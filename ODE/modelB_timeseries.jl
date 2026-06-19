using DifferentialEquations
using Plots
using Printf

# =========================================================
# Model B time-series simulation
#
# This script generates time-series plots for Model B under
# periodic antibiotic pressure.
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
# Post-invasion rule:
#   After invaders are introduced, every test is simulated for
#   exactly 200 complete antibiotic cycles.
#
#   The post-invasion simulation does not stop early, even if
#   the system appears to have converged.
#
# Resident preparation:
#   For resident-based tests, the resident community is first
#   simulated until convergence. Invaders are introduced only
#   after resident convergence.
#
# Plot structure:
#   For each selected (T_on, T_off) point, the script generates
#   one figure with four rows and two columns.
#
#   Each row corresponds to one initial-condition test:
#       Sresident
#       Cresident
#       President
#       DirectPC
#
#   Left column:
#       full 200-cycle post-invasion trajectory
#
#   Right column:
#       final zoomed time window
# =========================================================


# =========================================================
# User settings
# =========================================================

points_to_plot = [
    (600.0, 100.0),
    (600.0, 900.0)
]

test_kinds = (
    :Sresident,
    :Cresident,
    :President,
    :DirectPC
)

rare_density = 1e-4
extinction_threshold = 1e-10

# Resident-only simulations are run until convergence.
resident_ncheck_cycles = 5
resident_convergence_tol = 1e-6
max_resident_cycles = 300

# Post-invasion simulations always run for this fixed number
# of complete antibiotic cycles.
post_invasion_cycles = 200

# Diagnostic convergence check for the final post-invasion window.
# This is used only for reporting in plot titles.
# It does not stop the fixed 200-cycle simulation.
diagnostic_ncheck_cycles = 5
diagnostic_convergence_tol = 1e-6

# ODE solver output and tolerance settings.
saveat_value = 1.0
reltol_value = 1e-7
abstol_value = 1e-9

# Number of final cycles shown in the zoomed right-column panels.
zoom_cycles = 3


# =========================================================
# Parameters
# =========================================================

Base.@kwdef mutable struct ParamsB
    lambda0::Float64 = 1.0
    alpha::Float64   = 1e-6
    a_drug::Float64  = 1.0

    CR::Float64      = 0.05
    CRR::Float64     = 1.01 * 0.05

    phi::Float64     = 1e-8
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

    # Only non-lysogenic hosts are available for new infection.
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
    # Avalue = 0.0 for antibiotic-free phase.
    # Avalue = 1.0 for antibiotic-present phase.
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
# Convergence diagnostics
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


# =========================================================
# Resident-only preparation
#
# Resident communities are run until convergence before rare
# invaders are introduced.
# =========================================================

function run_resident_to_convergence(u0, T_off, T_on, p::ParamsB;
                                     saveat=1.0,
                                     reltol=1e-7,
                                     abstol=1e-9,
                                     max_cycles=300,
                                     ncheck=5,
                                     tol=1e-6)

    current_u = copy(u0)
    current_t = 0.0
    cycle_end_states = Vector{Vector{Float64}}()

    for cyc in 1:max_cycles
        current_u, current_t, _, _ = run_one_complete_cycle(
            current_u,
            current_t,
            T_off,
            T_on,
            p;
            saveat=saveat,
            reltol=reltol,
            abstol=abstol
        )

        push!(cycle_end_states, copy(current_u))

        if has_converged_last_cycles(
            cycle_end_states;
            ncheck=ncheck,
            tol=tol
        )
            return copy(current_u), true, cyc
        end
    end

    return copy(current_u), false, max_cycles
end


# =========================================================
# Fixed post-invasion simulation
#
# After invaders are introduced, the simulation always runs
# for exactly ncycles complete antibiotic cycles.
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
# The order of entries follows the state-variable order used
# by modelB!.
# =========================================================

seed_S_resident() = [
    0.1,    # S_S
    0.0,    # S_R
    0.1,    # L_SS
    0.0,    # L_SR
    0.0,    # L_RS
    0.0,    # L_RR
    1e-8,   # V_S
    0.0,    # V_R
    0.0     # A
]

seed_C_resident() = [
    0.0,    # S_S
    0.1,    # S_R
    0.0,    # L_SS
    0.0,    # L_SR
    0.1,    # L_RS
    0.0,    # L_RR
    1e-8,   # V_S
    0.0,    # V_R
    0.0     # A
]

seed_P_resident() = [
    0.1,    # S_S
    0.0,    # S_R
    0.0,    # L_SS
    0.1,    # L_SR
    0.0,    # L_RS
    0.0,    # L_RR
    0.0,    # V_S
    1e-8,   # V_R
    0.0     # A
]

seed_direct_PC() = [
    0.0,    # S_S
    0.05,   # S_R
    0.0,    # L_SS
    0.1,    # L_SR
    0.05,   # L_RS
    0.0,    # L_RR
    1e-8,   # V_S
    1e-8,   # V_R
    0.0     # A
]


# =========================================================
# Build post-invasion initial conditions
#
# Resident-based tests:
#   1. Run the resident community to convergence.
#   2. Add rare invaders.
#   3. Start the fixed post-invasion simulation.
#
# DirectPC:
#   Starts directly from the mixed initial condition.
# =========================================================

function build_post_invasion_initial_condition(kind::Symbol,
                                               T_off,
                                               T_on,
                                               p::ParamsB;
                                               rare_density=1e-4,
                                               saveat=1.0,
                                               reltol=1e-7,
                                               abstol=1e-9,
                                               max_resident_cycles=300,
                                               ncheck=5,
                                               tol=1e-6)

    if kind == :Sresident
        resident_u, converged, resident_cycles = run_resident_to_convergence(
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
            return nothing, false, resident_cycles
        end

        u0 = copy(resident_u)

        # Rare P invader: prophage-carried ARG strategy.
        u0[4] += rare_density   # L_SR
        u0[8] += rare_density   # V_R

        # Rare C invader: chromosome-carried ARG strategy.
        u0[2] += rare_density   # S_R
        u0[5] += rare_density   # L_RS

        u0[9] = 0.0
        return u0, true, resident_cycles

    elseif kind == :Cresident
        resident_u, converged, resident_cycles = run_resident_to_convergence(
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
            return nothing, false, resident_cycles
        end

        u0 = copy(resident_u)

        # Rare P invader: prophage-carried ARG strategy.
        u0[4] += rare_density   # L_SR
        u0[8] += rare_density   # V_R

        u0[9] = 0.0
        return u0, true, resident_cycles

    elseif kind == :President
        resident_u, converged, resident_cycles = run_resident_to_convergence(
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
            return nothing, false, resident_cycles
        end

        u0 = copy(resident_u)

        # Rare C invader: chromosome-carried ARG strategy.
        u0[2] += rare_density   # S_R
        u0[5] += rare_density   # L_RS

        # Tiny V_S background after introducing chromosome-carried resistance.
        u0[7] += 1e-8           # V_S

        u0[9] = 0.0
        return u0, true, resident_cycles

    elseif kind == :DirectPC
        u0 = seed_direct_PC()
        u0[9] = 0.0
        return u0, true, 0

    else
        error("Unknown test kind: $kind")
    end
end


# =========================================================
# Plot helpers
# =========================================================

function mask_below_threshold(y::AbstractVector; threshold=1e-10)

    # Values below the plotting threshold are hidden at those
    # time points only. A curve can reappear later if the same
    # population rises above the threshold again.
    yy = Float64.(copy(y))
    yy[yy .< threshold] .= NaN

    return yy
end


function prepare_plot_matrix(X::AbstractMatrix; threshold=1e-10)

    # Apply the threshold mask to the eight biological populations.
    # The antibiotic state A is not plotted.
    Y = Matrix{Float64}(undef, size(X, 1), 8)

    for j in 1:8
        Y[:, j] = mask_below_threshold(X[:, j]; threshold=threshold)
    end

    return Y
end


function finite_max(v; default=1e-10)

    # Compute a maximum while ignoring NaN values created by masking.
    vals = [x for x in v if isfinite(x)]

    return isempty(vals) ? default : maximum(vals)
end


function add_antibiotic_shading!(plt, xmin, xmax, T_off, T_on)

    # Gray vertical spans indicate antibiotic-present phases.
    cycle_len = T_off + T_on

    if cycle_len <= 0 || T_on <= 0
        return plt
    end

    k_start = floor(Int, xmin / cycle_len)
    k_end   = ceil(Int, xmax / cycle_len)

    for k in k_start:k_end
        on_left  = k * cycle_len + T_off
        on_right = on_left + T_on

        left  = max(on_left, xmin)
        right = min(on_right, xmax)

        if right > left
            vspan!(
                plt,
                [left, right];
                color=:grey85,
                alpha=0.45,
                label=false
            )
        end
    end

    return plt
end


function build_one_panel(t, Yplot, T_off, T_on, title_str;
                         xwindow=nothing,
                         threshold=1e-10,
                         show_legend=true)

    species_labels = [
        "S_S", "S_R", "L_SS", "L_SR",
        "L_RS", "L_RR", "V_S", "V_R"
    ]

    species_colors = [
        :royalblue,
        :orangered,
        :forestgreen,
        :blueviolet,
        :goldenrod,
        :teal,
        :deeppink,
        :saddlebrown
    ]

    species_styles = [
        :solid, :solid, :solid, :solid,
        :solid, :solid, :dash, :dash
    ]

    # Use the full trajectory unless a zoom window is specified.
    if xwindow === nothing
        idx = 1:length(t)
        xmin, xmax = t[begin], t[end]
    else
        xmin, xmax = xwindow
        idx = findall(tt -> xmin <= tt <= xmax, t)
    end

    tt = t[idx]
    YY = Yplot[idx, :]

    # Set the upper y-limit from finite plotted values.
    ymax = threshold

    for j in 1:size(YY, 2)
        ymax = max(ymax, finite_max(YY[:, j]; default=threshold))
    end

    ymax = max(ymax, threshold * 10)

    plt = plot(
        xlabel="Time",
        ylabel="Density",
        yscale=:log10,
        ylims=(threshold, ymax * 1.2),
        xlims=(xmin, xmax),
        legend=show_legend ? :right : false,
        framestyle=:box,
        grid=true,
        title=title_str,
        titlefontsize=10,
        guidefontsize=10,
        tickfontsize=8,
        legendfontsize=8
    )

    add_antibiotic_shading!(plt, xmin, xmax, T_off, T_on)

    for j in 1:8
        plot!(
            plt,
            tt,
            YY[:, j];
            label=species_labels[j],
            color=species_colors[j],
            linestyle=species_styles[j],
            linewidth=2.2
        )
    end

    return plt
end


function pretty_test_name(kind::Symbol)
    if kind == :Sresident
        return "Sresident"
    elseif kind == :Cresident
        return "Cresident"
    elseif kind == :President
        return "President"
    elseif kind == :DirectPC
        return "DirectPC"
    else
        return string(kind)
    end
end


# =========================================================
# Simulate one post-invasion test
#
# The resident phase may stop at convergence, but the
# post-invasion phase always runs for exactly 200 cycles.
# =========================================================

function simulate_one_test(kind::Symbol, T_on, T_off, p::ParamsB)
    u0, resident_ok, resident_cycles = build_post_invasion_initial_condition(
        kind,
        T_off,
        T_on,
        p;
        rare_density=rare_density,
        saveat=saveat_value,
        reltol=reltol_value,
        abstol=abstol_value,
        max_resident_cycles=max_resident_cycles,
        ncheck=resident_ncheck_cycles,
        tol=resident_convergence_tol
    )

    if !resident_ok
        return nothing, nothing, false, resident_cycles
    end

    t, X, cycle_end_states = run_fixed_post_invasion_cycles(
        u0,
        T_off,
        T_on,
        p;
        ncycles=post_invasion_cycles,
        saveat=saveat_value,
        reltol=reltol_value,
        abstol=abstol_value
    )

    # Diagnostic only: this reports whether the final window is
    # endpoint-stable. It does not affect the simulation length.
    final_window_converged = has_converged_last_cycles(
        cycle_end_states;
        ncheck=diagnostic_ncheck_cycles,
        tol=diagnostic_convergence_tol
    )

    return t, X, final_window_converged, resident_cycles
end


# =========================================================
# Plot one parameter point
#
# For one (T_on, T_off) pair, this function creates eight panels:
#   four full time-series panels
#   four zoomed final-window panels
# =========================================================

function plot_point_timeseries(T_on, T_off, p::ParamsB;
                               threshold=1e-10,
                               zoom_cycles=3)

    panels = Any[]
    cycle_len = T_on + T_off

    for kind in test_kinds
        t, X, final_window_converged, resident_cycles = simulate_one_test(
            kind,
            T_on,
            T_off,
            p
        )

        if t === nothing
            p1 = plot(
                title="$(pretty_test_name(kind)) | resident preparation did not converge",
                framestyle=:box,
                grid=false,
                legend=false
            )

            p2 = plot(
                title="$(pretty_test_name(kind)) | no post-invasion simulation",
                framestyle=:box,
                grid=false,
                legend=false
            )

            push!(panels, p1, p2)
            continue
        end

        Yplot = prepare_plot_matrix(X; threshold=threshold)

        title_full = @sprintf(
            "Model B | %s | full 200 post-invasion cycles | Ton=%.1f, Toff=%.1f | final-window stable=%s",
            pretty_test_name(kind),
            T_on,
            T_off,
            string(final_window_converged)
        )

        title_zoom = @sprintf(
            "Model B | %s | final %d of 200 post-invasion cycles | Ton=%.1f, Toff=%.1f",
            pretty_test_name(kind),
            zoom_cycles,
            T_on,
            T_off
        )

        p_full = build_one_panel(
            t,
            Yplot,
            T_off,
            T_on,
            title_full;
            xwindow=nothing,
            threshold=threshold,
            show_legend=true
        )

        zoom_start = max(t[begin], t[end] - zoom_cycles * cycle_len)

        p_zoom = build_one_panel(
            t,
            Yplot,
            T_off,
            T_on,
            title_zoom;
            xwindow=(zoom_start, t[end]),
            threshold=threshold,
            show_legend=true
        )

        push!(panels, p_full, p_zoom)
    end

    fig = plot(
        panels...;
        layout=(length(test_kinds), 2),
        size=(1800, 2200),
        dpi=250,
        margin=5Plots.mm
    )

    outfile = @sprintf(
        "modelB_timeseries_Ton_%d_Toff_%d.png",
        Int(round(T_on)),
        Int(round(T_off))
    )

    savefig(fig, outfile)
    println("Saved figure: ", outfile)

    return fig
end


# =========================================================
# Main
#
# Generate time-series figures for all parameter points listed
# in points_to_plot.
# =========================================================

function main()
    p = ParamsB()
    figs = Any[]

    for (T_on, T_off) in points_to_plot
        println("Running point: T_on = ", T_on, ", T_off = ", T_off)
        println("Post-invasion cycles fixed at: ", post_invasion_cycles)

        fig = plot_point_timeseries(
            T_on,
            T_off,
            p;
            threshold=extinction_threshold,
            zoom_cycles=zoom_cycles
        )

        push!(figs, fig)
    end

    return figs
end

main()