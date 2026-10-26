# Author: Zdenek Hurak 
# taken from his lectures 
## -- Plotting 

function plot_execution!(p, segments; i=1, colors=(:red, :blue), kwargs...)
    for (k, s) in enumerate(segments)
        plot!(p, s.t, getindex.(s.x, i); color=colors[s.q], label="", kwargs...)
        if k < length(segments)                         # Dashed line for the jump (if any) to the next segment.
            next = segments[k+1]
            plot!(p, [s.t[end], next.t[1]], [s.x[end][i], next.x[1][i]]; color=:gray, ls=:dash, label="", kwargs...)
        end
    end
    return p
end

## -- Generation 
T_on_max=20.0
T_off_min=16.0
T_guard_on2off=19.0
T_guard_off2on=17.0

automaton = GraphAutomaton(2)                       # Locations: 1 = on, 2 = off
    
add_transition!(automaton, 1, 2, 1)                 # Transition 1: on -> off
add_transition!(automaton, 2, 1, 2)                 # Transition 2: off -> on

# Flows and invariants. A halfspace {x | a'x ≤ b} is given as HalfSpace(a, b).
mode_on = @system(x' = -0.1*x + 5.0, x ∈ HalfSpace([1.0], T_on_max), dim=1)     # x ≤ T_on_max
mode_off = @system(x' = -0.1*x, x ∈ HalfSpace([-1.0], -T_off_min), dim=1)        # x ≥ T_off_min
modes = [mode_on, mode_off]

# Guards (and identity reset maps)
on2off = ConstrainedIdentityMap(1, HalfSpace([-1.0], -T_guard_on2off))           # x ≥ T_guard_on2off
off2on = ConstrainedIdentityMap(1, HalfSpace([1.0], T_guard_off2on))             # x ≤ T_guard_off2on

H = HybridSystem(automaton, modes, [on2off, off2on], [AutonomousSwitching()])

vectorfield(s::ConstrainedLinearContinuousSystem) = x -> s.A * x
vectorfield(s::ConstrainedAffineContinuousSystem) = x -> s.A * x + s.c

# Membership in a polyhedral set (given by halfspaces) with some tolerance (to tolerate numerical errors at the boundary).
isinside(x, X; tol=1e-8) = all(dot(c.a, x) <= c.b + tol for c in constraints_list(X))

guardset(H, tr) = stateset(resetmap(H, tr))