#= 
Author: Zdenek Hurak 
taken from his lectures 
modified a tiny bit by Josef Trdla
=#

#= 
There are two functions. Main differences are: 
    * one uses :lazy semantics 
    * one has semantics as argument
    * :lazy one  does not use interpolator (see 'sol(τ)') to determine if transitions was enabled -> less computation needed
=#

function simulate_hurak_lazy(H, q0, x0, tf; Δt=1e-3, maxjumps=1000, rng = Random.default_rng())

    # initial conds 
    q, x, t = q0[1], float.(collect(x0)), 0.0
    segments = NamedTuple[]
    
    for _ in 1:maxjumps
        sys = mode(H, q)
        # isinside(x, stateset(sys)) || error("The state $x at time $t is not in the invariant of the location $q.")

        # Flow within the location until the state is about to leave the invariant (or until the final time).
        f = vectorfield(sys)
        constraints = constraints_list(stateset(sys))
        condition(out, x, t, integrator) = (out .= [c.b - dot(c.a, x) for c in constraints])  # positive inside
        leaving_invariant = VectorContinuousCallback(condition, (integrator, idx) -> terminate!(integrator), length(constraints))
        sol = solve(ODEProblem((x, p, t) -> f(x), x, (t, tf)), Tsit5(); callback=leaving_invariant, abstol=1e-10, reltol=1e-10)
        t_end = sol.t[end]
        forced = t_end < tf                            # The state is about to leave the invariant before the final time.

        # Time instants (on a grid) at which some outgoing transition is enabled.
        enabled(τ) = [tr for tr in out_transitions(H, q) if isinside(sol(τ), guardset(H, tr))]
        trans = enabled(t_end)

        t_stop = min(tf, t_end)
        ts = [filter(<(t_stop), sol.t); t_stop]
        push!(segments, (q=q, t=ts, x=sol.(ts)))
        isnothing(trans) || isempty(trans) && break

        # Jump: choose one of the enabled transitions and apply its reset map.
        tr = rand(rng, trans)
        q, x, t = target(H, tr), apply(resetmap(H, tr), sol(t_stop)), t_stop
        t >= tf && break
    end

    return segments
end

function simulate_hurak_full(H, q0, x0, tf; semantics=:random, Δt=1e-3, maxjumps=1000, rng = Random.default_rng())
    # initial conds 
    q, x, t = q0[1], float.(collect(x0)), 0.0
    segments = NamedTuple[]
    
    for _ in 1:maxjumps
        sys = mode(H, q)
        # isinside(x, stateset(sys)) || error("The state $x at time $t is not in the invariant of the location $q.")

        # Flow within the location until the state is about to leave the invariant (or until the final time).
        f = vectorfield(sys)
        constraints = constraints_list(stateset(sys))
        condition(out, x, t, integrator) = (out .= [c.b - dot(c.a, x) for c in constraints])  # positive inside
        leaving_invariant = VectorContinuousCallback(condition, (integrator, idx) -> terminate!(integrator), length(constraints))
        sol = solve(ODEProblem((x, p, t) -> f(x), x, (t, tf)), Tsit5(); callback=leaving_invariant, abstol=1e-10, reltol=1e-10)
        t_end = sol.t[end]
        forced = t_end < tf                            # The state is about to leave the invariant before the final time.

        # Time instants (on a grid) at which some outgoing transition is enabled.
        enabled(τ) = [tr for tr in out_transitions(H, q) if isinside(sol(τ), guardset(H, tr))]
        grid = unique([collect(t:Δt:t_end); t_end])
        candidates = filter(τ -> !isempty(enabled(τ)), grid)

        # Choose the switching time (if any).
        τ = if isempty(candidates)
            forced && @warn "Blocking: no transition enabled when leaving the invariant of location $q at t = $t_end."
            nothing
        elseif semantics == :eager
            first(candidates)
        elseif semantics == :lazy
            forced ? last(candidates) : nothing
        elseif semantics == :random
            # If the final time is reached before being forced to leave, staying in the location is also an option.
            forced ? rand(rng, candidates) : rand(rng, [candidates; nothing])
        else
            error("Unknown semantics $semantics.")
        end

        t_stop = something(τ, t_end)
        ts = [filter(<(t_stop), sol.t); t_stop]
        push!(segments, (q=q, t=ts, x=sol.(ts)))
        isnothing(τ) && break

        # Jump: choose one of the enabled transitions and apply its reset map.
        tr = rand(rng, enabled(τ))
        q, x, t = target(H, tr), apply(resetmap(H, tr), sol(τ)), τ
        t >= tf && break
    end
    return segments
end