# Author: Josef Trdla

#= 
There are two functions. Main differences are: 
    * one uses :lazy semantics 
    * one uses :random semantics 
    * lazy one  does not use interpolator (see 'sol(τ)') to determine if transitions were enabled -> less computation needed
=#


function simulate_jt_lazy(H, q0, x0, tf; rng = Random.default_rng())

    sys_map = [sys for sys in H.modes]
    vect_map = vectorfield.(sys_map)
    H = Ref(H)
    
    f(x, p, t) = vect_map[p[1]](x)
    constraints_map = constraints_list.(stateset.(sys_map))
    
    function condition!(out, x, t, integrator) 
        q = integrator.p[1]
        [out[i] = c.b - dot(c.a, x) for (i,c) in enumerate(constraints_map[q])]
    end 

    function affect!(integrator, call_vect) 
        q = integrator.p[1]
        tr = rand(rng, [tr for tr in out_transitions(H[], integrator.p[1]) if isinside(integrator.u, guardset(H[], tr))])
        integrator.p[1] = target(H[], tr)
        integrator.u = apply(resetmap(H[], tr), integrator.u)
    end 

    leaving_invariant = VectorContinuousCallback(condition!, affect!, maximum(length.(constraints_map)))
        
    return solve(ODEProblem(f, x0, (0, tf), q0), Tsit5(); callback=leaving_invariant, abstol=1e-10, reltol=1e-10)
end

function simulate_jt_random(H, q0, x0, tf; rng = Random.default_rng(), rate_var = 5)

    #=
    we use two user params: 
    p[1] ... current mode 
    p[2] ... enabled mode (if no mode is enabled, it equals p[1]!)
    =#
    sys_map = [sys for sys in H.modes]
    vect_map = vectorfield.(sys_map)
    H = Ref(H)
    
    f(x, p, t) = vect_map[p[1]](x)
    constraints_map = constraints_list.(stateset.(sys_map))
    
    function slvr_condition!(out, x, t, integrator) 
        # SOME transition is enabled: out[1]
        trs = [tr for tr in out_transitions(H[], integrator.p[1]) if isinside(x, guardset(H[], tr))]
        if length(trs) >= 1
            out[1] = 1
            integrator.p[2] = target(H[], rand(rng, trs))
        else 
            out[1] = -1
        end 

        # SOME invariant changes: out[1+1]
        [out[i+1] = c.b - dot(c.a, x) for (i,c) in enumerate(constraints_map[integrator.p[1]])]

    end 

    # Not needed here 
    function slvr_affect!(integrator, call_vect) 
        
        # if invariant, switch! 
        if !all(call_vect[2:end].==0)
            integrator.p[1] = integrator.p[2]
        end 

        # else, let jump rate and affect handle the switch! 
    end 

    rate(u, p, t) = rate_var
    
    function jump_affect!(integrator)
        # println("jump aff")
        integrator.p[1] = integrator.p[2]
    end 

    crj = ConstantRateJump(rate, jump_affect!)
    jump_prob = JumpProblem(ODEProblem(f, x0, (0, tf), [q0[1],q0[1]]), Direct(), crj)
    being_enabled = VectorContinuousCallback(slvr_condition!, slvr_affect!, maximum(length.(constraints_map))+1)

    return solve(jump_prob, Tsit5(); callback=being_enabled, abstol=1e-10, reltol=1e-10)
end

