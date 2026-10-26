# Author: Josef Trdla

using HybridSystems, MathematicalSystems, LazySets, OrdinaryDiffEq, LinearAlgebra, Random, Plots

using JumpProcesses
using BenchmarkTools

include("hurak_func.jl")
include("jt_func.jl")

## Problem generation 

# note, H is obtained from utils! 
include("utils.jl")
H 

x0 = [22.0]
tf = 30.0
rng = Random.default_rng()
Δt=1e-3
maxjumps=1000
q0 = [2]

## -- SEMANTIC: LAZY

## HURAK FUNC 
@benchmark segments = simulate_hurak_lazy(H, q0, x0, tf; Δt=1e-3, maxjumps=1000, rng = Random.default_rng());
p1 = plot(xlabel="Time (s)", ylabel="Temperature (°C)", title = "simulate_hurak_lazy");
plot_execution!(p1, segments; lw=2, guidefontsize=10, tickfontsize=10, size=(1200, 600), margins = (5, :mm));
hline!(p1,[T_on_max, T_off_min]; color = :red, linestyle = :dash, label = "Invariants");
hline!(p1,[T_guard_on2off, T_guard_off2on]; color = :green, linestyle = :dash, label = "Guards")

## JT FUNC 
@benchmark sol = simulate_jt_lazy(H, q0, x0, tf; rng = Random.default_rng()) 

p1 = plot(xlabel="Time (s)", ylabel="Temperature (°C)", title = "simulate_jt_lazy", size=(1200, 600), margins = (5, :mm));
plot!(p1,sol; lw=2, guidefontsize=10, tickfontsize=10);
hline!(p1,[T_on_max, T_off_min]; color = :red, linestyle = :dash, label = "Invariants");
hline!(p1,[T_guard_on2off, T_guard_off2on]; color = :green, linestyle = :dash, label = "Guards")
## -- SEMANTIC: RANDOM

## HURAK FUNC 
@benchmark segments = simulate_hurak_full(H, q0, x0, tf; semantics = :random, Δt=1e-3, maxjumps=1000, rng = Random.default_rng())
p1 = plot(xlabel="Time (s)", ylabel="Temperature (°C)", title = "simulate_hurak_full");
plot_execution!(p1, segments; lw=2, guidefontsize=10, tickfontsize=10, size=(1200, 600), margins = (5, :mm));
hline!(p1,[T_on_max, T_off_min]; color = :red, linestyle = :dash, label = "Invariants");
hline!(p1,[T_guard_on2off, T_guard_off2on]; color = :green, linestyle = :dash, label = "Guards")

## JT FUNC
@benchmark sol = simulate_jt_random(H, q0, x0, tf; rng = Random.default_rng(), rate_var = 1)

p1 = plot(xlabel="Time (s)", ylabel="Temperature (°C)", title = "simulate_jt_random", size=(1200, 600), margins = (5, :mm));
plot!(p1,sol; lw=2, guidefontsize=10, tickfontsize=10);
hline!(p1,[T_on_max, T_off_min]; color = :red, linestyle = :dash, label = "Invariants");
hline!(p1,[T_guard_on2off, T_guard_off2on]; color = :green, linestyle = :dash, label = "Guards")