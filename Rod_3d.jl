using WaterLily
using StaticArrays
function spin(L=2^5; Re=1000, U=1, ϵ=0.5, thk=2ϵ+√2, period=200)
    # Line segment SDF
    function sdf(x, t)
        y = x .- SA[0, clamp(x[2], -L/2, L/2)]
        √sum(abs2, y) - thk/2
    end
    # Oscillating motion and rotation
    function map(x, t)
        β = 2π * t / period  # 360 degrees rotation over the given period
        R = SA[
            cos(β) -sin(β) 0;
            sin(β) cos(β) 0;
            0 0 1
            ]
        R * (x - SA[3L, 3L])
    end
    Simulation((6L, 6L, L), (0, 0, 0), L; U, ν=U*L/Re, body=AutoBody(sdf, map), ϵ)
end

include("2D.jl")

sim_gif!(spin(), duration=8, clims=(-10, 10), plotbody=true, remeasure=true)