using WaterLily
using StaticArrays
function spin(L=2^5; Re=1000, U=1, period=200, mem=Array)
    # Line segment SDF
    function sdf(x)
        center,r = SA[3L, 3L, 0], L/2
        x, y, z = xyz - center
        √sum(abs2, SA[x, y, 0]) - r
    end
    # Oscillating motion and rotation
    function map(x, t)
        β = 2π * t / period  # 360 degrees rotation over the given period
        R = SA[
            cos(β) -sin(β) 0;
            sin(β) cos(β) 0;
            0 0 1
            ]
        R * (x - SA[3L, 3L, L])
    end
    Simulation((8L, 6L, 16), (0, 0, 0), L; U, ν=U*L/Re, body=AutoBody(sdf, map), mem)
end

sim= spin();

include("3D.jl")

# sim_gif!(spin(), duration=8, clims=(-10, 10), plotbody=true, remeasure=true)
makie_video!(sim, remeeasure = true, name = "3D_Spin.mp4", duration=6)