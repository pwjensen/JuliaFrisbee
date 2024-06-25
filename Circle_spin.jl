using WaterLily
using StaticArrays
function hover(L=2^5;Re=1000,U=1,amp=π/4,ϵ=0.5,thk=2ϵ+√2, period=200)
    # Line segment SDF
    function sdf(x,t)
        y = x .- SA[0,clamp(x[2],-L/2,L/2)]
        √sum(abs2,y)-thk/2
    end
    # Oscillating motion and rotation
    function map(x,t)
        β = 2π * t / period  # 360 degrees rotation over the given period
        R = SA[cos(β) -sin(β); sin(β) cos(β)]
        R * (x - SA[3L, 3L])
    end
    Simulation((6L,6L),(0,0),L;U,ν=U*L/Re,body=AutoBody(sdf,map),ϵ)
end

include("2D.jl")

sim_gif!(hover(),duration=8,clims=(-10,10),plotbody=true, remeasure=true)