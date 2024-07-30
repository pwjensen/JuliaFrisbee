using WaterLily
using StaticArrays
using WGLMakie
norm(x::StaticArray) = sqrt(x'*x)

function make_sim(;n=2^6, U = 1, Re =1000, mem=Array) # mem=CUDA.CuArray
	R = n/2
	function sdf(xyz,t)
			x,y,z = xyz
			r = norm(SA[y,z]); r-R # cylinder
			norm(SA[r-min(r,R), y])-1.5 #disk
	end
	#map(xyz, t) = xyz # -SA[2n/3,0,0]
	function map(x, t) # missing y and z components
        β = 2π * t / period
        R = SA[
            cos(β) -sin(β);
            sin(β) cos(β)
        ]
        R * (x - SA[3L, 3L])
    end
	return Simulation((2n,2n,2n), (U,0,0), R;
		ν=U*R/Re, body=AutoBody(sdf,map), mem)
end

sim = make_sim();
sim_step!(sim, 0.1)

include("3D.jl")
#using CUDA
#@assert CUDA.functional()

# WGLMakie.mesh(body_mesh(sim))
#GLMakie.mesh(normal_mesh(sim))
# WGLMakie.mesh(normal_mesh(sim))

#=
# Ripped from examples
using GLMakie
Makie.inline!(false)
#CUDA.allowscalar(false)
begin
    # Define geometry and motion on GPU
    sim_step!(sim,sim_time(sim)+0.05);

    # Create CPU buffer arrays for geometry flow viz
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(2,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.1, color=:red)

    #Set up flow viz
    ω = ω!(md,d,sim) |> Observable;
    volume!(ω, algorithm=:mip, colormap=:algae, colorrange=(1,10))
    fig
end

# Loop in time
record(fig,"rod.mp4",1:100) do frame
#foreach(1:10) do frame
    @show frame
    sim_step!(sim,sim_time(sim)+0.05);
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
end
=#