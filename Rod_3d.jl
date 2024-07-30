using WaterLily
using StaticArrays

function spin(p=5; Re=1000, mem=Array, U=1)
    # Period, length, viscosity
    n=2^p; period = 200; L = 2^5;  ν=U*L/Re
    # Motion Function
    α = 2π / period

    # Lots of fixes needed
    A(t) = SA[cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1] * (ω - SA[3L, 3L, L])

    # Issues with what variables are doing and how to convert only to rotating the jellyfish
    sphere = AutoBody((x,t)->abs(√sum(abs2,x)-L)-1,
            (x,t)->A(t).*x)
    plane =  AutoBody((x,t)->x[3]-h,
            (x,t)->A(t).*x)
    body = sphere-plane

    Simulation((8L, 6L, 16), (0, 0, 0), L; U, ν, body, mem)
end


using Meshing, GeometryBasics
function geom!(md,d,sim,t=WaterLily.time(sim))
    a = sim.flow.σ
    WaterLily.measure_sdf!(a,sim.body,t)
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
    normal_mesh(GeometryBasics.Mesh(md,Meshing.MarchingCubes(),origin=Vec(0,0,0),widths=size(md)))
end

function ω!(md,d,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
end

function mirrorto!(a,b)
    n = size(b,1)
    a[reverse(1:n),reverse(1:n),:].=b
    a[reverse(n+1:2n),1:n,:].=a[1:n,1:n,:]
    a[:,reverse(n+1:2n),:].=a[:,1:n,:]
    return a
end

using GLMakie
Makie.inline!(false)
#CUDA.allowscalar(false)
begin
    # Define geometry and motion on GPU
    sim = spin();
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