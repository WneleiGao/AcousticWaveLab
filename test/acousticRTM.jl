using PyPlot, AcousticWave
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
v  = 2500. * ones(nz, nx);
v[101:end, 1:end] = 3000.;
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# initialize source term
f0 = 30.0; ot = [0.0]; pos_src = [5 101]
src = InitSources(pos_src, f0, ot, dt);
# receiver location
irx = collect(1:2:nx); irz = 5*ones(Int64, length(irx));
pos_rec = convert(Array{Int64,2}, hcat(irz, irx))
# position of receiver
shot  = MultiStepForward(pos_rec, src[1], fidMtx, tmax=1.0)
SeisPlot(shot.d, pclip=95, cmap="seismic", name="/Users/wenlei/Desktop/shot.pdf")
dp = dpDirect(shot, v[1], src[1], dx, dt)
shot.d = dp .* shot.d
SeisPlot(shot.d, cmap="seismic", name="/Users/wenlei/Desktop/shot_rd.pdf")

# compute source side wavefield with smooth velocity model
v1 = modSmooth(v, 10)
fidMtx1 = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v1);
pathsrc = "/Users/wenlei/Desktop/src.bin"
MultiStepForward(pathsrc, src, fidMtx1, tmax=1.0)
pathout  = "/Users/wenlei/Desktop/srcWave"
waveAnim(pathsrc, pathout, interval=100)

# compute receiver side wavefield with smooth velocity model
pathrec = "/Users/wenlei/Desktop/rec.bin"
MultiStepAdjoint(pathrec, shot, fidMtx1)
animRec = "/Users/wenlei/Desktop/recWave"
waveAnim(pathrec, animRec, interval=100)
(Ivz, Ivx, Ipz, Ipx) = RTMimaging(pathsrc, pathrec, image_type="all");
SeisPlot(Ipz, cmap="seismic")
Ip = RTMimaging(pathsrc, pathrec, image_type="p");
SeisPlot(Ip, cmap="seismic") 

# ==============================================================================
using PyPlot, AcousticWave
nz = 101; nx = 101; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
v  = 2500. * ones(nz, nx);
v[101:end, 1:end] = 3000.;
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
m31 = fidMtx.MpzBvz * fidMtx.MvzBvz;
m33 = fidMtx.MpzBvz * fidMtx.MvzBp + fidMtx.MpzBpz;
m34 = fidMtx.MpzBvz * fidMtx.MvzBp;

m42 = fidMtx.MpxBvx * fidMtx.MvxBvx;
m43 = fidMtx.MpxBvx * fidMtx.MvxBp;
m44 = fidMtx.MpxBvx * fidMtx.MvxBp + fidMtx.MpxBpx;

N = (nz+2*ext)*(nx+2*ext);
M = vcat(hcat(fidMtx.MvzBvz, spzeros(N,N) , fidMtx.MvzBp, fidMtx.MvzBp),
         hcat(spzeros(N,N) , fidMtx.MvxBvx, fidMtx.MvxBp, fidMtx.MvxBp),
         hcat(m31          , spzeros(N,N) , m33         , m34         ),
         hcat(spzeros(N,N) , m42          , m43         , m44         ))

# ==============================================================================
# using PyPlot, AcousticWave
# nz = 101; nx = 101; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
# v  = 3000. * ones(nz, nx);
# v[101:end, 1:end] = 3000.;
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# fidMtx1= InitFidMtx_back(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
#
# function OneStepBackward(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
#     spt2.it = spt1.it - 1
#     spt2.pz = fidMtx1.MpzBpz * spt1.pz + fidMtx1.MpzBvz * spt1.vz;
#     spt2.px = fidMtx1.MpxBpx * spt1.px + fidMtx1.MpxBvx * spt1.vx;
#     spt2.vz = fidMtx1.MvzBvz * spt1.vz + fidMtx1.MvzBp * (spt2.pz+spt2.px);
#     spt2.vx = fidMtx1.MvxBvx * spt1.vx + fidMtx1.MvxBp * (spt2.pz+spt2.px);
#     return nothing
# end
#
# spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 2);
# spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 2);
# spt3 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 2);
# spt1.vz = rand(length(spt1.vz)); spt1.vx = rand(length(spt1.vx));
# spt1.pz = rand(length(spt1.pz)); spt1.px = rand(length(spt1.px));
# OneStepBackward(spt2, spt1, fidMtx);
# OneStepAdjoint!(spt3, spt1, fidMtx);
#
# p2 = spt2.pz + spt2.px
# p3 = spt3.pz + spt3.px
#
# tmp = spt3-spt2
