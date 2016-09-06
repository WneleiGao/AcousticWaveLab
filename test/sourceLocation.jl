using PyPlot, AcousticWave
# =============generate synthetic data set================================
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.0];
src = InitSources(pos_src, f0, ot, dt);
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix);
shot = MultiStepForward(pos_rec, src, fidMtx, tmax=0.5);
wt = Ricker(f0, dt); ntw = length(wt);
mu = 0.1; lambda = 120.0; tmax = 0.5; maxit = 10;


(dis, obj) =CG_dis(shot, fidMtx, pos_rec, wt, maxit=10, mu=0.001, tmax=0.5)
dis = reshape(dis, nz, nx)
SeisPlot(dis, name="/Users/wenlei/Desktop/CG.pdf", dx=10, dy=10, xlabel="Length", xunits="(m)", ylabel="Depth", yunits="(m)")
dis = zeros(Float64,nz,nx)
dis[pos_src[1,1], pos_src[1,2]] = 1.; dis[pos_src[2,1], pos_src[2,2]] = 1.;
dis = vec(dis)
ntw = 5*ntw;
(w, obj) = CG_w(ntw, shot, fidMtx, pos_rec, dis, maxit=10, mu=0.001, tmax=0.5)
plot(w);
plot(wt);
xlable = "Time (s)"
ylable = "Amplitude"
savefig("/Users/wenlei/Desktop/wavelet.pdf")

dis = fista(shot, wt, mu, lambda, fidMtx, pos_rec, tmax, maxit)

# (dis, w) = alterMini(w, shot, fidMtx, pos_rec, mu, lambda, outer_it=6, fistait=10, CGit=5, tmax=0.5)
# (dis, w) = alternateCG(w, shot, fidMtx, pos_rec, outer_it=6)




# ============================================================================================================
# using PyPlot, AcousticWave
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
# v  = 2500. * ones(nz, nx);
# v[101:end, 1:end] = 3000.;
# vmax = maximum(v); vmin = minimum(v);
#
# Nz = nz + 2*ext; Nx = nx + 2*ext
# w = Ricker(f0, dt); dis = randn(nz*nx); it = floor(Int64, length(w)/2);
# spt = dis2spt(w, dis, nz, nx, ext, iflag, dt, it);
# # SeisPlot(reshape(spt.pz, Nz, Nx));
# # SeisPlot(reshape(spt.px, Nz, Nx));
#
# tmp = reshape(spt.pz, Nz, Nx) + reshape(spt.px, Nz, Nx);
# tmp = tmp[ext+1:ext+nz, ext+1:ext+nx];
# norm(w[it]*dis-vec(tmp))
#
# spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
# dis1 = spt2dis(w, spt1);
# tmp_spt = dot(spt.vz, spt1.vz) + dot(spt.vx, spt1.vx) + dot(spt.pz, spt1.pz) + dot(spt.px, spt1.px)
# tmp_dis = dot(dis, dis1)

# ===========================================================================================================
# using PyPlot, AcousticWave
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
# v  = 2500. * ones(nz, nx);
# v[101:end, 1:end] = 3000.;
# vmax = maximum(v); vmin = minimum(v);
# w = Ricker(f0, dt); it = 33;
# spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
# dis  = zeros(nz, nx); dis[101,75] = 1.0; dis[101,125] = 1.0; dis = vec(dis);
# AddSources!(spt2, w, dis)
# spt3 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
# f0 = 30.0; ot = [0.0,0.0]; pos_src = [101 75;101 125]
# src = InitSources(pos_src, f0, ot, dt);
# AddSources!(spt3, src)
# norm(spt2.vx-spt3.vx)
# norm(spt2.vz-spt3.vz)
# norm(spt2.pz-spt3.pz)
# norm(spt2.px-spt3.px)

# =============================================================================================================
using PyPlot, AcousticWave
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.0];
src = InitSources(pos_src, f0, ot, dt);
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix);
w = Ricker(f0, dt);
dis = zeros(nz, nx); dis[pos_src[1,1], pos_src[1,2]] = 1.0;
dis[pos_src[2,1], pos_src[2,2]] = 1.0; dis = vec(dis);
shot = MultiStepForward(pos_rec, w, dis, fidMtx, tmax=0.5);
dis1 = MultiStepAdjoint(w, shot, fidMtx)

# ============================================================================================================
# dot product test
# using PyPlot, AcousticWave
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.0];
# src = InitSources(pos_src, f0, ot, dt);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
# pos_rec = hcat(iz, ix);
# w = Ricker(f0, dt);
# dis = zeros(nz, nz);
# for ix = 1: nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            dis[iz, ix] = randn()
#         end
#     end
# end
# dis[pos_src[2,1], pos_src[2,2]] = 1.0; dis = vec(dis);
# shot = MultiStepForward(pos_rec, w, dis, fidMtx, tmax=0.5);
#
# nt = shot.nt
# shot1 = InitShot(0, 0, pos_rec, 0.0, nt, dt)
# nrec  = size(pos_rec, 1)
# ls    = floor(Int64, length(w)/2)
# for irec = 1 : nrec
#     tmp  = conv(rand(nt), w)
#     shot1.d[:, irec] = tmp[ls+1:end-ls]
# end
# dis1 = MultiStepAdjoint(w, shot1, fidMtx)
#
# tmp  = dot(dis, dis1)
# tmp1 = dot(vec(shot.d), vec(shot1.d))

# ============================================================================================================
# using PyPlot, AcousticWave
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3; f0 = 30.0
# v  = 2500. * ones(nz, nx);
# v[101:end, 1:end] = 3000.;
# vmax = maximum(v); vmin = minimum(v);
# Nz = nz + 2*ext; Nx = nx + 2*ext
# w = Ricker(f0, dt); dis = randn(nz*nx);
# nt = length(w)
# spts = Array(SnapShot, nt)
# for it = 1 : nt
#     spts[it] = dis2spt(w, dis, nz, nx, ext, iflag, dt, it)
# end
# it = 36
# a = reshape(spts[it].pz, 261, 261);
# b = a[ext+1:ext+nz, ext+1:ext+nx];
# tmp = b ./ reshape(dis, nz, nx)
#
#
# spts1 = Array(SnapShot, nt)
# for it = 1 : nt
#     spts1[it] = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
#     spts1[it].vz = randn(length(spts1[it].vz))
#     spts1[it].vx = randn(length(spts1[it].vx))
#     spts1[it].pz = randn(length(spts1[it].pz))
#     spts1[it].px = randn(length(spts1[it].px))
# end
# w1 = zeros(nt)
# for it = 1 : nt
#     spt2wlet(w1, dis, spts1[it])
# end
#
# tmp = 0.0
# for it = 1 : nt
#     tmp = tmp + dot(spts[it].vz, spts1[it].vz)
#     tmp = tmp + dot(spts[it].vx, spts1[it].vx)
#     tmp = tmp + dot(spts[it].pz, spts1[it].pz)
#     tmp = tmp + dot(spts[it].px, spts1[it].px)
# end
# tmp1 = dot(w, w1)


# ============================================================================================================
# dot product test
# using PyPlot, AcousticWave
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.0];
# src = InitSources(pos_src, f0, ot, dt);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
# pos_rec = hcat(iz, ix);
# w = Ricker(f0, dt);
# dis = zeros(nz, nz);
# for ix = 1: nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            dis[iz, ix] = randn()
#         end
#     end
# end
# dis = vec(dis)
# shot = MultiStepForward(pos_rec, w, dis, fidMtx, tmax=0.5);
#
#
# nt = shot.nt
# shot1 = InitShot(0, 0, pos_rec, 0.0, nt, dt)
# nrec  = size(pos_rec, 1)
# ls    = floor(Int64, length(w)/2)
# for irec = 1 : nrec
#     tmp  = conv(rand(nt), w)
#     shot1.d[:, irec] = tmp[ls+1:end-ls]
# end
# ntw = length(w)
# w1  = MultiStepAdjoint(ntw, dis, shot1, fidMtx)
#
# tmp  = dot(w, w1)
# tmp1 = dot(vec(shot.d), vec(shot1.d))

# ===========================================================================================================
using PyPlot, AcousticWave
# =============generate synthetic data set================================
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.0];
src = InitSources(pos_src, f0, ot, dt);
vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix); tmax = 0.5; maxit=20;
w = Ricker(f0, dt)
lambda = power(maxit, fidMtx, w, pos_rec, tmax)
