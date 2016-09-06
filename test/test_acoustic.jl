# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos = [101 101]; f0 = 30.0; ot = [0.0];
# src = InitSources(pos, f0, ot, dt);
# # discretize spatial derivative operator
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# path = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path, src, fidMtx, tmax=1.0);
#
# iz = 5 * ones(Int64, 201); ix = collect(1:201);
# pos = hcat(iz, ix)
# shot = fullWfd2Shot(pos, path);
# shot1 = MultiStepForward(pos, src[1], fidMtx, tmax=1.0)
#
# path = "/Users/wenlei/Desktop/adjoint.bin"
# MultiStepAdjoint!(path, shot, fidMtx)
#
#
# using PyPlot
# imshow(read(path, 500))
#
# iz = 5 * ones(Int64, 201); ix = collect(1:201);
# pos = hcat(iz, ix)
# shot = wfd2shot(pos, path);
# figure()
# plot(shot.d[:,71])


# ==============================================================================
# dot product test for one step
# spt1 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 1);
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
# spt2 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 1);
# OneStepForward!(spt2, spt1, fidMtx);
#
# tmp2 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 2)
# tmp2.vz = randn(length(tmp2.vz))
# tmp2.vx = randn(length(tmp2.vx))
# tmp2.pz = randn(length(tmp2.pz))
# tmp2.px = randn(length(tmp2.px))
# tmp1 = InitSnapShot(101, 101, nz, nx, ext, 1, dt, 2)
# OneStepAdjoint!(tmp1, tmp2, fidMtx)
#
# x  = vcat(spt1.vz, spt1.vx, spt1.pz, spt1.px)
# x1 = vcat(tmp1.vz, tmp1.vx, tmp1.pz, tmp1.px)
# y  = vcat(spt2.vz, spt2.vx, spt2.pz, spt2.px)
# y1 = vcat(tmp2.vz, tmp2.vx, tmp2.pz, tmp2.px)
#
# x'*x1 - y'*y1



# ==============================================================================
# nz = 201; nx = 201; dt = 2e-3; iflag = 1; ext = 20
# if iflag ==1
#    zupper = ext
#    Nz = nz + 2*ext
# elseif iflag == 2
#    Nz = nx +   ext
#    zupper = ext
# end
# Nx = nx + 2*ext
# spt1 = InitSnapShot(101, 101, nz, nx, ext, iflag, dt, 1);
# spt1.vz = randn(length(spt1.vz));
# spt1.vx = randn(length(spt1.vx));
# spt1.pz = randn(length(spt1.pz));
# spt1.px = randn(length(spt1.px));
#
# iz = collect(1:5:201); ix = collect(1:3:201);
# nrz= length(iz); nrx = length(ix);
# iz = vec(repmat(iz, nrx, 1));
# ix = vec(repmat(ix', (nrz), 1));
# pos = hcat(iz, ix);
# shot1 = InitShot(101, 101, pos, 0.0, 1, dt)
# spt2shot!(shot1, spt1, Nz, Nx, ext, zupper)
#
#
# shot2 = InitShot(101, 101, pos, 0.0, 1, dt)
# shot2.d = randn(1, size(shot2.d, 2))
# spt2 = InitSnapShot(101, 101, nz, nx, ext, iflag, dt, 1);
# AddShot2SnapShot!(spt2, shot2, Nz, Nx, ext, zupper)
#
# x  = vcat(spt1.vz, spt1.vx, spt1.pz, spt1.px)
# x1 = vcat(spt2.vz, spt2.vx, spt2.pz, spt2.px)
# y  = vec(shot1.d)
# y1 = vec(shot2.d)
# x'*x1 - y'*y1


# ==============================================================================
# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos = [101 101]; f0 = 30.0; ot = [0.0];
# src = InitSources(pos, f0, ot, dt);
# # discretize spatial derivative operator
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# path = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path, src, fidMtx, tmax=1.0);
#
# iz = 5 * ones(Int64, 201); ix = collect(1:201);
# pos = hcat(iz, ix)
# shot = fullWfd2Shot(pos, path);
# shot1 = MultiStepForward(pos, src[1], fidMtx, tmax=1.0)
#
# path = "/Users/wenlei/Desktop/adjoint.bin"
# MultiStepAdjoint!(path, shot, fidMtx)

# ==============================================================================
# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3; f0 = 30.0
# v  = 3000. * ones(nz, nx);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# # initialize source term
# ix  = collect(1:2:201); iz = 101*ones(Int64, length(ix));
# pos = hcat(iz, ix)
# f0 = 30.0; ot = zeros(size(pos,1));
# src = InitSources(pos, f0, ot, dt);
# path = "/Users/wenlei/Desktop/src2spts.bin"
# src2SnapShots(path, src, nz, nx, ext, iflag);
#
# # position of receiver
# irz = []; irx = []
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            irz = push!(irz, iz)
#            irx = push!(irx, ix)
#         end
#     end
# end
# pos_rec = convert(Array{Int64,2}, hcat(irz, irx))
# shot  = MultiStepForward(pos_rec, path, fidMtx, tmax=1.0)
# path_wfd = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path_wfd, src , fidMtx, tmax=1.0)
# shot1 = fullWfd2Shot(pos_rec, path_wfd)

# ==============================================================================
# using AcousticWave
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3; f0 = 30.0
# v  = 3000. * ones(nz, nx);
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# # initialize source term
# isz = []; isx = []
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() <= 0.1
#            isz = push!(isz, iz)
#            isx = push!(isx, ix)
#         end
#     end
# end
# pos_src = convert(Array{Int64,2}, hcat(isz, isx))
# f0 = 30.0; ot = zeros(size(pos_src,1));
# src = InitSources(pos_src, f0, ot, dt);
# tmax = 1; nt = round(Int64, tmax/dt+1);
# w = Ricker(f0, dt); ns = nt + 1 - length(w)
# for isrc = 1 : length(src)
#     src[isrc].nt= nt
#     src[isrc].p = conv(randn(ns),w)
# end
# path = "/Users/wenlei/Desktop/src2spts.bin"
# src2SnapShots(path, src, nz, nx, ext, iflag);
# # position of receiver
# irx = collect(1:2:nx); irz = 5*ones(Int64,length(irx));
# pos_rec = hcat(irz, irx)
# shot  = MultiStepForward(pos_rec, path, fidMtx, tmax=1.0)
#
#
# path_adj = "/Users/wenlei/Desktop/adjoint.bin"
# shot1 = InitShot(0, 0, pos_rec, 0.0, nt, dt)
# for irec = 1 : length(shot1.irz)
#     shot1.d[:,irec] = conv(randn(ns), w)
# end
# MultiStepAdjoint_test(path_adj, shot1, fidMtx)
#
# L = 0.0
# for ir = 1 : length(shot.irx)
#     L = L + dot(shot.d[:,ir], shot1.d[:,ir])
# end
#
# T = 0.0
# for it = 1 : nt
#     spt1 = readSnapShot(path    , it)
#     spt2 = readSnapShot(path_adj, it)
#     T = T + dot(spt1.vz, spt2.vz) + dot(spt1.vx, spt2.vx) + dot(spt1.pz, spt2.pz) + dot(spt1.px, spt2.px)
# end

# ==============================================================================
# using AcousticWave
# path = "/Users/wenlei/Desktop/acoustic.bin"
# nz = 201; nx = 201; ext = 20; iflag = 1; dz = 10.; dx = 10.; dt = 2e-3;
# v  = 3000. * ones(nz, nx);
# # initialize source term
# pos = [101 101]; f0 = 30.0; ot = [0.0];
# src = InitSources(pos, f0, ot, dt);
# # discretize spatial derivative operator
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
# path = "/Users/wenlei/Desktop/acoustic.bin"
# MultiStepForward(path, src, fidMtx, tmax=1.0);
#
# using PyPlot, AcousticWave
# path = "/Users/wenlei/Desktop/acoustic.bin"
# pathout = "/Users/wenlei/Desktop/homo"
# waveAnim(path, pathout)
# sptPlot(path, 100)
# irx = collect(1:2:nx); irz = 5*ones(Int64,length(irx));
# pos_rec = hcat(irz, irx)
# shot = fullWfd2Shot(pos_rec, path)
