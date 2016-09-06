# using PyPlot, AcousticWave
# # =============dot product test for joint inversion================================
# nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
# v  = 3000. * ones(nz, nx);
#
# # initialize source term
# isx = []; isz = []; pic = zeros(nz, nx)
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() < 0.001
#            push!(isx, ix); push!(isz, iz)
#            pic[iz, ix] = 1.0
#         end
#     end
# end
# pos_src = convert(Array{Int64,2}, hcat(isz, isx)); f0 = 30.0; ot = zeros(size(pos_src,1));
# src = InitSources(pos_src, f0, ot, dt);
# (tl, tu) = SrcRange(src)
# ntw = round(Int64, (tu-tl)/dt+1)
# path = "/Users/wenlei/Desktop/srcSpt.bin"
# src2srcSpts(path, src, nz, nx)
#
# # finite difference coefficients
# vmax = maximum(v); vmin = minimum(v);
# fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
#
# # receiver positions
# irx = []; irz = []; pic_rec = zeros(nz, nx)
# for ix = 1 : nx
#     for iz = 1 : nz
#         if rand() < 0.01
#            push!(irx, ix); push!(irz, iz)
#            pic_rec[iz, ix] = 1.0
#         end
#     end
# end
# pos_rec = convert(Array{Int64,2}, hcat(irz, irx));
#
# # modeling
# shot = MultiStepForward(pos_rec, path, fidMtx, tmax=0.5);
# shot1= MultiStepForward(pos_rec, src, fidMtx, tmax=0.5);
#
# nt = shot.nt
# shot2 = InitShot(0, 0, pos_rec, 0.0, nt, dt)
# nrec  = size(pos_rec, 1)
# w = Ricker(40.0, dt)
# ls = floor(Int64, length(w)/2)
# for irec = 1 : nrec
#     tmp  = conv(rand(nt), w)
#     shot2.d[:, irec] = tmp[ls+1:end-ls]
# end
# path_adj = "/Users/wenlei/Desktop/adjoint.bin"
# MultiStepAdjoint(ntw, path_adj, shot2, fidMtx)
#
# tmp = dot(vec(shot.d), vec(shot2.d))
# tmp1 = 0.0
# for it = 1 : ntw
#     p = readSrcSpt(path, it)
#     p1 = readSrcSpt(path_adj, it)
#     tmp1 = tmp1 + dot(vec(p), vec(p1))
# end


using PyPlot, AcousticWave
# =============the result of adjoint================================
nz = 201; nx = 201; ext = 30; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = zeros(size(pos_src,1));
src = InitSources(pos_src, f0, ot, dt);
(tl, tu) = SrcRange(src)
ntw = round(Int64, (tu-tl)/dt+1)

vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix);
shot = MultiStepForward(pos_rec, src, fidMtx, tmax=0.5);
# path_adj = "/Users/wenlei/Desktop/adj.bin"
# MultiStepAdjoint(ntw, path_adj, shot, fidMtx)

path_x = "/Users/wenlei/Desktop/x.bin"
path_s = "/Users/wenlei/Desktop/s.bin"
path_p = "/Users/wenlei/Desktop/p.bin"
CGjoint(ntw, path_x, path_s, path_p, shot, fidMtx, pos_rec, tmax=0.5, mu=0.1, niter=30)
SeisPlot(readSrcSpt(path_x, 10), dx=10, dy=10, xlabel="Length", xunits="(m)", ylabel="Depth", yunits="(m)")

shot_pre = MultiStepForward(pos_rec, path_x, fidMtx, tmax=0.5)
it = 40
plot(shot.d[:,it])
plot(shot_pre.d[:,it])



nz = 100; nx = 77; ot = 0.0; dt = 0.001; nt = 33;
path1 = "/Users/wenlei/Desktop/tmp1.bin"
path2 = "/Users/wenlei/Desktop/tmp2.bin"
fid1 = open(path1, "w");
fid2 = open(path2, "w");
write(fid1, nz, nx, ot, dt, nt)
write(fid2, nz, nx, ot, dt, nt)
for it = 1 : nt
    write(fid1,     ones(nz*nx))
    write(fid2, 2.0*ones(nz*nx))
end
close(fid1)
close(fid2)

W = ones(nz*nx)
mu = -0.5
obtainS(path1, mu, W, path2)
