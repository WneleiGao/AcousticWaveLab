using PyPlot, AcousticWave
# =============generate synthetic data set================================
nz = 201; nx = 201; ext = 20; iflag = 1; dz = 5.; dx = 5.; dt = 1e-3;
v  = 3000. * ones(nz, nx);
# initialize source term
pos_src = [101 75; 101 125]; f0 = 30.0; ot = [0.0, 0.01];
src = InitSources(pos_src, f0, ot, dt);
path = "/Users/wenlei/Desktop/srcSpt.bin"
src2srcSpts(path, src, nz, nx)
p = readSrcSpt(path, 47); p = reshape(p, nz, nx);



vmax = maximum(v); vmin = minimum(v);
fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v);
iz = collect(5:5:201); ix = 20 * ones(Int64, length(iz));
pos_rec = hcat(iz, ix);
shot = MultiStepForward(pos_rec, src, fidMtx, tmax=0.5);
