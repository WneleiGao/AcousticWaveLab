type FidMtx
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dz :: Float64
     dx :: Float64
     dt :: Float64
     MvzBvz :: SparseMatrixCSC{Float64,Int64}
     MvzBp  :: SparseMatrixCSC{Float64,Int64}
     MvxBvx :: SparseMatrixCSC{Float64,Int64}
     MvxBp  :: SparseMatrixCSC{Float64,Int64}
     MpzBpz :: SparseMatrixCSC{Float64,Int64}
     MpzBvz :: SparseMatrixCSC{Float64,Int64}
     MpxBpx :: SparseMatrixCSC{Float64,Int64}
     MpxBvx :: SparseMatrixCSC{Float64,Int64}
end


function DispStable!(vmax::Float64, vmin::Float64, f0::Float64, dz::Float64, dx::Float64, dt::Float64)
    h = vmin / 5 / f0
    tt = 6*dx / (7*sqrt(3)*vmax)
    if dx >= h || dz >= h || dt >= tt
       println("maximum acceptable grid size: $h")
       println("maximum acceptable time step: $tt")
       error("unstable or frequency dispersion")
    end
    return nothing
end

function InitFidMtx(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmax::Float64, vmin::Float64, f0::Float64, v::Array{Float64,2})
    DispStable!(vmax, vmin, f0, dz, dx, dt)
    Cpml = InitCpml(nz, nx, ext, vmax, dz, dx, iflag)
    # Cpml = InitCpml(nz, nx, ext, vmax, rho, dz, dx, iflag)
    (MvzBvz, MvzBp)  = Mvz(Cpml.vz,    dz, dt, ext, iflag)
    (MvxBvx, MvxBp)  = Mvx(Cpml.vx,    dx, dt, ext, iflag)
    (MpzBpz, MpzBvz) = Mpz(Cpml.pz, v, dz, dt, ext, iflag)
    (MpxBpx, MpxBvx) = Mpx(Cpml.px, v, dx, dt, ext, iflag)
    fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
    return fidMtx
end

function InitFidMtx_back(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmax::Float64, vmin::Float64, f0::Float64, v::Array{Float64,2})
    DispStable!(vmax, vmin, f0, dz, dx, dt)
    Cpml = InitCpml(nz, nx, ext, vmax, dz, dx, iflag)
    # Cpml = InitCpml(nz, nx, ext, vmax, rho, dz, dx, iflag)
    (MvzBvz, MvzBp)  = Mvz_back(Cpml.vz,    dz, dt, ext, iflag)
    (MvxBvx, MvxBp)  = Mvx_back(Cpml.vx,    dx, dt, ext, iflag)
    (MpzBpz, MpzBvz) = Mpz_back(Cpml.pz, v, dz, dt, ext, iflag)
    (MpxBpx, MpxBvx) = Mpx_back(Cpml.px, v, dx, dt, ext, iflag)
    fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
    return fidMtx
end



# function InitFidMtx(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmax::Float64, vmin::Float64, f0::Float64, rho::Array{Float64,2}, lambda::Array{Float64,2})
#     DispStable!(vmax, vmin, f0, dz, dx, dt)
#     Cpml = InitCpml(nz, nx, ext, vmax, rho, iflag)
#     # Cpml = InitCpml(nz, nx, ext, vmax, rho, dz, dx, iflag)
#     (MvzBvz, MvzBp)  = Mvz(Cpml.vz, rho,    dz, dt, ext, iflag)
#     (MvxBvx, MvxBp)  = Mvx(Cpml.vx, rho,    dx, dt, ext, iflag)
#     (MpzBpz, MpzBvz) = Mpz(Cpml.pz, lambda, dz, dt, ext, iflag)
#     (MpxBpx, MpxBvx) = Mpx(Cpml.px, lambda, dx, dt, ext, iflag)
#     fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
#     return fidMtx
# end

function vp2lambda(vp::Array{Float64,2}, rho::Array{Float64,2})
    (m , n ) = size(vp)
    (m1, n1) = size(rho)
    if m1 != m || n1 != n
       error("size not match")
    end
    lambda = Array(typeof(vp[1]), m, n)
    for j = 1 : n
        for i = 1 : m
            lambda[i,j] = rho[i,j] * vp[i,j]^2
        end
    end
    return lambda
end

function lambda2vp(lambda::Array{Float64, 2}, rho::Array{Float64, 2})
    (m, n) = size(lambda)
    (m1, n1) = size(rho)
    if m1 != m || n1 != n
       error("size not match")
    end
    vp = Array(typeof(lambda[1]), m, n)
    for j = 1 : n
        for i = 1 : m
            vp[i,j] = sqrt(lambda[i,j] / rho[i,j])
        end
    end
    return vp
end
