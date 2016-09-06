# function Mvx(vx_pml::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     rho = modExpand(rho, ext, iflag)
#     (m,n) = size(rho)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dx; c2 = a2/dx;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             if ix < n
#                a = (rho[iz,ix+1]+rho[iz,ix]) / (2*dt)
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vx_pml[iz,ix] / 2
#             denum[(ix-1)*m+iz] = 1 / (a+b)
#             C3[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#     tmp = spzeros(n, n)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for ix = 2: n-2
#         tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
#         tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
#     end
#     tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
#     tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
#     MvxBvx = spdiagm(C3)
#     MvxBp = spdiagm(denum) * kron(tmp, speye(m))
#     return MvxBvx, MvxBp
# end

function Mvx(vx_pml::SparseMatrixCSC{Float64,Int64}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vx_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vx_pml[iz,ix]/2
            b = 1/dt - vx_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for ix = 2: n-2
        tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
        tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
    end
    tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
    tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
    MvxBvx = spdiagm(C3)
    MvxBp  = spdiagm(denum) * kron(tmp, speye(m))
    return MvxBvx, MvxBp
end


function Mvx_back(vx_pml::SparseMatrixCSC{Float64,Int64}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vx_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt - vx_pml[iz,ix]/2
            b = 1/dt + vx_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for ix = 2: n-2
        tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
        tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
    end
    tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
    tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
    MvxBvx = spdiagm(C3)
    MvxBp  = -spdiagm(denum) * kron(tmp, speye(m))
    return MvxBvx, MvxBp
end
