# function Mvz(vz_pml::SparseMatrixCSC{Float64,Int64}, rho::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     rho = modExpand(rho, ext, iflag)
#     (m,n) = size(rho)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             if iz < m
#                a = (rho[iz+1,ix]+rho[iz,ix]) / (2*dt)
#             else
#                a = rho[iz,ix] / dt
#             end
#             b = vz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = 1 / (a+b)
#             C3[(ix-1)*m+iz]= (a-b) / (a+b)
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
#     for iz = 2: m-2
#         tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
#         tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
#     end
#     tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
#     tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
#     MvzBvz = spdiagm(C3)
#     MvzBp  = spdiagm(denum) * kron(speye(n), tmp)
#     return MvzBvz, MvzBp
# end


function Mvz(vz_pml::SparseMatrixCSC{Float64,Int64}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vz_pml[iz,ix]/2
            b = 1/dt - vz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for iz = 2: m-2
        tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
        tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
    end
    tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
    tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
    MvzBvz = spdiagm(C3)
    MvzBp  = spdiagm(denum) * kron(speye(n), tmp)
    return MvzBvz, MvzBp
end

function Mvz_back(vz_pml::SparseMatrixCSC{Float64,Int64}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    (m,n) = size(vz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt - vz_pml[iz,ix]/2
            b = 1/dt + vz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for iz = 2: m-2
        tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
        tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
    end
    tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
    tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
    MvzBvz = spdiagm(C3)
    MvzBp  = -spdiagm(denum) * kron(speye(n), tmp)
    return MvzBvz, MvzBp
end
