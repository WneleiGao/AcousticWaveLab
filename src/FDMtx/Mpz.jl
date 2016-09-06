# function Mpz(pz_pml::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     lambda = modExpand(lambda, ext, iflag)
#     (m, n) = size(lambda)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dz; c2 = a2/dz;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + pz_pml[iz,ix]/2
#             b = 1/dt - pz_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(m, m)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#     for iz = 3: m-1
#         tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
#         tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
#     end
#     tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
#     MpzBpz = spdiagm(C3)
#     MpzBvz = spdiagm(denum) * kron(speye(n), tmp)
#     return MpzBpz, MpzBvz
# end

function Mpz(pz_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(pz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + pz_pml[iz,ix]/2
            b = 1/dt - pz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for iz = 3: m-1
        tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
        tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
    end
    tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
    MpzBpz = spdiagm(C3)
    MpzBvz = spdiagm(denum) * kron(speye(n), tmp)
    return MpzBpz, MpzBvz
end

function Mpz_back(pz_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(pz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt - pz_pml[iz,ix]/2
            b = 1/dt + pz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for iz = 3: m-1
        tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
        tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
    end
    tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
    MpzBpz = spdiagm(C3)
    MpzBvz = -spdiagm(denum) * kron(speye(n), tmp)
    return MpzBpz, MpzBvz
end
