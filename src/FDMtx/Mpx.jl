# function Mpx(px_pml::SparseMatrixCSC{Float64,Int64}, lambda::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
#     lambda = modExpand(lambda, ext, iflag)
#     (m, n) = size(lambda)
#     a1 = 9/8;   a2 = -1/24;
#     c1 = a1/dx; c2 = a2/dx;
#     C3 = zeros(m*n)
#     denum = zeros(m*n)
#     for ix = 1 : n
#         for iz = 1 : m
#             a = 1/dt + px_pml[iz,ix]/2
#             b = 1/dt - px_pml[iz,ix]/2
#             denum[(ix-1)*m+iz] = lambda[iz,ix] / a
#             C3[(ix-1)*m+iz]= b / a
#         end
#     end
#     tmp = spzeros(n, n)
#     tmp[1,1] =  c1; tmp[1,2] = c2;
#     tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
#     for ix = 3: n-1
#         tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
#         tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
#     end
#     tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
#     MpxBpx = spdiagm(C3)
#     MpxBvx = spdiagm(denum) * kron(tmp, speye(m))
#     return  MpxBpx, MpxBvx
# end

function Mpx(px_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(px_pml)
    a1 = 9/8  ; a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + px_pml[iz,ix]/2
            b = 1/dt - px_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for ix = 3: n-1
        tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
        tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
    end
    tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
    MpxBpx = spdiagm(C3)
    MpxBvx = spdiagm(denum) * kron(tmp, speye(m))
    return MpxBpx, MpxBvx
end

function Mpx_back(px_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(px_pml)
    a1 = 9/8  ; a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt - px_pml[iz,ix]/2
            b = 1/dt + px_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for ix = 3: n-1
        tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
        tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
    end
    tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
    MpxBpx = spdiagm(C3)
    MpxBvx = -spdiagm(denum) * kron(tmp, speye(m))
    return MpxBpx, MpxBvx
end
