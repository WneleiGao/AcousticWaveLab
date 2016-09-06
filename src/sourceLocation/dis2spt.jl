function dis2spt(w::Array{Float64,1}, dis::Array{Float64,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if length(w) < it
       error("out of source time range")
    end
    if length(dis) != nz*nx
       error("length of dis does not much model size")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    spt = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
    tmp = zeros(Nz, Nx); dis_tmp = reshape(dis, nz, nx)
    tmp[zupper+1: zupper+nz, ext+1:ext+nx] = 1/2 * w[it] * dis_tmp
    tmp = vec(tmp)
    spt.pz[:] = tmp[:]
    spt.px[:] = tmp[:]
    return spt
end


function spt2dis(w::Array{Float64,1}, spt::SnapShot)
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    it = spt.it
    if length(w) < spt.it
       error("out of source time range")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    pz = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    px = reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    dis = 1/2 * w[it] * pz + 1/2 * w[it] * px
    return vec(dis)
end


function spt2wlet(w::Array{Float64,1}, dis::Array{Float64,1}, spt::SnapShot)
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    it = spt.it
    if length(w) < it
       error("out of source time range")
    end
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    p = reshape(spt.pz, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx] + reshape(spt.px, Nz, Nx)[zupper+1:zupper+nz, ext+1:ext+nx]
    p = vec(p)
    w[it] = 1/2 * dot(dis, p)
    return nothing
end
