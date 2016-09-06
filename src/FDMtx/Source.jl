type Source
     iz :: Int64
     ix :: Int64
     ot :: Float64
     dt :: Float64
     nt :: Int64
     p  :: Array{Float64,1}
end

function Ricker(f0::Float64, dt::Float64)
	  nw = 2.2/f0/dt
	  nw = 2*floor(Int,nw/2)+1
	  nc = floor(Int,nw/2)
	  w  = zeros(nw)
	  k  = collect(1:nw)
    k  = vec(k)
	  alpha = (nc-k+1)*f0*dt*pi
	  beta = alpha.^2
	  w = (1.-beta.*2).*exp(-beta)
	  return w
end

function InitSources(pos::Array{Int64, 2}, f0::Float64, ot::Array{Float64, 1}, dt::Float64)
    nsrc = size(pos, 1)
    if nsrc != length(ot)
       error("length of position and ot dismatch")
    end
    src = Array(Source, nsrc)
    for isrc = 1 : nsrc
        p = Ricker(f0, dt)
        nt= length(p)
        src[isrc] = Source(pos[isrc,1], pos[isrc,2], ot[isrc], dt, nt, p)
    end
    return src
end

function AddSources!(spt::SnapShot, src::Array{Source,1})
    dt = spt.dt
    if spt.dt != src[1].dt
       error("time sampling dismatch")
    end
    nz = spt.nz
    nx = spt.nx
    ext= spt.ext
    iflag = spt.iflag
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
       Nx  = nx + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
       Nx  = nx + 2*ext
    end
    it  = spt.it
    nsrc = length(src)
    for isrc = 1: nsrc
        if src[isrc].ot <= (it-1)*dt <= src[isrc].ot+(src[isrc].nt-1)*dt
           indt = it - round(Int64, src[isrc].ot/dt)
           indz = src[isrc].iz + zupper
           indx = src[isrc].ix + ext
           pos = (indx-1) * Nz + indz
           spt.pz[pos] = spt.pz[pos] + src[isrc].p[indt]/2
           spt.px[pos] = spt.px[pos] + src[isrc].p[indt]/2
        end
    end
    return nothing
end

function AddSources!(spt::SnapShot, w::Array{Float64,1}, dis::Array{Float64,1})
    it = spt.it
    if it > length(w)
       error("it is out the range of wavelet")
    end
    nz = spt.nz ; nx = spt.nx;
    ext= spt.ext; iflag = spt.iflag;
    if length(dis) != nz*nx
       error("size of dis dismatch with the size of SnapShot")
    end
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
       Nx  = nx + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
       Nx  = nx + 2*ext
    end
    tmp = zeros(Nz, Nx); dis_tmp = reshape(dis, nz, nx)
    tmp[zupper+1: zupper+nz, ext+1:ext+nx] = 1/2 * w[it] * dis_tmp
    tmp = vec(tmp)
    spt.pz[:] = spt.pz[:] + tmp[:]
    spt.px[:] = spt.px[:] + tmp[:]
    return nothing
end

function AddSource!(spt::SnapShot, src::Source)
    dt = spt.dt
    if spt.dt != src.dt
       error("time sampling dismatch")
    end
    nz = spt.nz
    nx = spt.nx
    ext= spt.ext
    iflag = spt.iflag
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
       Nx  = nx + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
       Nx  = nx + 2*ext
    end
    it  = spt.it
    if src.ot <= (it-1)*dt <= src.ot+(src.nt-1)*dt
       indt = it - round(Int64, src.ot/dt)
       indz = src.iz + zupper
       indx = src.ix + ext
       pos = (indx-1) * Nz + indz
      #  spt.vz[pos] = spt.vz[pos] + src.p[indt]
       spt.pz[pos] = spt.pz[pos] + (src.p[indt])/2
       spt.px[pos] = spt.px[pos] + (src.p[indt])/2
    end
    return nothing
end


function SrcRange(src::Array{Source, 1})
    tmax = 0.0
    tmin = 0.0
    nsrc = length(src)
    dt   = src[1].dt
    for isrc = 1: nsrc
        tmp = src[isrc].ot
        tmp1= tmp + (src[isrc].nt-1)*dt
        if tmp < tmin
           tmin = tmp
        end
        if tmp1 > tmax
           tmax = tmp1
        end
    end
    return tmin, tmax
end

function src2spt(src::Array{Source,1}, it::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64)
    dt = src[1].dt
    spt = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, it)
    AddSources!(spt, src)
    return spt
end

function src2SnapShots(path::ASCIIString, src::Array{Source,1}, nz::Int64, nx::Int64, ext::Int64, iflag::Int64)
    dt = src[1].dt
    (tl, tu) = SrcRange(src)
    nt  = round(Int64, (tu-tl)/dt) + 1
    fid = open(path, "w")
    write(fid, src[1].iz, src[1].ix, nz, nx, ext, iflag, tl, dt)
    for it = 1 : nt
        it1 = round(Int64, (tl+(it-1)*dt)/dt+1)
        spt = src2spt(src, it1, nz, nx, ext, iflag)
        writeSnapShot(fid, spt)
    end
    close(fid)
    return nothing
end

type SrcSpt
     nz :: Int64
     nx :: Int64
     ot :: Float64
     dt :: Float64
     it :: Int64
     p  :: Array{Float64, 1}
end

function InitSrcSpt(nz::Int64, nx::Int64, ot::Float64, dt::Float64, it::Int64)
    p = zeros(nz*nx)
    srcspt = SrcSpt(nz, nx, ot, dt, it, p)
    return srcspt
end

function src2srcSpts(path::ASCIIString, src::Array{Source,1}, nz::Int64, nx::Int64)
    ns = length(src)
    dt = src[1].dt
    (tl, tu) = SrcRange(src)
    nt  = round(Int64, (tu-tl)/dt) + 1
    fid = open(path, "w")
    write(fid, nz, nx, tl, src[1].dt, nt)
    for it = 1 : nt
        p   = zeros(nz*nx)
        for is = 1 : ns
            if src[is].ot <= tl+(it-1)*dt <= src[is].ot + dt*(src[is].nt-1)
               indt = it - round(Int64, (src[is].ot-tl)/dt)
               ind = (src[is].ix-1)*nz + src[is].iz
               p[ind] = p[ind] + src[is].p[indt]
            end
        end
        write(fid, p)
    end
    close(fid)
    return nothing
end

function readSrcSpt(path::ASCIIString, it::Int64)
    fid = open(path, "r")
    nz = read(fid, Int64)  ; nx = read(fid, Int64)  ;
    ot = read(fid, Float64); dt = read(fid, Float64);
    position = sizeof(Int64)*5 + sizeof(Float64)*nz*nx*(it-1)
    seek(fid, position)
    p = read(fid, Float64, nz*nx)
    p = reshape(p, nz, nx)
    close(fid)
    return p
end

function srcSpt2wlet(path::ASCIIString, iz::Int64, ix::Int64)
    fid = open(path, "r")
    (nz, nx, ot, dt, nt) = InfoSrcSpt(path)
    w = zeros(nt)
    for it = 1 : nt
        p = readSrcSpt(path, it)
        w[it] = p[iz, ix]
    end
    return w
end

function InfoSrcSpt(path::ASCIIString)
    fid = open(path, "r")
    nz = read(fid, Int64)
    nx = read(fid, Int64)
    ot = read(fid, Float64)
    dt = read(fid, Float64)
    nt = read(fid, Int64)
    close(fid)
    return nz, nx, ot, dt, nt
end


function AddSrcSpt2Spt(spt::SnapShot, path::ASCIIString)
    it = spt.it
    nz = spt.nz; nx = spt.nx;
    ext = spt.ext; iflag = spt.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    p = readSrcSpt(path, it) * 1/2
    tmp = reshape(spt.pz, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
    spt.pz = vec(tmp)
    tmp = reshape(spt.px, Nz, Nx)
    tmp[zupper+1:nz+zupper, ext+1:ext+nx] = tmp[zupper+1:nz+zupper, ext+1:ext+nx] + p
    spt.px = vec(tmp)
    return nothing
end


function saveSpt2SrcSpt(path::ASCIIString, spt::SnapShot, nz::Int64, nx::Int64, nt::Int64, ext::Int64, Nz::Int64, Nx::Int64, zupper::Int64)
    fid = open(path, "w")
    write(fid, nz, nx, 0.0, spt.dt, nt)
    tmp = 1/2*(reshape(spt.pz, Nz, Nx) + reshape(spt.px, Nz, Nx))
    tmp = vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])
    write(fid, tmp)
    flush(fid)
    return fid
end


function saveSpt2SrcSpt(fid::IOStream, spt::SnapShot, nz::Int64, nx::Int64, ext::Int64, Nz::Int64, Nx::Int64, zupper::Int64)
    tmp = 1/2*(reshape(spt.pz, Nz, Nx) + reshape(spt.px, Nz, Nx))
    tmp = vec(tmp[zupper+1:zupper+nz, ext+1:ext+nx])
    write(fid, tmp)
    flush(fid)
    return nothing
end

function reverseSrcSptOrder(path::ASCIIString, path_tmp::ASCIIString)
    (nz, nx, ot, dt, nt) = InfoSrcSpt(path_tmp)
    fid = open(path, "w")
    write(fid, nz, nx, ot, dt, nt)
    fid_tmp = open(path_tmp, "r")
    for it = nt : -1 : 1
        position = sizeof(Float64)*5 + sizeof(Float64)*nz*nx*(it-1)
        seek(fid_tmp, position)
        d = read(fid_tmp, Float64, nz*nx);
        write(fid, d)
    end
    close(fid_tmp); close(fid)
    rm(path_tmp)
    return nothing
end
