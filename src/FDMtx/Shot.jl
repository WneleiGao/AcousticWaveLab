type Shot3
     isz :: Int64
     isx :: Int64
     irz :: Array{Int64, 1}
     irx :: Array{Int64, 1}
     ot  :: Float64
     nt  :: Int64
     dt  :: Float64
     vz  :: Array{Float64, 2}
     vx  :: Array{Float64, 2}
     p   :: Array{Float64, 2}
end

function InitShot3(isz::Int64, isx::Int64, pos::Array{Int64,2}, ot::Float64, nt::Int64, dt::Float64)
    nrec = size(pos, 1)
    irz  = pos[:, 1]
    irx  = pos[:, 2]
    d    = zeros(nt, nrec)
    shot = Shot3(isz, isx, irecz, irecx, ot, nt, dt, d, d, d)
    return shot
end


"""
copy(s3)

copy three components shot gather 

"""
function Base.copy(s::Shot3)
    s1 = Shot3(s.isz, s.isx, s.irz, s.irx, s.ot, s.nt, s.dt, s.vz, s.vx, s.p)
    return s1
end


type Shot
     isz :: Int64
     isx :: Int64
     irz :: Array{Int64, 1}
     irx :: Array{Int64, 1}
     ot  :: Float64
     nt  :: Int64
     dt  :: Float64
     d   :: Array{Float64, 2}
end


function InitShot(isz::Int64, isx::Int64, pos::Array{Int64,2}, ot::Float64, nt::Int64, dt::Float64)
    nrec = size(pos, 1)
    irecz= pos[:, 1]
    irecx= pos[:, 2]
    d    = zeros(nt, nrec)
    shot = Shot(isz, isx, irecz, irecx, ot, nt, dt, d)
    return shot
end


function spt2shot3!(shot::Shot, spt::SnapShot, Nz::Int64, Nx::Int64, ext::Int64, zupper::Int64)
    nrec = length(shot.irx)
    it = spt.it
    for irec = 1 : nrec
        iz = shot.irz[irec] + zupper
        ix = shot.irx[irec] + ext
        ind= (ix-1)*Nz + iz
        shot.d[it, irec] = spt.pz[ind] + spt.px[ind]
    end
    return nothing
end

function AddShot32SnapShot!(spt::SnapShot, shot::Shot, Nz::Int64, Nx::Int64, ext::Int64, zupper::Int64)
    nrec = length(shot.irx)
    it = spt.it
    for irec = 1 : nrec
        iz = shot.irz[irec] + zupper
        ix = shot.irx[irec] + ext
        ind= (ix-1)*Nz + iz
        spt.pz[ind] = spt.pz[ind] + shot.d[it, irec]
        spt.px[ind] = spt.px[ind] + shot.d[it, irec]
    end
end





import Base.copy
function copy(s::Shot)
    s1 = Shot(s.isz, s.isx, s.irz, s.irx, s.ot, s.nt, s.dt, s.d)
    return s1
end

function spt2shot!(shot::Shot, spt::SnapShot, Nz::Int64, Nx::Int64, ext::Int64, zupper::Int64)
    nrec = length(shot.irx)
    it = spt.it
    for irec = 1 : nrec
        iz = shot.irz[irec] + zupper
        ix = shot.irx[irec] + ext
        ind= (ix-1)*Nz + iz
        shot.d[it, irec] = spt.pz[ind] + spt.px[ind]
    end
    return nothing
end

function AddShot2SnapShot!(spt::SnapShot, shot::Shot, Nz::Int64, Nx::Int64, ext::Int64, zupper::Int64)
    nrec = length(shot.irx)
    it = spt.it
    for irec = 1 : nrec
        iz = shot.irz[irec] + zupper
        ix = shot.irx[irec] + ext
        ind= (ix-1)*Nz + iz
        spt.pz[ind] = spt.pz[ind] + shot.d[it, irec]
        spt.px[ind] = spt.px[ind] + shot.d[it, irec]
    end
end

function wfd2shot(pos::Array{Int64,2}, path::ASCIIString)
    (isz, isx, nz, nx, ext, iflag, dt, ot, nt) = InfoWfd(path)
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    shot= InitShot(isz, isx, pos, ot, nt, dt)
    nrec= size(pos, 1)
    fid = open(path, "r"); seek(fid, sizeof(Int64)*8);
    for it = 1 : nt
        d = reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        for irec = 1 : nrec
            iz = shot.irz[irec] + zupper;
            ix = shot.irx[irec] + ext;
            shot.d[it, irec] = d[iz, ix]
        end
    end
    return shot
end

function fullWfd2Shot(pos::Array{Int64,2}, path::ASCIIString)
    (isz, isx, nz, nx, ext, iflag, ot, dt, nt) = InfoFullWfd(path)
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    shot= InitShot(isz, isx, pos, ot, nt, dt)
    nrec= size(pos, 1)
    fid = open(path, "r")
    head = sizeof(Float64) * (8+Nz*Nx*2)
    sizeSpt = sizeof(Float64)*Nz*Nx*4
    for it = 1 : nt
        position = head + sizeSpt*(it-1); seek(fid, position);
        pz = reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        px = reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        p  = pz + px
        for irec = 1 : nrec
            iz = shot.irz[irec] + zupper;
            ix = shot.irx[irec] + ext;
            shot.d[it, irec] = p[iz, ix]
        end
    end
    close(fid)
    return shot
end


function writeShot!(path::ASCIIString, s::Shot)
    fid = open(path, "w")
    write(fid, s.isz); write(fid, s.isx);
    nrec = length(s.irz); write(fid, nrec);
    write(fid, s.irz); write(fid, s.irx);
    write(fid, s.ot); write(fid, s.nt); write(fid, s.dt);
    write(fid, vec(s.d))
    close(fid)
end

function readShot(path::ASCIIString)
    fid = open(path, "r")
    isz = read(fid, Int64); isx = read(fid, Int64);
    nr  = read(fid, Int64)
    irz = read(fid, Int64, nr); irx = read(fid, Int64, nr);
    ot  = read(fid, Float64); nt = read(fid, Int64); dt = read(fid, Float64);
    d   = reshape(read(fid, Float64, nt*nr), nt, nr)
    s   = Shot(isz, isx, irz, irx, ot, nt, dt, d)
    return s
end

import Base.-
function -(shot1::Shot, shot2::Shot)
    if shot1.irz != shot2.irz || shot1.irx != shot2.irx || shot1.nt != shot2.nt
       error("size dismatch")
    end
    pos_rec = hcat(shot1.irz, shot1.irx)
    shot = InitShot(0, 0, pos_rec, 0.0, shot1.nt, shot1.dt)
    shot.d = shot1.d - shot2.d
    return shot
end
