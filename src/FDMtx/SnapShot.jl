type SnapShot
     isz:: Int64
     isx:: Int64
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dt :: Float64
     it :: Int64
     vz :: Array{Float64, 1}
     vx :: Array{Float64, 1}
     pz :: Array{Float64, 1}
     px :: Array{Float64, 1}
end

function InitSnapShot(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if iflag == 1
       Nz = nz + 2*ext
       Nx = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       Nx = nx + 2*ext
    end
    spt = SnapShot(isz, isx, nz, nx, ext, iflag, dt, it, zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx))
    return spt
end

import Base.copy
function copy(snapShot1::SnapShot, snapShot2::SnapShot)
    if snapShot2.nz != snapShot1.nz || snapShot2.nx != snapShot1.nx || snapShot2.ext != snapShot1.ext || snapShot2.dt != snapShot1.dt || snapShot2.iflag != snapShot1.iflag
       error("the two snapShot are different")
    end
    snapShot1.isz = snapShot2.isz
    snapShot1.isx = snapShot2.isx
    snapShot1.it  = snapShot2.it
    snapShot1.vz[:] = snapShot2.vz[:]
    snapShot1.vx[:] = snapShot2.vx[:]
    snapShot1.pz[:] = snapShot2.pz[:]
    snapShot1.px[:] = snapShot2.px[:]
    return nothing
end

import Base.norm
function norm(snapShot::SnapShot)
    L = norm(snapShot.vz)
    L = L + norm(snapShot.vx)
    L = L + norm(snapShot.vz)
    L = L + norm(snapShot.pz)
    L = L + norm(snapShot.px)
    return L
end

import Base.-
function -(snapShot1::SnapShot, snapShot2::SnapShot)
    isz= snapShot1.isz
    isx= snapShot1.isx
    nz = snapShot1.nz
    nx = snapShot1.nx
    ext = snapShot1.ext
    iflag = snapShot1.iflag
    dt = snapShot1.dt
    it = snapShot1.it
    if nz != snapShot2.nz || nx != snapShot2.nx || ext != snapShot2.ext || iflag != snapShot2.iflag
       error("size dismatch")
    end
    snapShot3 = InitSnapShot(isz, isx, nz, nx, ext, iflag, dt, it)
    snapShot3.vz = snapShot1.vz - snapShot2.vz
    snapShot3.vx = snapShot1.vx - snapShot2.vx
    snapShot3.pz = snapShot1.pz - snapShot2.pz
    snapShot3.px = snapShot1.px - snapShot2.px
    return snapShot3
end

import Base.+
function +(snapShot1::SnapShot, snapShot2::SnapShot)
    isz= snapShot1.isz
    isx= snapShot1.isx
    nz = snapShot1.nz
    nx = snapShot1.nx
    ext = snapShot1.ext
    iflag = snapShot1.iflag
    dt = snapShot1.dt
    it = snapShot1.it
    if nz != snapShot2.nz || nx != snapShot2.nx || ext != snapShot2.ext || iflag != snapShot2.iflag
       error("size dismatch")
    end
    snapShot3 = InitSnapShot(isz, isx, nz, nx, ext, iflag, dt, it)
    snapShot3.vz = snapShot1.vz + snapShot2.vz
    snapShot3.vx = snapShot1.vx + snapShot2.vx
    snapShot3.pz = snapShot1.pz + snapShot2.pz
    snapShot3.px = snapShot1.px + snapShot2.px
    return snapShot3
end

import Base.write
function write(path::ASCIIString, spt::SnapShot, ot::Float64, dt::Float64)
    fid = open(path, "w")
    write(fid, spt.isz); write(fid, spt.isx  );
    write(fid, spt.nz ); write(fid, spt.nx   );
    write(fid, spt.ext); write(fid, spt.iflag);
    write(fid, ot     ); write(fid, dt       );
    write(fid, spt.pz+spt.px)
    flush(fid)
    return fid
end

function write(fid::IOStream, spt::SnapShot)
    write(fid, vec(spt.pz+spt.px))
    return nothing
end

import Base.read
function read(path::ASCIIString, it::Int64)
    fid = open(path, "r") ; seek(fid, sizeof(Int64)*2);
    nz  = read(fid, Int64); nx = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    if iflag == 1
       Nz = nz + 2*ext
       Nx = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       Nx = nx + 2*ext
    end
    position = sizeof(Int64)*8 + (it-1)*sizeof(Float64)*Nz*Nx
    seek(fid, position)
    d = read(fid, Float64, Nz*Nx)
    d = reshape(d, Nz, Nx)
    close(fid)
    return d
end

function readSnapShot(path::ASCIIString, it::Int64)
    fid = open(path, "r") ;
    isz = read(fid, Int64); isx   = read(fid, Int64)
    nz  = read(fid, Int64); nx    = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    ot= read(fid, Float64); dt  = read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    position = sizeof(Int64)*8 + (it-1)*sizeof(Float64)*Nz*Nx*4
    seek(fid, position)
    vz = read(fid, Float64, Nz*Nx)
    vx = read(fid, Float64, Nz*Nx)
    pz = read(fid, Float64, Nz*Nx)
    px = read(fid, Float64, Nz*Nx)
    spt = SnapShot(isz, isx, nz, nx, ext, iflag, dt, it, vz, vx, pz, px)
    close(fid)
    return spt
end

function InfoWfd(path::ASCIIString; print_flag = false)
    fid = open(path, "r")
    isz = read(fid, Int64); isx = read(fid, Int64);
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    ot  = read(fid, Float64); dt =read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*8) / (Nz*Nx*sizeof(Float64)))
    close(fid)
    if print_flag
       println("shot z: $isz, shot x: $isx")
       println("zlength: $nz, xlength: $nx")
       println("padding: $ext, surface: $iflag")
       println("ot: $ot, dt: $dt, nt: $nt")
    end
    return isz, isx, nz, nx, ext, iflag, ot, dt, nt
end

function InfoFullWfd(path::ASCIIString; print_flag = false)
    fid = open(path, "r")
    isz = read(fid, Int64); isx = read(fid, Int64);
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    ot  = read(fid, Float64); dt =read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*8) / (Nz*Nx*4*sizeof(Float64)))
    close(fid)
    if print_flag
       println("shot z: $isz, shot x: $isx")
       println("zlength: $nz, xlength: $nx")
       println("padding: $ext, surface: $iflag")
       println("ot: $ot, dt: $dt, nt: $nt")
    end
    return isz, isx, nz, nx, ext, iflag, ot, dt, nt
end

function numberOfSnapShots(path::ASCIIString)
    fid = open(path, "r")
    isz = read(fid, Int64); isx = read(fid, Int64);
    nz  = read(fid, Int64); nx  = read(fid, Int64);
    ext = read(fid, Int64); iflag=read(fid, Int64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    nt = round(Int64, (filesize(fid) - sizeof(Int64)*8) / (Nz*Nx*4*sizeof(Float64)))
    close(fid)
    return nt
end

function writeSnapShot(path::ASCIIString, spt::SnapShot)
    fid = open(path, "w")
    write(fid, spt.isz); write(fid, spt.isx  );
    write(fid, spt.nz ); write(fid, spt.nx   );
    write(fid, spt.ext); write(fid, spt.iflag);
    write(fid, 0.0    ); write(fid, spt.dt   );
    write(fid, spt.vz ); write(fid, spt.vx   );
    write(fid, spt.pz ); write(fid, spt.px   );
    flush(fid)
    return fid
end

function writeSnapShot(fid::IOStream, spt::SnapShot)
    write(fid, spt.vz ); write(fid, spt.vx   );
    write(fid, spt.pz ); write(fid, spt.px   );
    flush(fid)
    return nothing
end

function reverseSnapShotsOrder(path::ASCIIString, path_tmp::ASCIIString)
    (isx, isy, nz, nx, ext, iflag, ot, dt, nt) = InfoFullWfd(path_tmp)
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nx = nx +   ext
    end
    Nx = nx + 2*ext
    fid = open(path, "w")
    write(fid, isx, isy, nz, nx, ext, iflag, ot, dt)
    fid_tmp = open(path_tmp, "r")
    for it = nt:-1:1
        position = sizeof(Float64)*8 + sizeof(Float64)*Nz*Nx*4*(it-1)
        seek(fid_tmp, position)
        d = read(fid_tmp, Float64, Nz*Nx*4);
        write(fid, d)
    end
    close(fid_tmp); close(fid)
    rm(path_tmp)
    return nothing
end
