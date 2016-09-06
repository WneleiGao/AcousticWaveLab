doc"""
    OneStepForward!(spt2, spt1, fidMtx)

One time step forward modeling of acoustic wave field, the current snapshot is
kept in spt2. $ p^{n+1} = F \times p^{n} + s^{n+1} $

# Arguments
* `spt2 :: SnapShot`: composite type of SnapShot
* `spt1 :: SnapShot`: composite type of SnapShot
* `fidMtx :: FidMtx`: composite type of sparse partial differential matrix
"""
function OneStepForward!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
    spt2.it  = spt1.it + 1
    spt2.vz = fidMtx.MvzBvz * spt1.vz + fidMtx.MvzBp  * (spt1.pz+spt1.px)
    spt2.vx = fidMtx.MvxBvx * spt1.vx + fidMtx.MvxBp  * (spt1.pz+spt1.px)
    spt2.pz = fidMtx.MpzBpz * spt1.pz + fidMtx.MpzBvz *  spt2.vz
    spt2.px = fidMtx.MpxBpx * spt1.px + fidMtx.MpxBvx *  spt2.vx
    return nothing
end

"""
    MultiStepForward(pos_rec, srcs, fidMtx)

Multiple time step forward modeling of acoustic wave field, generate one shot gather, support
inject multiple sources simutaneously

# Arguments
* `pos_rec :: Array{Int64,2}`: receivers' location, first column specify vertical index of receivers
second column specify horizontal index.
* `srcs :: Array{Source,1}`: composite type for point source injection.
* `fidMtx :: FidMtx`: composite type of sparse partial differential matrix

# Output
* `shot :: Shot`: composite type of common shot gather
"""
function MultiStepForward(pos::Array{Int64,2}, src::Array{Source,1}, fidMtx::FidMtx; tmax=2.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    if iflag == 1
       zupper =    ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper =    0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    dt = fidMtx.dt
    (tl, tu) = SrcRange(src)
    nt = round(Int64, tmax/dt)+1
    shot = InitShot(0, 0, pos, 0.0, nt, dt)
    spt1 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 2)
    AddSources!(spt1, src)
    spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu   #Add sources
           AddSources!(spt2, src)
        end
        copy(spt1, spt2)
        spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return shot
end

"""
    MultiStepForward(pos_rec, src, fidMtx)

Multiple time step forward modeling of acoustic wave field, generate one shot gather

# Arguments
* `pos_rec :: Array{Int64,2}`: receivers' location, first column specify vertical index of receivers
second column specify horizontal index.
* `src :: Source`: composite type for point source injection.
* `fidMtx :: FidMtx`: composite type of sparse partial differential matrix

# Output
* `shot :: Shot`: composite type of common shot gather
"""
function MultiStepForward(pos::Array{Int64,2}, src::Source, fidMtx::FidMtx; tmax=2.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    dt = fidMtx.dt
    nt = round(Int64, tmax/dt)+1
    isz= src.iz; isx= src.ix;
    shot = InitShot(isz, isx, pos, 0.0, nt, dt)
    tl = src.ot; tu = tl + (src.nt-1)*dt;
    spt1 = InitSnapShot(src.iz, src.ix, nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(src.iz, src.ix, nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu   #Add sources
           AddSource!(spt2, src)
        end
        copy(spt1, spt2)
        spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return shot
end

"""
    MultiStepForward(path, srcs, fidMtx)

Multiple time step forward modeling of acoustic wave field, generate one shot gather, support
inject multiple sources simutaneously

# Arguments
* `path :: ASCIIString`: the directory for writing binary wavefiled
* `srcs :: Array{Source,1}`: array of composite source type
* `fidMtx :: FidMtx`: composite type of sparse partial differential matrix

# Output
* `shot :: Shot`: composite type of common shot gather
"""
function MultiStepForward(path::ASCIIString, src::Array{Source,1}, fidMtx::FidMtx; tmax=2.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    dt = fidMtx.dt ;
    (tl, tu) = SrcRange(src)
    nt = round(Int64, tmax/dt)+1
    spt1 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 2)
    AddSources!(spt1, src)
    fid = writeSnapShot(path, spt1)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu   #Add sources
           AddSources!(spt2, src)
        end
        copy(spt1, spt2)
        writeSnapShot(fid, spt1)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    close(fid)
    return nothing
end

# inject sources as SrcSpt, return shot
function MultiStepForward(pos::Array{Int64,2}, path::ASCIIString, fidMtx::FidMtx; tmax=1.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    dt = fidMtx.dt
    nt = round(Int64, tmax/dt+1)
    (nz1, nx1, ot, dt1, nt1) = InfoSrcSpt(path)
    if nz != nz1 || nx != nx1 || dt != dt1
       error("size dismatch")
    end
    shot = InitShot(0, 0, pos, ot, nt, dt)
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 1)
    AddSrcSpt2Spt(spt1, path)
    spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if it <= nt1
           AddSrcSpt2Spt(spt2, path)
        end
        copy(spt1, spt2)
        spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return shot
end

#  forward modeling for source location
function MultiStepForward(pos::Array{Int64,2}, w::Array{Float64,1}, dis::Array{Float64,1}, fidMtx::FidMtx; tmax=1.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    dt = fidMtx.dt
    nt = round(Int64, tmax/dt+1)
    ntsrc = length(w)
    shot = InitShot(0, 0, pos, 0.0, nt, dt)
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 1)
    AddSources!(spt1, w, dis)
    spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 1)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if it <= ntsrc
           AddSources!(spt2, w, dis)
        end
        copy(spt1, spt2)
        spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return shot
end

# function MultiStepForward!(path::ASCIIString, src::Array{Source,1}, fidMtx::FidMtx; tmax=2.0, samp=1)
#     nz = fidMtx.nz ; nx = fidMtx.nx;
#     ext= fidMtx.ext; iflag = fidMtx.iflag;
#     dt = fidMtx.dt ;
#     (tl, tu) = SrcRange(src)
#     nt = round(Int64, tmax/dt)+1
#     spt1 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 1)
#     spt2 = InitSnapShot(src[1].iz, src[1].ix, nz, nx, ext, iflag, dt, 2)
#     AddSources!(spt1, src)
#     fid = write(path, spt1, 0.0, dt)
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, fidMtx)
#         if tl <= (it-1)*dt <= tu   #Add sources
#            AddSources!(spt2, src)
#         end
#         copy(spt1, spt2)
#         if mod(it-1, samp) == 0
#            write(fid, spt1)
#         end
#         if mod(it, 500) == 0
#            println("the current step: $it")
#         end
#     end
#     close(fid)
#     return nothing
# end

# function MultiStepForward(pos::Array{Int64,2}, path::ASCIIString, fidMtx::FidMtx; tmax=1.0, interval=500, print_flag=false)
#     nz = fidMtx.nz ; nx    = fidMtx.nx;
#     ext= fidMtx.ext; iflag = fidMtx.iflag;
#     if iflag == 1
#        zupper = ext
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        zupper = 0
#        Nz = nx +   ext
#     end
#     Nx = nx + 2*ext
#     dt = fidMtx.dt
#     nt = round(Int64, tmax/dt+1)
#     ntsrc = numberOfSnapShots(path)
#     shot = InitShot(0, 0, pos, 0.0, nt, dt)
#     spt1 = readSnapShot(path, 1)
#     spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
#     spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, 1)
#     for it = 2: nt
#         OneStepForward!(spt2, spt1, fidMtx)
#         if it <= ntsrc
#            spt1 = readSnapShot(path, it)
#            spt2 = spt2 + spt1
#         end
#         copy(spt1, spt2)
#         spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
#         if print_flag
#            if mod(it, interval) == 0
#               println("the current step: $it")
#            end
#         end
#     end
#     return shot
# end
