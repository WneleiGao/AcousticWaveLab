function OneStepAdjoint!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
    spt2.it = spt1.it - 1
    spt2.vz = spt1.vz + (fidMtx.MpzBvz)' * spt1.pz
    spt2.vx = spt1.vx + (fidMtx.MpxBvx)' * spt1.px
    spt2.pz =           (fidMtx.MpzBpz)' * spt1.pz
    spt2.px =           (fidMtx.MpxBpx)' * spt1.px

    spt2.pz = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.pz
    spt2.px = (fidMtx.MvzBp )' * spt2.vz + (fidMtx.MvxBp)' * spt2.vx + spt2.px
    spt2.vz = (fidMtx.MvzBvz)' * spt2.vz
    spt2.vx = (fidMtx.MvxBvx)' * spt2.vx
    return nothing
end

function MultiStepAdjoint(path::ASCIIString, shot::Shot, fidMtx::FidMtx; interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx   ;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; nt    = shot.nt     ; dt = shot.dt;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot, Nz, Nx, ext, zupper)
    path_tmp = join([path "_tmp"])
    fid  = writeSnapShot(path_tmp, spt1)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
        copy(spt1, spt2)
        writeSnapShot(fid, spt1)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    close(fid)
    reverseSnapShotsOrder(path, path_tmp)
    return nothing
end

function MultiStepAdjoint(ntw::Int64, path::ASCIIString, shot::Shot, fidMtx::FidMtx)
    nz = fidMtx.nz ; nx    = fidMtx.nx   ;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; nt    = shot.nt     ; dt = shot.dt;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot, Nz, Nx, ext, zupper)
    path_tmp = join([path "_tmp"])
    if ntw > nt
       error("length of src is larger than records")
    elseif ntw == nt
       fid = saveSpt2SrcSpt(path_tmp, spt1, nz, nx, ntw, ext, Nz, Nx, zupper)
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
           copy(spt1, spt2)
           saveSpt2SrcSpt(fid, spt1, nz, nx, ext, Nz, Nx, zupper)
       end
    elseif ntw < nt
       for it = nt-1: -1: ntw
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
           copy(spt1, spt2)
       end
       fid = saveSpt2SrcSpt(path_tmp, spt1, nz, nx, ntw, ext, Nz, Nx, zupper)
       for it = ntw-1 : -1 : 1
           OneStepAdjoint!(spt2, spt1, fidMtx)
           AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
           copy(spt1, spt2)
           saveSpt2SrcSpt(fid, spt1, nz, nx, ext, Nz, Nx, zupper)
       end
    end
    close(fid)
    reverseSrcSptOrder(path, path_tmp)
    return nothing
end

function MultiStepAdjoint(w::Array{Float64,1}, shot::Shot, fidMtx::FidMtx; interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx   ;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt    = shot.dt     ;
    nt = shot.nt   ; ntsrc = length(w)   ;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    dis= zeros(nz*nx)
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot, Nz, Nx, ext, zupper)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
        copy(spt1, spt2)
        if it <= ntsrc
           dis = dis + spt2dis(w, spt1)
        end
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return dis
end


function MultiStepAdjoint(ntw::Int64, dis::Array{Float64,1}, shot::Shot, fidMtx::FidMtx; interval=500, print_flag=false)
    nz = fidMtx.nz ; nx    = fidMtx.nx   ;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    ot = shot.ot   ; dt    = shot.dt     ;
    nt = shot.nt   ;
    if iflag == 1
       zupper = ext
       Nz = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       zupper = 0
    end
    Nx = nx + 2*ext
    w  = zeros(ntw)
    spt1 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    spt2 = InitSnapShot(0, 0, nz, nx, ext, iflag, dt, nt)
    AddShot2SnapShot!(spt1, shot, Nz, Nx, ext, zupper)
    for it = nt-1: -1: 1
        OneStepAdjoint!(spt2, spt1, fidMtx)
        AddShot2SnapShot!(spt2, shot, Nz, Nx, ext, zupper)
        copy(spt1, spt2)
        if it <= ntw
           spt2wlet(w, dis, spt1)
        end
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return w
end
