function waveAnim(path::ASCIIString, pathout="NULL"; inter=5, attribute="p", wbox=6, hbox=6, ox=0.0, dx=5.0, oz=0.0, dz=5.0, clip=1.0, cmap="seismic", aspect="auto", interval=40)
    fid = open(path, "r")
    isz = read(fid, Int64); isx   = read(fid, Int64);
    nz  = read(fid, Int64); nx    = read(fid, Int64);
    ext = read(fid, Int64); iflag = read(fid, Int64);
    ot  = read(fid, Float64); dt  = read(fid, Float64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx  = nx + 2*ext
    head    = sizeof(Int64) * 8
    cptSize = sizeof(Float64) * Nz * Nx
    sptSize = cptSize * 4
    nt =  round(Int64, (filesize(fid)-head)/sptSize)
    fig = plt.figure(figsize=(wbox, hbox))
    ims = PyCall.PyObject[]
    for it = 1 : inter : nt
        if attribute == "vz"
           position = head + (it-1)*sptSize; seek(fid, position);
           d = reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        elseif attribute == "vx"
           position = head + (it-1)*sptSize + cptSize; seek(fid, position);
           d = reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        elseif attribute == "p"
           position = head + (it-1)*sptSize + cptSize*2; seek(fid, position);
           d = reshape(read(fid, Float64, Nz*Nx), Nz, Nx) + reshape(read(fid, Float64, Nz*Nx), Nz, Nx)
        else
           error("wrong attribute")
        end
        vmax = quantile(vec(abs(d)), clip)
        vmin = -vmax
        im = plt.imshow(d, cmap=cmap, vmin=vmin,vmax=vmax,extent=[ox, ox+(size(d,2)-1)*dx, oz+(size(d,1)-1)*dz,oz], aspect=aspect)
        push!(ims, PyCall.PyObject[im])
    end
    close(fid)
    ani = anim.ArtistAnimation(fig, ims, interval=interval, blit=true)
    pathout = join([pathout ".mp4"])
    ani[:save](pathout)
    return nothing
end
