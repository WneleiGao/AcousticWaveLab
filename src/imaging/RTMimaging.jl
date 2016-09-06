# function RTMimaging(pathsrc, pathrec)
#     fidsrc = open(pathsrc, "r");
#     fidrec = open(pathrec, "r");
#     seek(fidsrc, sizeof(Float64)*2)
#     nz = read(fidsrc, Int64); nx    = read(fidsrc, Int64);
#     ext= read(fidsrc, Int64); iflag = read(fidsrc, Int64);
#     if iflag == 1
#        Nz = nz + 2*ext
#     elseif iflag == 2
#        Nz = nz +   ext
#     end
#     Nx = nz + 2*ext
#     head = sizeof(Int64)*8; sptsize = sizeof(Float64)*Nz*Nx*4
#     nt = round(Int64, (filesize(fidsrc)-head)/sptsize)
#     Ip = zeros(Nz, Nx)
#     head = sizeof(Int64)*8 + sizeof(Float64)*Nz*Nx*2
#     for it = 1 : nt
#         position = head + (it-1)*sptsize
#         seek(fidsrc, position); seek(fidrec, position)
#         psrc = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx) + reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
#         prec = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx) + reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
#         Ip = Ip + psrc .* prec
#     end
#     close(fidsrc); close(fidrec);
#     return Ip
# end

function RTMimaging(pathsrc, pathrec; image_type="p")
    fidsrc = open(pathsrc, "r");
    fidrec = open(pathrec, "r");
    seek(fidsrc, sizeof(Float64)*2)
    nz = read(fidsrc, Int64); nx    = read(fidsrc, Int64);
    ext= read(fidsrc, Int64); iflag = read(fidsrc, Int64);
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nz + 2*ext
    head = sizeof(Int64)*8; sptsize = sizeof(Float64)*Nz*Nx*4
    nt = round(Int64, (filesize(fidsrc)-head)/sptsize)
    if image_type == "p"
       Ip = zeros(Nz, Nx);
       head = sizeof(Int64)*8 + sizeof(Float64)*Nz*Nx*2
       for it = 1 : nt
           position = head + (it-1)*sptsize
           seek(fidsrc, position); seek(fidrec, position)
           psrc = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx) + reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
           prec = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx) + reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
           Ip = Ip + psrc .* prec
       end
       close(fidsrc); close(fidrec);
       return Ip
    elseif image_type == "all"
       Ivz = zeros(Nz, Nx); Ivx = zeros(Nz, Nx);
       Ipz = zeros(Nz, Nx); Ipx = zeros(Nz, Nx);
       for it = 1 : nt
           position = head + (it-1)*sptsize
           seek(fidsrc, position); seek(fidrec, position)
           Svz = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
           Svx = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
           Spz = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
           Spx = reshape(read(fidsrc, Float64, Nz*Nx), Nz, Nx)
           Rvz = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
           Rvx = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
           Rpz = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
           Rpx = reshape(read(fidrec, Float64, Nz*Nx), Nz, Nx)
           Ivz = Ivz + Svz .* Rvz
           Ivx = Ivx + Svx .* Rvx
           Ipz = Ipz + Spz .* Rpz
           Ipx = Ipx + Spx .* Rpx
       end
       close(fidsrc); close(fidrec);
       return Ivz, Ivx, Ipz, Ipx
    end
end
