function sptPlot(path::ASCIIString, it::Int64; name="NULL", attribute="p", xlabel="Length", xunits="(m)", ylabel="Depth", yunits="(m)", ox=0, dx=5, oy=0, dy=5, cmap="seismic", clip=1.0, dpi=100, aspect="auto", wbox=6, hbox=6, cbar=false, labelsize=15, ticksize=15)
    spt = readSnapShot(path, it)
    nz = spt.nz ; nx    = spt.nx   ;
    ext= spt.ext; iflag = spt.iflag;
    if iflag == 1
       Nz = nz + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    fig = plt.figure(figsize=(wbox, hbox), dpi=dpi, facecolor="w", edgecolor="k")
    if attribute == "p"
       d = reshape(spt.pz+spt.px, Nz, Nx)
    elseif attribute == "vz"
       d = reshape(spt.vz, Nz, Nx)
    elseif attribute == "vx"
       d = reshape(spt.vx, Nz, Nx)
    else
       error("wrong attribute")
    end
    vmax = quantile(vec(abs(d)), clip); vmin = -vmax
    im = plt.imshow(d, cmap=cmap, vmin=vmin, vmax=vmax, extent=[ox,ox+(size(d,2)-1)*dx,oy+(size(d,1)-1)*dy,oy], aspect=aspect)
    plt.xlabel(join([xlabel " " xunits]), fontsize=labelsize)
    plt.ylabel(join([ylabel " " yunits]), fontsize=labelsize)
    if cbar == true
      plt.colorbar()
    end
    ax = plt.gca()
    plt.setp(ax[:get_xticklabels](), fontsize=ticksize)
    plt.setp(ax[:get_yticklabels](), fontsize=ticksize)
    if (name == "NULL")
  		plt.show()
  	else
  		plt.savefig(name,dpi=dpi)
  		plt.close()
  	end
    return nothing
end
