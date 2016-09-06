type PMLCoef
     vz::SparseMatrixCSC{Float64, Int64}
     vx::SparseMatrixCSC{Float64, Int64}
     pz::SparseMatrixCSC{Float64, Int64}
     px::SparseMatrixCSC{Float64, Int64}
end

# function InitCpml(nz::Int64, nx::Int64, ext::Int64, vmax::Float64, rho::Array{Float64,2}, iflag::Int64)
#     rho = modExpand(rho, ext, iflag)
#     if ext == 5
#        R = 0.01
#     elseif ext == 10
#        R = 0.001
#     elseif ext == 20
#        R = 0.0001
#     else
#        error("unsupported damping layers")
#     end
#     m = 0.25; n = 0.75; ll = 4;
#     if iflag == 1
#        Nz = nz + 2*ext
#        Nx = nx + 2*ext
#        vz = spzeros(Nz, Nx)
#        vx = spzeros(Nz, Nx)
#        pz = spzeros(Nz, Nx)
#        px = spzeros(Nz, Nx)
#        #  upper and left
#        for i = 1 : ext
#            vz[i,:] = -vmax/ext * log(R) * (m*(ext+1-i)/ext + n*((ext+1-i)/ext)^ll) * rho[i,:]
#            vx[:,i] = -vmax/ext * log(R) * (m*(ext+1-i)/ext + n*((ext+1-i)/ext)^ll) * rho[:,i]
#            pz[i,:] = -vmax/ext * log(R) * (m*(ext+1-i)/ext + n*((ext+1-i)/ext)^ll)
#            px[:,i] = -vmax/ext * log(R) * (m*(ext+1-i)/ext + n*((ext+1-i)/ext)^ll)
#        end
#        #  lower
#        for iz = nz+ext+1 : Nz
#            vz[iz,:] = -vmax/ext * log(R) * (m*(iz-nz-ext)/ext + n*((iz-nz-ext)/ext)^ll) * rho[iz,:]
#            pz[iz,:] = -vmax/ext * log(R) * (m*(iz-nz-ext)/ext + n*((iz-nz-ext)/ext)^ll)
#        end
#        #  right
#        for ix = nx+ext+1 : Nx
#            vx[:,ix] = -vmax/ext * log(R) * (m*(ix-nx-ext)/ext + n*((ix-nx-ext)/ext)^ll) * rho[:,ix]
#            px[:,ix] = -vmax/ext * log(R) * (m*(ix-nx-ext)/ext + n*((ix-nx-ext)/ext)^ll)
#        end
#     elseif iflag == 2
#        Nz = nz +   ext
#        Nx = nx + 2*ext
#        vz = spzeros(Nz, Nx)
#        vx = spzeros(Nz, Nx)
#        pz = spzeros(Nz, Nx)
#        px = spzeros(Nz, Nx)
#        #  lower
#        for iz = nz+1 : Nz
#            vz[iz,:] = -vmax/ext * log(R) * (m*(iz-nz)/ext + n*((iz-nz)/ext)^ll) * rho[iz,:]
#            pz[iz,:] = -vmax/ext * log(R) * (m*(iz-nz)/ext + n*((iz-nz)/ext)^ll)
#        end
#        #  left
#        for ix = 1 : ext
#            vx[:,ix] = -vmax/ext * log(R) * (m*(ext+1-ix)/ext + n*((ext+1-ix)/ext)^ll) * rho[:,ix]
#            px[:,ix] = -vmax/ext * log(R) * (m*(ext+1-ix)/ext + n*((ext+1-ix)/ext)^ll)
#        end
#        #  right
#        for ix = nx+ext+1 : Nx
#            px[:,ix] = -vmax/ext * log(R) * (m*(ix-nx-ext)/ext + n*((ix-nx-ext)/ext)^ll) * rho[:,ix]
#            px[:,ix] = -vmax/ext * log(R) * (m*(ix-nx-ext)/ext + n*((ix-nx-ext)/ext)^ll)
#        end
#     end
#     Cpml = PMLCoef(vz, vx, pz, px)
#     return Cpml
# end

function modExpand(par::Array{Float64,2}, ext::Int64, iflag::Int64)
    temp1 = repmat(par[:,  1], 1, ext)
    temp2 = repmat(par[:,end], 1, ext)
    par   = hcat(temp1, par, temp2)
    if iflag == 1
       temp1 = repmat(par[1,  :], ext, 1)
       temp2 = repmat(par[end,:], ext, 1)
       par   = vcat(temp1, par, temp2)
    elseif iflag == 2
       temp2 = repmat(par[end,:], ext, 1)
       par   = vcat(par, temp2)
    end
    return par
end

function modSmooth(par, L)
    n = 2*L + 1
    par1 = modExpand(par, L, 1)
    s = convert(Array{typeof(par1[1]), 1}, hanning(n))
    s = s / sum(s)
    par1 = conv2(s,s,par1)
    par1 = par1[2*L+1:end-2*L, 2*L+1:end-2*L]
    return par1
end


function InitCpml(nz::Int64, nx::Int64, ext::Int64, vmax::Float64, dz::Float64, dx::Float64, iflag::Int64)
      if ext == 10
         R = 0.01
      elseif ext == 20
         R = 0.001
      elseif ext == 30
         R = 0.0001
      else
         error("unsupported damping layers")
      end
      ll = 2;
      if iflag == 1
         Nz = nz + 2*ext
         Nx = nx + 2*ext
         vz = spzeros(Nz, Nx)
         vx = spzeros(Nz, Nx)
         pz = spzeros(Nz, Nx)
         px = spzeros(Nz, Nx)
         #  upper and left
         for i = 1 : ext
             vz[i,:] = -1.5 * vmax/(ext*dz) * log(R) *  ((ext+1-i)/ext)^ll
             vx[:,i] = -1.5 * vmax/(ext*dx) * log(R) *  ((ext+1-i)/ext)^ll
             pz[i,:] = -1.5 * vmax/(ext*dz) * log(R) *  ((ext+1-i)/ext)^ll
             px[:,i] = -1.5 * vmax/(ext*dx) * log(R) *  ((ext+1-i)/ext)^ll
         end
         #  lower
         for iz = nz+ext+1 : Nz
             vz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz-ext)/ext)^ll
             pz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz-ext)/ext)^ll
         end
         #  right
         for ix = nx+ext+1 : Nx
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
         end
      elseif iflag == 2
         Nz = nz +   ext
         Nx = nx + 2*ext
         vz = spzeros(Nz, Nx)
         vx = spzeros(Nz, Nx)
         pz = spzeros(Nz, Nx)
         px = spzeros(Nz, Nx)
         #  lower
         for iz = nz+1 : Nz
             vz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz)/ext)^ll
             pz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz)/ext)^ll
         end
         #  left
         for ix = 1 : ext
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ext+1-ix)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ext+1-ix)/ext)^ll
         end
         #  right
         for ix = nx+ext+1 : Nx
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
         end
      end
      Cpml = PMLCoef(vz, vx, pz, px)
      return Cpml
end
