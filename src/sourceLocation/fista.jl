function power(maxit::Int64, fidMtx::FidMtx, w::Array{Float64,1}, pos::Array{Int64,2}, tmax::Float64)
    N = fidMtx.nz *fidMtx.nx
    x = randn(N)
    lambda = 0.0
    for k = 1 : maxit
        b = MultiStepForward(pos, w, x, fidMtx, tmax=tmax)
        y = MultiStepAdjoint(w, b, fidMtx)
        n = norm(x)
        x = y / n
        lambda = n
        println("iteration: $k, maximum eig: $lambda")
    end
    return lambda
end


function fista(shot::Shot, w::Array{Float64,1}, mu::Float64, lambda::Float64, fidMtx::FidMtx, pos::Array{Int64,2}, tmax::Float64, maxit::Int64)
    # J = zeros(maxit)
    N = fidMtx.nz *fidMtx.nx
    x = zeros(N)
    T = mu / (2*lambda)
    t = 1
    yk = copy(x)
    for k = 1 : maxit
        tmpx = copy(x)
        bshot = MultiStepForward(pos, w, x, fidMtx, tmax=tmax)
        bshot.d = bshot.d - shot.d
        x = MultiStepAdjoint(w, bshot, fidMtx)
        x = yk - x/lambda
        x = softThresh(x, T)
        tmpt = t
        t = (1 + sqrt(1+4*t^2)) / 2
        yk = x + (tmpt-1)/t * (x-tmpx)
        println("iteration $k")
    end
    return x
end


function softThresh(x::Array{Float64,1}, t::Float64)
    tmp = abs(x) - t
    tmp = (tmp + abs(tmp)) / 2
    y   = sign(x) .* tmp
    return y
end
