function alterMini(w::Array{Float64,1}, shot::Shot, fidMtx::FidMtx, pos::Array{Int64,2}, mu::Float64, lambda::Float64; outer_it=20, fistait=10, CGit=5, tmax=0.5)
    nz = fidMtx.nz; nx = fidMtx.nx
    ntw = length(w)
    for k = 1 : outer_it
        println("outer loop: $k")
        dis = fista(shot, w, mu, lambda, fidMtx, pos_rec, tmax, fistait)
        (w, obj_w) = CG_w(ntw, shot, fidMtx, pos, dis, maxit=CGit, tmax=tmax)
    end
    return dis, w
end

function alterCG(w::Array{Float64,1}, shot::Shot, fidMtx::FidMtx, pos::Array{Int64,2}; outer_it=20, inner_it=5, tmax=0.5)
    nz = fidMtx.nz; nx = fidMtx.nx
    ntw = length(w)
    for k = 1 : outer_it
        println("outer loop: $k")
        (dis, obj_dis) = CG_dis(shot, fidMtx, pos, w, maxit=inner_it, tmax=tmax)
        (w  , obj_w  ) = CG_w(ntw, shot, fidMtx, pos, dis, maxit=inner_it, tmax=tmax)
    end
    return dis, w
end

function CG_w(ntw::Int64, shot::Shot, fidMtx::FidMtx, pos::Array{Int64,2}, dis::Array{Float64,1}; maxit=10, mu=0.1, tmax=0.5) # the default initial guess is zeros
    obj = zeros(maxit+1)
    x = zeros(ntw)
    r = copy(shot)
    s = MultiStepAdjoint(ntw, dis, r, fidMtx)
    p = copy(s)
    gamma = dot(s, s)
    obj[1] = gamma
    println("iteration 0, objective $gamma")
    for k = 1 : maxit
        q = MultiStepForward(pos, p, dis, fidMtx, tmax=tmax)
        delta = dot(vec(q.d), vec(q.d)) + mu*dot(p,p)
        alpha = gamma / delta
        x     = x + alpha*p
        r.d   = r.d - alpha*(q.d)
        s     = MultiStepAdjoint(ntw, dis, r, fidMtx) - mu*x
        gamma0= gamma
        gamma = dot(s,s)
        println("iteration $k, objective $gamma")
        obj[k+1] = gamma
        beta  = gamma / gamma0
        p     = s + beta*p
    end
    return x, obj
end

function CG_dis(shot::Shot, fidMtx::FidMtx, pos::Array{Int64,2}, w::Array{Float64,1}; maxit=10, mu=0.1, tmax=0.5) # the default initial guess is zeros
    obj = zeros(maxit+1)
    N = fidMtx.nz * fidMtx.nx
    x = zeros(N)
    r = copy(shot)
    s = MultiStepAdjoint(w, r, fidMtx)
    p = copy(s)
    gamma = dot(s, s)
    obj[1] = gamma
    println("iteration 0, objective $gamma")
    for k = 1 : maxit
        q = MultiStepForward(pos, w, p, fidMtx, tmax=tmax)
        delta = dot(vec(q.d), vec(q.d)) + mu*dot(p,p)
        alpha = gamma / delta
        x     = x + alpha*p
        r.d   = r.d - alpha*(q.d)
        s     = MultiStepAdjoint(w, r, fidMtx) - mu*x
        gamma0= gamma
        gamma = dot(s,s)
        println("iteration $k, objective $gamma")
        beta  = gamma / gamma0
        obj[k+1] = gamma
        p     = s + beta*p
    end
    return x, obj
end
