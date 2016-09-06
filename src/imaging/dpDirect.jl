function dpDirect(shot::Shot, v::Float64, src::Source, dx::Float64, dt::Float64)
    nr = length(shot.irz)
    nt = size(shot.d, 1 )
    isz = src.iz; isx = src.ix;
    lsrc= length(src.p)
    damp= hanning(11)[2:6]
    if shot.irz[1] != isz
       error("source and receivers are not on same level")
    end
    dp  = ones(nt, nr)
    for ir = 1 : nr
        l  = abs((shot.irx[ir]-isx)*dx)
        it = ceil(Int64,l/(v*dt))+1+lsrc*2
        dp[1:it, ir] = zeros(it)
        dp[it+1:it+5, ir] = damp
    end
    return dp
end
