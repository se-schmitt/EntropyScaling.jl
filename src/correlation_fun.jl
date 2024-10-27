# Correlation function of the ES framework
function fun_es(s, p; prop::String)
    # Global parameters
    glob = merge(Dict(
        "vis" => [-1.6386, 1.3923],
        "tcn" => [-1.9107, 1.0725],),
        Dict(("dif","selfdif","mutdif") .=> Ref([ 0.6632, 9.4714])))
    g = glob[prop]

    # Compute the correlation function
    Y = (p[1] .+ p[2].*log.(s.+1.) .+ p[3].*s .+ p[4].*s.^2 .+ p[5].*s.^3) ./ (1. .+ g[1].*log.(s.+1.) .+ g[2].*s)

    return Y
end

# Transition function
W(s, sₓ=0.5, κ=20.0) = 1. ./ (1. .+ exp.(κ .* (s .- sₓ)))