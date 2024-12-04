export FrameworkModel, FrameworkParams

"""
    FrameworkParams

Structure to store the parameters of the framework model. The parameters are:
"""
struct FrameworkParams{T,P} <: AbstractEntropyScalingParams
    α::Array{T,2}
    m::Vector{Float64}
    σ::Vector{Float64}
    ε::Vector{Float64}
    Y₀⁺min::Vector{Float64}
    base::BaseParam{P}
end 

function FrameworkParams(prop::AbstractTransportProperty, eos, σ, ε, Y₀⁺min, data, what_fit; solute=nothing)
    what_fit = isempty(what_fit) ? [false;ones(Bool,4)] : what_fit
    α0 = get_α0_framework(prop)
    return FrameworkParams(α0, get_m(eos), σ, ε, Y₀⁺min, BaseParam(prop, eos, data, what_fit; solute=solute))
end

get_α0_framework(prop) = any(typeof(prop) .<: [Viscosity, DiffusionCoefficient]) ? zeros(Real,5,1) : [1.;zeros(Real,4,1);]

"""
    FrameworkModel{T} <: AbstractEntropyScalingModel

A generic entropy scaling model.
"""
struct FrameworkModel <: AbstractEntropyScalingModel
    components::Vector{String}
    params::Vector{FrameworkParams}
    info::EOSInfo
    references::Vector{Reference}
    FrameworkModel(eos,params) = new(get_components(eos),params,get_eos_info(eos))
end

function FrameworkModel(eos, datasets::Vector{TPD}; opts::FitOptions=FitOptions(), solute=nothing) where TPD <: TransportPropertyData
    # Check eos and solute
    length(eos) == 1 || error("Only one component allowed for fitting.")

    params = FrameworkParams[]
    for prop in [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient(), InfDiffusionCoefficient()]
        data = collect_data(datasets, prop)
        if data.N_dat > 0
            if typeof(prop) == InfDiffusionCoefficient 
                isnothing(solute) && error("Solute EOS model must be provided for diffusion coefficient at infinite dilution.")
                solute_ = solute
            else
                solute_ = nothing
            end

            # Init
            σ, ε, Ymin = init_framework_model(eos, prop; solute=solute_)
            what_fit = prop in keys(opts.what_fit) ? opts.what_fit[prop] : Bool[]
            param = FrameworkParams(prop, eos, σ, ε, Ymin, data, what_fit; solute=solute_)

            #TODO make this a generic function
            # Calculate density 
            for k in findall(isnan.(data.ϱ))
                data.ϱ[k] = molar_density(eos, data.p[k], data.T[k])
            end

            # Scaling
            s = entropy_conf.(eos, data.ϱ, data.T)
            sˢ = reduced_entropy.(Ref(param),s)
            Yˢ = scaling.(Ref(param), eos, data.Y, data.T, data.ϱ, s)

            # Fit
            function resid!(du, p, xy)
                (xs,ys) = xy
                param.α[param.base.what_fit] .= p
                du .= scaling_model.(Ref(param),xs) .- ys
                return nothing
            end
            Yˢ_fit = prop == ThermalConductivity() ? Yˢ : log.(Yˢ)
            prob = NonlinearLeastSquaresProblem(
                NonlinearFunction(resid!, resid_prototype=similar(Yˢ)),
                randn(sum(param.base.what_fit)), (sˢ, Yˢ_fit),
            )
            sol = solve(prob,reltol=1e-6)
            α_fit = get_α0_framework(prop)
            α_fit[param.base.what_fit] .= sol.u
            
            push!(params, FrameworkParams(float.(α_fit), param.m, param.σ, param.ε, param.Y₀⁺min, param.base))
        end
    end
    return FrameworkModel(eos, params)
end

get_prop_type(::FrameworkParams{T,P}) where {T, P <: AbstractTransportProperty} = P

function Base.getindex(model::FrameworkModel, prop::P) where P <: AbstractTransportProperty
    tprop = typeof(prop)
    k = get_prop_type.(model.params) .== tprop
    if sum(k) == 1
        i = findfirst(k)
    elseif sum(k) == 0
        @info "No parameters for $tprop."
        i = k
    else
        @info "Multiple parameters for $tprop. Returning the first one."
        i = findfirst(k)
    end
    return model.params[i]
end

function init_framework_model(eos, prop; solute=nothing, inv_dil=Int64[])
    # Calculation of σ and ε
    eos_pure = vcat(split_model(eos),isnothing(solute) ? [] : solute)
    cs = crit_pure.(eos_pure)
    (Tc, pc) = [getindex.(cs,i) for i in 1:2]
    σε = correspondence_principle.(Tc,pc)
    (σ, ε) = [getindex.(σε,i) for i in 1:2]

    if !isnothing(solute)
        length(eos_pure) != 2 && error("Solvent and solute must each contain one component.")
        σ = [mean(σ)]
        ε = [geomean(ε)]
    end

    for i in inv_dil
        σ[i] = mean(σ)
        ε[i] = geomean(ε)
    end

    # Calculation of Ymin
    Ymin = Vector{Float64}(undef,length(Tc))
    for i in 1:length(eos)
        optf = OptimizationFunction((x,p) -> property_CE_plus(prop,eos_pure[i], x[1], σ[i], ε[i]), AutoForwardDiff())
        prob = OptimizationProblem(optf, [2*Tc[i]])
        sol = solve(prob, Optimization.LBFGS(), reltol=1e-6)
        Ymin[i] = sol.objective[1]
    end

    return σ, ε, Ymin
end



# Scaling model (correlation: Yˢ = Yˢ(sˢ,α,g))
function scaling_model(param::FrameworkParams{T,Viscosity}, s, x=[1.]) where T
    g = [-1.6386, 1.3923]
    return generic_scaling_model(param, s, x, g)
end
function scaling_model(param::FrameworkParams{T,ThermalConductivity}, s, x=[1.]) where T
    g = [-1.9107, 1.0725]
    return generic_scaling_model(param, s, x, g)
end
function scaling_model(param::FrameworkParams{T,P}, s, x=[1.]) where {T, P<:DiffusionCoefficient}
    g = [ 0.6632, 9.4714]
    return generic_scaling_model(param, s, x, g)
end

function generic_scaling_model(param::FrameworkParams, s, x, g)
    α = param.α * x
    return (α' * [1., log(s+1.), s, s^2, s^3]) / (1. + g' * [log(s+1.), s])
end

function scaling(param::FrameworkParams, eos, Y, T, ϱ, s, z=[1.]; inv=false)
    k = !inv ? 1 : -1

    # Entropy
    sˢ = reduced_entropy(param,s)

    # Transport property scaling
    Y₀⁺ = sum(property_CE_plus.(Ref(param.base.prop), eos, T, param.σ, param.ε) .* z)
    Y₀⁺min = sum(param.Y₀⁺min .* z)

    W(x, sₓ=0.5, κ=20.0) = 1.0/(1.0+exp(κ*(x-sₓ)))

    Yˢ = (W(sˢ)/Y₀⁺ + (1.0-W(sˢ))/Y₀⁺min)^k * plus_scaling(param.base, Y, T, ϱ, s; inv=inv)

    return Yˢ
end

function reduced_entropy(param::FrameworkParams, s, z=[1.])
    return -s / R / sum(param.m .* z)
end


function viscosity(es_model::FrameworkModel, eos, p, T, z=[1.]; phase=:unknown)
    ϱ = molar_density(eos, p, T, z; phase=phase)
    return ϱT_viscosity(es_model, eos, T, ϱ, z)
end
function ϱT_viscosity(es_model::FrameworkModel, eos, ϱ, T, z=[1.])
    s = entropy_conf(eos, ϱ, T, z)
    sˢ = reduced_entropy(es_model.params[1],s)
    Yˢ = scaling(es_model.params[1], Y, T, ϱ, s)
end

function thermal_conductivity(es_model::FrameworkModel, eos, p, T, z=[1.]; phase=:unknown)
end

function self_diffusion_coefficient(es_model::FrameworkModel, eos, p, T, z=[1.]; phase=:unknown)
end

function MS_diffusion_coefficient(es_model::FrameworkModel, eos, p, T, z; phase=:unknown)
end