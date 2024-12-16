export FrameworkModel, FrameworkParams

"""
    FrameworkParams{P<:AbstractTransportProperty,T<:Number}

Structure to store the parameters of the framework model. The parameters are:
- `α`: a matrix of size `(nparams,ncomponents)` containing the parameters of the corresponding transport property
- `m`: a vector of length `ncomponents` containing segment information
- `σ`: a vector of length `ncomponents` containing the molecular size parameters
- `ε`: a vector of length `ncomponents` containing the dispersion energies
- `Y₀⁺min`: vector of length `ncomponents` containing the minimum scaled property
- `base`: a `BaseParam` containing molecular weight, transport property and fitting information.
"""
struct FrameworkParams{P,T} <: AbstractEntropyScalingParams
    α::Matrix{T}
    m::Vector{Float64}
    Y₀⁺min::Vector{Float64}
    CE_model::ChapmanEnskogModel
    base::BaseParam{P}
end

# Constructor for fitting
function FrameworkParams(prop::AbstractTransportProperty, eos, data; solute=nothing)
    α0 = get_α0_framework(prop)
    Mw = get_Mw(eos)
    CE_model, Y₀⁺min = init_framework_params(eos, prop; Mw=Mw, solute=solute)
    base = BaseParam(prop, Mw, data; solute=solute)
    return FrameworkParams(α0, get_m(eos), Y₀⁺min, CE_model, base)
end

# Constructor for existing parameters
function FrameworkParams(prop::AbstractTransportProperty, eos, α::Array{T,2};
                         solute=nothing) where {T}
    size(α,1) == 5 || throw(DimensionMismatch("Parameter array 'α' must have 5 rows."))
    size(α,2) == length(eos) || throw(DimensionMismatch("Parameter array 'α' doesn't fit EOS model."))

    Mw = get_Mw(eos)
    CE_model, Y₀⁺min = init_framework_params(eos, prop; Mw=Mw, solute=solute)    
    base = BaseParam(prop, Mw)
    return FrameworkParams(α, convert(Vector{Float64},get_m(eos)), Y₀⁺min, CE_model, base)
end

# Constructor for merging multiple parameter sets
function FrameworkParams(self::FrameworkParams{<:SelfDiffusionCoefficient},
                         inf::FrameworkParams{<:InfDiffusionCoefficient}, idiff)

    what_inf = 1:length(self.base) .!= idiff
    new_self = deepcopy(self)
    new_self.α[:,what_inf] = inf.α[:,what_inf]
    for k in [:m,:Y₀⁺min]
        getfield(new_self,k)[what_inf] = getfield(inf,k)[what_inf]
    end
    for k in [:σ,:ε,:Mw]
        getfield(new_self.CE_model,k)[what_inf] = getfield(inf.CE_model,k)[what_inf]
    end
    # new_self.base.Mw[what_inf] = inf.base.Mw[what_inf] #TODO ensure used correctly
    
    return new_self
end

get_α0_framework(prop::Union{Viscosity,DiffusionCoefficient}) = zeros(Real,5,1)
get_α0_framework(prop) = [ones(Real,1);zeros(Real,4,1);]

function init_framework_params(eos, prop; Mw, solute = nothing)
    # Calculation of σ and ε
    eos_pure = vcat(split_model(eos),isnothing(solute) ? [] : solute)
    cs = crit_pure.(eos_pure)
    (Tc, pc) = [getindex.(cs,i) for i in 1:2]
    σε = correspondence_principle.(Tc,pc)
    (σ, ε) = [getindex.(σε,i) for i in 1:2]

    !isnothing(solute) ? append!(Mw,get_Mw(solute)) : nothing
    if typeof(prop) == InfDiffusionCoefficient
        length(eos_pure) != 2 && error("Solvent and solute must each contain one component.")
        _1 = ones(Float64, length(eos))
        σ = mean(σ)*_1
        ε = geomean(ε)*_1
        Mw = calc_M_CE(Mw)*_1
    end
    CE_model = ChapmanEnskogModel(repeat([""],length(eos)),σ,ε,Mw)

    # Calculation of Ymin
    Ymin = Vector{Float64}(undef,length(eos))
    for i in 1:length(eos)
        optf = OptimizationFunction((x,p) -> property_CE_plus(prop, CE_model, eos, x[1]; i=i), AutoForwardDiff())
        prob = OptimizationProblem(optf, [2*Tc[i]])
        sol = solve(prob, Optimization.LBFGS(), reltol=1e-8)
        Ymin[i] = sol.objective[1]
    end
    
    return CE_model, Ymin
end

"""
    FrameworkModel{T} <: AbstractEntropyScalingModel

A generic entropy scaling model.
"""
struct FrameworkModel{E,FP} <: AbstractEntropyScalingModel
    components::Vector{String}
    params::FP
    eos::E
end

@modelmethods FrameworkModel FrameworkParams

function cite_model(::FrameworkModel)
    print("Entropy Scaling Framework:\n---\n" *
          "(1) Schmitt, S.; Hasse, H.; Stephan, S. Entropy Scaling Framework for " *
          "Transport Properties Using Molecular-Based Equations of State. Journal of " *
          "Molecular Liquids 2024, 395, 123811. DOI: " *
          "https://doi.org/10.1016/j.molliq.2023.123811")
    return nothing
end


function FrameworkModel(eos, datasets::Vector{TPD}; opts::FitOptions=FitOptions(),
                        solute=nothing) where TPD <: TransportPropertyData
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
            what_fit = prop in keys(opts.what_fit) ? opts.what_fit[prop] : [false;ones(Bool,4)]
            param = FrameworkParams(prop, eos, data; solute=solute_)

            #TODO make this a generic function
            # Calculate density
            for k in findall(isnan.(data.ϱ))
                data.ϱ[k] = molar_density(eos, data.p[k], data.T[k])
            end

            # Scaling
            s = entropy_conf.(eos, data.ϱ, data.T)
            sˢ = scaling_variable.(param,s)
            Yˢ = scaling.(param, eos, data.Y, data.T, data.ϱ, s)

            # Fit
            function resid!(du, p, xy)
                (xs,ys) = xy
                param.α[what_fit] .= p
                du .= scaling_model.(param,xs) .- ys
                return nothing
            end
            Yˢ_fit = prop == ThermalConductivity() ? Yˢ : log.(Yˢ)
            prob = NonlinearLeastSquaresProblem(
                NonlinearFunction(resid!, resid_prototype=similar(Yˢ)),
                randn(sum(what_fit)), (sˢ, Yˢ_fit),
            )
            sol = solve(prob, SimpleGaussNewton(), reltol=1e-8)
            α_fit = get_α0_framework(prop)
            α_fit[what_fit] .= sol.u

            push!(params, FrameworkParams(float.(α_fit), param.m, param.Y₀⁺min, param.CE_model, param.base))
        end
    end
    return FrameworkModel(eos, params)
end


# Scaling model (correlation: Yˢ = Yˢ(sˢ,α,g))
function scaling_model(param::FrameworkParams{<:AbstractViscosity}, s, x=[1.])
    g = (-1.6386, 1.3923)
    return generic_scaling_model(param, s, x, g)
end
function scaling_model(param::FrameworkParams{<:AbstractThermalConductivity}, s, x=[1.])
    g = (-1.9107, 1.0725)
    return generic_scaling_model(param, s, x, g)
end
function scaling_model(param::FrameworkParams{<:DiffusionCoefficient}, s, x=[1.])
    g =  (0.6632, 9.4714)
    return generic_scaling_model(param, s, x, g)
end

function generic_scaling_model(param::FrameworkParams, s, x, g)
    α = param.α
    g1,g2 = g[1],g[2]
    num = zero(Base.promote_eltype(α,s,x,g))
    num += _dot(@view(α[1,:]),x) + _dot(@view(α[2,:]),x)*log1p(s)
    denom = 1 + g1*log1p(s) + g2*s
    si = s
    for i in 3:size(α,1)
        num += _dot(@view(α[i,:]),x)*si
        si *= s
    end
    return num/denom
end

#sigmoid function with bias
W(x, sₓ=0.5, κ=20.0) = 1.0/(1.0+exp(κ*(x-sₓ)))

function scaling(param::FrameworkParams, eos, Y, T, ϱ, s, z=[1.]; inv=false)
    k = !inv ? 1 : -1

    # Transport property scaling
    prop = transport_property(param)
    Y₀⁺ = property_CE_plus(prop, param.CE_model, eos, T, z)
    Y₀⁺min = mix_CE(prop, param.CE_model, param.Y₀⁺min, z)
    
    # Entropy
    sˢ = scaling_variable(param,s,z)
    Yˢ = (W(sˢ)/Y₀⁺ + (1.0-W(sˢ))/Y₀⁺min)^k * plus_scaling(param.base, Y, T, ϱ, s, z; inv=inv)
    return Yˢ
end

function scaling_variable(param::FrameworkParams, s, z = Z1)
    return -s / R / _dot(get_m(param),z)
end

#TODO generalize
function ϱT_self_diffusion_coefficient(model::FrameworkModel, ϱ, T, z)
    param_self = model[SelfDiffusionCoefficient()]
    param_inf = model[InfDiffusionCoefficient()]
    s = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param_self, s, z)
    Di = similar(z)

    for i in 1:length(model.eos)
        param = FrameworkParams(param_self, param_inf, i)
        Dˢ = exp(scaling_model(param, sˢ, z))
        Di[i] = scaling(param, model.eos, Dˢ, T, ϱ, s, z; inv=true)
    end

    return Di
end