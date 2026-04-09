export ESFramework

struct ESFrameworkParam{P,T} <: AbstractEntropyScalingParam{P}
    α0::CL.SingleParam{T}
    α1::CL.SingleParam{T}
    α2::CL.SingleParam{T}
    α3::CL.SingleParam{T}
    αln::CL.SingleParam{T}
    m::CL.SingleParam{T}
    Y₀⁺min::CL.SingleParam{T}
    ce::ChapmanEnskog
    prop::P
end

struct ESFrameworkDiffParam{T} <: AbstractEntropyScalingParam{DiffusionCoefficient}
    α0::CL.PairParam{T}
    α1::CL.PairParam{T}
    α2::CL.PairParam{T}
    α3::CL.PairParam{T}
    αln::CL.PairParam{T}
    m::CL.SingleParam{T}
    Y₀⁺min::CL.PairParam{T}
    ce::Matrix{ChapmanEnskog}
end
Base.length(param::P) where {P<:ESFrameworkDiffParam} = length(param.m)
transport_property(::Type{<:ESFrameworkDiffParam}) = DiffusionCoefficient()

struct ESFramework{E,P} <: ESFrameworkModel
    components::Vector{String}
    params::P
    eos::E
    sources::Vector{String}
end

"""
    ESFramework{E,P} <: ESFrameworkModel

    ESFramework(components, eos=nothing; userlocations=String[], collision_integral=KimMonroe(), verbose=false)
    ESFramework(components, eos, datasets::Vector{<:TransportPropertyData}; opts::FitOptions=FitOptions(), solute=nothing, verbose=false)

Entropy scaling framework [schmitt_entropy_2024,schmitt_entropy_2025](@cite).

The entropy scaling framework provides a physical way to model transport properties
(viscosity, thermal conductivity, diffusion coefficients) based on molecular-based EOS.
It enables fitting new models using only a few experimental data.

# Parameters

- `α0::CL.SingleParam`
- `α1::CL.SingleParam`
- `α2::CL.SingleParam`
- `α3::CL.SingleParam`
- `αln::CL.SingleParam`

# Example

```julia
using EntropyScaling, Clapeyron

# Load experimental sample data for n-butane
(T_exp, ϱ_exp, η_exp) = EntropyScaling.load_sample_data()
data = ViscosityData(T_exp, nothing, ϱ_exp, η_exp, :unknown)

# Create EOS model
eos_model = PCSAFT("butane")

# Create entropy scaling model (fit of parameters)
model = ESFramework(eos_model, [data])

# Calculation of viscosity at a state point
η = viscosity(model, 0.1e6, 300.)
```
"""
ESFramework

db_model_path(::Type{ESFramework}) = joinpath("ESFramework", "ESFramework_[PROP]_[EOS].csv")
const PARAMS_FRAMEWORK = ["α0","α1","α2","α3","αln"]
const REF_FRAMEWORK = ["10.1016/j.molliq.2023.123811","10.1038/s41467-025-57780-z"]

# Parameter function
function ESFramework(components, eos; userlocations=Dict(), collision_integral=KimMonroe(), verbose=false)
    _components = CL.format_components(components)
    
    params_dict = _get_empty_params_dict()
    for prop in [Viscosity(), ThermalConductivity(), DiffusionCoefficient()]
        _userlocations = get(userlocations, prop, String[])
        if prop == DiffusionCoefficient() && length(eos) > 1
            length(eos) == 1 && break
            filepaths = [get_db_path(ESFramework, prop, eos) for prop in [InfDiffusionCoefficient(), SelfDiffusionCoefficient()]]
            _params = CL.getparams(_components, filepaths; userlocations=_userlocations, asymmetricparams=PARAMS_FRAMEWORK, ignore_missing_singleparams=PARAMS_FRAMEWORK)
            for par in PARAMS_FRAMEWORK
                if !(par in keys(_params))
                    _params[par] = CL.PairParam(par, _components)
                end
            end
            _prop = prop
        else
            _prop = (prop == DiffusionCoefficient()) ? SelfDiffusionCoefficient() : prop
            filepaths = get_db_path(ESFramework, _prop, eos)
            _params = CL.getparams(_components, [filepaths]; userlocations=_userlocations, ignore_missing_singleparams=PARAMS_FRAMEWORK)
        end
        components_missing = [all(_v.ismissingvalues[i] for (_,_v) in _params) for i in eachindex(_components)]
    
        if any(components_missing)
            verbose && @info "No RefpropRES $(name(_prop)) parameters found for components: $(join(_components[components_missing],','))."
        else
            α0 = _params["α0"]
            if _prop == ThermalConductivity()
                α0.values[α0.ismissingvalues] .= 1.
            end
            α1 = _params["α1"]
            α2 = _params["α2"]
            α3 = _params["α3"]
            αln = _params["αln"]

            m = eos.params.segment

            if _prop isa DiffusionCoefficient
                N = length(eos)
                ce = Matrix{ChapmanEnskog}(undef,2,2)
                Y₀⁺min = CL.PairParam("minimum(Y₀⁺)", _components)
                for i in 1:N, j = 1:N 
                    if i == j
                        _eos = first(CL.split_model(eos, [i]))
                        _ce, _Y₀⁺min = init_framework_params(_eos, SelfDiffusionCoefficient(); collision_integral)
                    else
                        _ce, _Y₀⁺min = init_framework_params(eos, InfDiffusionCoefficient(i => j); collision_integral)
                    end
                    ce[i,j] = _ce
                    Y₀⁺min.values[i,j] = only(_Y₀⁺min)
                end
                param = ESFrameworkDiffParam(α0,α1,α2,α3,αln,m,Y₀⁺min,ce)
            else
                ce, _Y₀⁺min = init_framework_params(eos, _prop; collision_integral)
                Y₀⁺min = CL.SingleParam("minimum(Y₀⁺)", _components, _Y₀⁺min)
                param = ESFrameworkParam(α0,α1,α2,α3,αln,m,Y₀⁺min,ce,_prop)
            end

            params_dict[prop] = param
        end
    end

    params = ParamVector(params_dict)
    ismissing(params) && error("No parameters found for components: $(join(components, ',')).")

    return ESFramework(_components, params, eos, REF_FRAMEWORK)
end

# Fitting function
function ESFramework(components, eos, datasets::Vector{<:TransportPropertyData}; tofit=Dict(), collision_integral=KimMonroe(), verbose=false)
    _components = CL.format_components(components)

    params = ESFrameworkParam[]
    _dataprops = getproperty.(datasets, :prop) |> unique
    for prop in _dataprops
        if prop in [Viscosity(), ThermalConductivity(), SelfDiffusionCoefficient()]
            length(eos) != 1 && error("Only one component allowed for fitting.")
            _eos = eos
            _components = _components
        else
            _eos = split_model(eos)[prop.solvent]
            _components = ["$(_components[prop.solute]) in $(_components[prop.solvent])"]
        end
        data = collect_data(datasets, prop)

        _tofit = prop in keys(tofit) ? tofit[prop] : get_tofit(ESFramework, prop)

        ce, _Y₀⁺min = init_framework_params(eos, prop; collision_integral)
        α0, α1, α2, α3, αln = get_α0(prop, _components)

        m = _eos.params.segment
        Y₀⁺min = CL.SingleParam("minimum(Y₀⁺)", _components, _Y₀⁺min)

        param = ESFrameworkParam(α0, α1, α2, α3, αln, m, Y₀⁺min, ce, prop)

        for k in findall(isnan.(data.ϱ))
            data.ϱ[k] = molar_density(_eos, data.p[k], data.T[k])
        end

        s   = entropy_conf.(_eos, data.ϱ, data.T)
        sˢ  = scaling_variable.(param, s)
        Yˢ  = scaling.(param, _eos, data.Y, data.T, data.ϱ, s)

        f_scale(x) = prop isa ThermalConductivity ? x : log(x)
        fit_fun(xs, p) = begin
            _T = Base.promote_eltype(xs,p)
            _param = CL.promote_model_struct(_T, param)
            for (i,_n) in enumerate(_tofit)
                getproperty(_param, _n) .= p[i]
            end
            return f_scale.(scaling_model.(_param, xs))
        end
        sol = curve_fit(fit_fun, sˢ, f_scale.(Yˢ), randn(length(_tofit)); autodiff=:forwarddiff)
        α_fitted = coef(sol)
        for (i,_n) in enumerate(_tofit)
            getproperty(param, _n) .= α_fitted[i]
        end

        push!(params, param)
    end

    return ESFramework(_components, params, eos, REF_FRAMEWORK)
end

# Fitting utils
get_tofit(::Type{ESFramework}, prop::AbstractTransportProperty) = [:α1,:α2,:α3,:αln]
get_tofit(::Type{ESFramework}, prop::ThermalConductivity) = [:α0,:α1,:α2,:α3,:αln]

get_α0(::AbstractTransportProperty, components) = begin 
    N = length(components)
    return (CL.SingleParam(_n,components,zeros(N)) for _n in ("α0","α1","α2","α3","αln"))
end
get_α0(::ThermalConductivity, components) = begin 
    N = length(components)
    return (CL.SingleParam(_n,components,(_n == "α0") ? ones(N) : zeros(N)) for _n in ("α0","α1","α2","α3","αln"))
end

function init_framework_params(eos, prop; collision_integral)
    eos_pure = split_model(eos)
    cs = crit_pure.(eos_pure)
    (Tc, pc) = [getindex.(cs, i) for i in 1:2]
    σε = correspondence_principle.(Tc, pc)
    (σ, ε) = [getindex.(σε, i) for i in 1:2]

    Mw = get_Mw(eos)
    if typeof(prop) == InfDiffusionCoefficient
        idx_sol = [prop.solvent,prop.solute]
        σ = [mean(σ[idx_sol])]
        ε = [geomean(ε[idx_sol])]
        Mw = [calc_M_CE(Mw[idx_sol])]
        __components = get_components(eos)
        _components = ["$(__components[prop.solute]) in $(__components[prop.solvent])"]
    else
        _components = get_components(eos)
    end

    ce = ChapmanEnskog(_components; userlocations=(; sigma=σ, epsilon=ε, Mw=Mw), collision_integral)

    idx_Y = (prop isa InfDiffusionCoefficient) ? (prop.solvent:prop.solvent) : 1:length(eos)
    Y₀⁺min = zeros(length(idx_Y))
    for (i,idx) in enumerate(idx_Y)
        optf(x) = property_CE_plus(prop, ce, eos, x[1]; i=idx)
        sol = optimize(optf, [2*Tc[idx]], LBFGS(), Optim.Options(f_reltol=1e-8); autodiff=AutoForwardDiff())
        Y₀⁺min[i] = Optim.minimum(sol)[1]
    end

    return ce, Y₀⁺min
end

# Scaling model (correlation: Yˢ = Yˢ(sˢ,α,g))
function scaling_model(param::ESFrameworkParam{<:AbstractViscosity}, s, x=[1.])
    g = (-1.6386, 1.3923)
    return exp(generic_scaling_model(param, s, x, g))
end
function scaling_model(param::ESFrameworkParam{<:AbstractThermalConductivity}, s, x=[1.])
    g = (-1.9107, 1.0725)
    return generic_scaling_model(param, s, x, g)
end
function scaling_model(param::ESFrameworkParam{<:AbstractDiffusionCoefficient}, s, x=[1.])
    g = (0.6632, 9.4714)
    return exp(generic_scaling_model(param, s, x, g))
end
function scaling_model(param::ESFrameworkDiffParam, s, x=[1.])
    g = (0.6632, 9.4714)
    return exp.(generic_scaling_model(param, s, x, g))
end

function generic_scaling_model(param::ESFrameworkParam, s, x, g)
    num = _dot(param.α0,x) + _dot(param.αln,x)*log1p(s) + _dot(param.α1,x)*s + _dot(param.α2,x)*s^2 + _dot(param.α3,x)*s^3
    denom = 1 + g[1]*log1p(s) + g[2]*s
    return num / denom
end

function generic_scaling_model(param::ESFrameworkDiffParam, s, x, g)
    num = param.α0.values*x .+ (param.αln.values*x).*log1p(s) .+ (param.α1.values*x).*s .+ (param.α2.values*x).*s^2 .+ (param.α3.values*x).*s^3
    denom = 1 + g[1]*log1p(s) + g[2]*s
    return num ./ denom
end

#sigmoid function with bias
W(x, sₓ=0.5, κ=20.0) = 1.0/(1.0+exp(κ*(x-sₓ)))

function scaling(param::ESFrameworkParam{P,TT}, eos, Y, T, ϱ, s, z=Z1; inv=false) where {P,TT}
    k = !inv ? 1 : -1

    prop   = transport_property(param)
    Y₀⁺    = property_CE_plus(prop, param.ce, eos, T, z)
    Y₀⁺min = mix_CE(prop, param.ce, param.Y₀⁺min.values, z)

    sˢ  = scaling_variable(param, s, z)
    Ws  = W(sˢ)
    base = BaseParam(P(),param.ce.Mw)
    Yˢ  = (Ws/Y₀⁺ + (1-Ws)/Y₀⁺min)^k * plus_scaling(base, Y, T, ϱ, s, z; inv=inv)
    return Yˢ
end

function scaling_variable(param::Union{ESFrameworkParam, ESFrameworkDiffParam}, s, z=Z1)
    return -s / R / _dot(param.m, z)
end

function ϱT_self_diffusion_coefficient(model::ESFramework, ϱ, T, z::AbstractVector)
    params_diff = model.params[DiffusionCoefficient()]
    param = _init_selfdiff_param(params_diff)
    s  = entropy_conf(model.eos, ϱ, T, z)
    sˢ = scaling_variable(param, s, z)

    D = zero(z)
    for i in eachindex(z)
        _set_selfdiff_param!(param, params_diff, i)
        Dˢ = scaling_model(param, sˢ, z)
        D[i] = scaling(param, model.eos, Dˢ, T, ϱ, s, z; inv=true)
    end
    
    return D
end

_init_selfdiff_param(params::ESFrameworkDiffParam) = _init_param(params, SelfDiffusionCoefficient())
_init_msdiff_param(params::ESFrameworkDiffParam) = _init_param(params, MaxwellStefanDiffusionCoefficient())

function _init_param(params::ESFrameworkDiffParam, prop; idx=nothing)
    _idx = isnothing(idx) ? (1:length(params)) : idx
    N = length(_idx)
    components = params.α0.components[_idx]
    params_new = ESFrameworkParam(
        CL.SingleParam("α0", components, zeros(N)),
        CL.SingleParam("α1", components, zeros(N)),
        CL.SingleParam("α2", components, zeros(N)),
        CL.SingleParam("α3", components, zeros(N)),
        CL.SingleParam("αln", components, zeros(N)),
        CL.SingleParam(params.m.name, components, getindex.(Ref(params.m), _idx)),
        CL.SingleParam("Y₀⁺min", components, zeros(N)),
        ChapmanEnskog(
            components, 
            CL.SingleParam("sigma", components, zeros(N)),
            CL.SingleParam("epsilon", components, zeros(N)),
            CL.SingleParam("Mw", components, zeros(N)),
            params.ce[1,1].collision,
            String[]
        ),
        prop
    )
    return params_new
end

function _set_selfdiff_param!(param::ESFrameworkParam, param_all::ESFrameworkDiffParam, i)
    for k in (:α0, :α1, :α2, :α3, :αln, :Y₀⁺min)
        getproperty(param, k).values .= getproperty(param_all, k).values[i,:]
    end
    for k in (:Mw, :sigma, :epsilon)
        getproperty(param.ce, k).values .= first.(getproperty.(getproperty.(param_all.ce[i,:], k), :values))
    end
    return nothing
end

function _set_msdiff_param!(param::ESFrameworkParam, param_all::ESFrameworkDiffParam, i, j)
    for k in (:α0, :α1, :α2, :α3, :αln, :Y₀⁺min)
        vals = getproperty(param_all, k).values
        getproperty(param, k).values[1] = vals[j,i]
        getproperty(param, k).values[2] = vals[i,j]
    end
    for k in (:Mw, :sigma, :epsilon)
        vals = first.(getproperty.(getproperty.(param_all.ce, k), :values))
        getproperty(param.ce, k).values[1] = vals[j,i]
        getproperty(param.ce, k).values[2] = vals[j,i]
    end
    return nothing
end
