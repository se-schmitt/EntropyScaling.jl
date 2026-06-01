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
    _eos = _build_eos(_components, eos)
    
    params_dict = _get_empty_params_dict()
    for prop in [Viscosity(), ThermalConductivity(), DiffusionCoefficient()]
        _userlocations = get(userlocations, prop, String[])
        if length(_eos) == 1 && prop == DiffusionCoefficient() && isempty(_userlocations)
            _userlocations = get(userlocations, SelfDiffusionCoefficient(), String[])
        end
        if prop == DiffusionCoefficient()
            _userlocations = !isempty(_userlocations) ? map(_ensure_matrix, _userlocations) : _userlocations
            filepaths = [get_db_path(ESFramework, prop, _eos) for prop in [InfDiffusionCoefficient(), SelfDiffusionCoefficient()]]
            _params = CL.getparams(_components, filepaths; userlocations=_userlocations, asymmetricparams=PARAMS_FRAMEWORK, ignore_missing_singleparams=PARAMS_FRAMEWORK)
            for par in PARAMS_FRAMEWORK
                if !(par in keys(_params))
                    _params[par] = CL.PairParam(par, _components)
                end
            end
            _prop = prop
        else
            _prop = (prop == DiffusionCoefficient()) ? SelfDiffusionCoefficient() : prop
            filepaths = get_db_path(ESFramework, _prop, _eos)
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

            m = _eos.params.segment

            ce, Y₀⁺min = init_framework_params(_eos, _prop; collision_integral)
            if _prop isa AbstractDiffusionCoefficient
                param = ESFrameworkDiffParam(α0,α1,α2,α3,αln,m,Y₀⁺min,ce)
            else
                param = ESFrameworkParam(α0,α1,α2,α3,αln,m,Y₀⁺min,ce,_prop)
            end

            params_dict[prop] = param
        end
    end

    params = ParamVector(params_dict)
    ismissing(params) && error("No parameters found for components: $(join(components, ',')).")

    return ESFramework(_components, params, _eos, REF_FRAMEWORK)
end

# Fitting function
function ESFramework(components, eos, datasets::Vector{<:TransportPropertyData}; tofit=Dict(), α_initial=nothing, collision_integral=KimMonroe(), verbose=false)
    _components = CL.format_components(components)
    _eos = _build_eos(_components, eos)

    params_dict = _get_empty_params_dict()
    _dataprops = getproperty.(datasets, :prop) |> unique
    D_nan = ones(Bool, length(_eos), length(_eos))
    for prop in _dataprops
        if prop in [Viscosity(), ThermalConductivity()]
            length(_eos) != 1 && error("Only one component allowed for fitting.")
            _eos_pure = _eos
            idx = 1
        elseif prop isa SelfDiffusionCoefficient
            idx = (prop.component, prop.component)
            _eos_pure = CL.split_model(_eos)[prop.component]
        elseif prop isa InfDiffusionCoefficient
            idx = (prop.solute, prop.solvent)
            _eos_pure = CL.split_model(_eos)[prop.solvent]
        end
        data = collect_data(datasets, prop)

        _tofit = prop in keys(tofit) ? tofit[prop] : get_tofit(ESFramework, prop)

        if prop isa AbstractDiffusionCoefficient
            _prop = DiffusionCoefficient()
            D_nan[idx...] = false
        else
            _prop = prop
        end

        if !(prop isa AbstractDiffusionCoefficient) || ismissing(params_dict[DiffusionCoefficient()])
            ce, Y₀⁺min = init_framework_params(_eos, prop; collision_integral)
            α0, α1, α2, α3, αln = get_α0(prop, _components)
            m = _eos.params.segment
            if prop isa AbstractDiffusionCoefficient
                param = ESFrameworkDiffParam(α0, α1, α2, α3, αln, m, Y₀⁺min, ce)
            else
                param = ESFrameworkParam(α0, α1, α2, α3, αln, m, Y₀⁺min, ce, prop)
            end
        else
            param = params_dict[DiffusionCoefficient()]
        end

        for k in findall(isnan.(data.ϱ))
            data.ϱ[k] = inv(CL.volume(_eos_pure, data.p[k], data.T[k]))
        end

        if prop isa InfDiffusionCoefficient
            _param = _init_param(param, prop; idx=prop.solvent:prop.solvent)
            _set_ij_diff_param!(_param, param, prop.solute, prop.solvent)
        elseif prop isa SelfDiffusionCoefficient
            _param = _init_param(param, prop; idx=prop.component:prop.component)
            _set_ij_diff_param!(_param, param, prop.component, prop.component)
        else
            _param = param
        end

        s = CL.VT_entropy_res.(_eos_pure, inv.(data.ϱ), data.T)
        sˢ  = scaling_variable.(_param, s)
        Yˢ = scaling.(_param, _eos_pure, data.Y, data.T, data.ϱ, s)

        f_scale(x) = prop isa ThermalConductivity ? x : log(x)
        fit_fun(xs, p) = begin
            _T = Base.promote_eltype(xs,p)
            __param = CL.promote_model_struct(_T, _param)
            for (i,_n) in enumerate(_tofit)
                getproperty(__param, _n) .= p[i]
            end
            return f_scale.(scaling_model.(__param, xs))
        end
        _α_initial = isnothing(α_initial) ? zeros(length(_tofit)) : α_initial
        sol = curve_fit(fit_fun, sˢ, f_scale.(Yˢ), _α_initial; autodiff=:forwarddiff)
        α_fitted = coef(sol)

        if !(prop isa AbstractDiffusionCoefficient) || ismissing(params_dict[DiffusionCoefficient()])
            for (i,_n) in enumerate(_tofit)
                getproperty(param, _n).values[idx...] = α_fitted[i]
            end
            
            params_dict[_prop] = param
        else
            for (i,_n) in enumerate(_tofit)
                getproperty(params_dict[_prop], _n).values[idx...] = α_fitted[i]
            end
        end
    end

    if !ismissing(params_dict[DiffusionCoefficient()])
        for _idx in findall(D_nan), _n in [:α0, :α1, :α2, :α3, :αln]
            getproperty(params_dict[DiffusionCoefficient()], _n).values[_idx] = NaN
        end
    end

    params = ParamVector(params_dict)

    return ESFramework(_components, params, _eos, REF_FRAMEWORK)
end

_ensure_matrix(x::AbstractMatrix) = x
_ensure_matrix(x::AbstractVector) = permutedims(x)

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
get_α0(::AbstractDiffusionCoefficient, components) = begin 
    N = length(components)
    return (CL.PairParam(_n,components,zeros(N,N)) for _n in ("α0","α1","α2","α3","αln"))
end

function init_framework_params(eos, prop; collision_integral)
    N = length(eos)

    eos_pure = CL.split_model(eos)
    cs = CL.crit_pure.(eos_pure)
    (Tc, pc) = (getindex.(cs, 1), getindex.(cs, 2))
    σε = correspondence_principle.(Tc, pc)
    (σ, ε) = (getindex.(σε, 1), getindex.(σε, 2))

    Mw = CL.mw(eos)
    _components = get_components(eos)
    
    ce = ChapmanEnskog(_components; userlocations=(; sigma=σ, epsilon=ε, Mw=Mw), collision_integral)

    _Y₀⁺min = zeros(N)
    for i in 1:N
        optf(x) = property_CE_plus(prop, ce, eos, x[1]; i)
        sol = optimize(optf, [2*Tc[i]], LBFGS(), Optim.Options(f_reltol=1e-8); autodiff=AutoForwardDiff())
        _Y₀⁺min[i] = Optim.minimum(sol)[1]
    end
    Y₀⁺min = CL.SingleParam("minimum(Y₀⁺)", _components, _Y₀⁺min)

    return ce, Y₀⁺min
end

function init_framework_params(eos, ::AbstractDiffusionCoefficient; collision_integral)
    N = length(eos)

    eos_pure = CL.split_model(eos)
    cs = CL.crit_pure.(eos_pure)
    (Tc, pc) = (getindex.(cs, 1), getindex.(cs, 2))
    σε = correspondence_principle.(Tc, pc)
    (σ, ε) = (getindex.(σε, 1), getindex.(σε, 2))

    Mw = CL.mw(eos)
    _components = get_components(eos)

    ce = Matrix{ChapmanEnskog}(undef,N,N)
    Y₀⁺min = CL.PairParam("minimum(Y₀⁺)", _components)

    for i in 1:N, j = 1:N 
        if i == j
            _ul = (; sigma=σ[i:j], epsilon=ε[i:j], Mw=Mw[i:j])
            _ce = ChapmanEnskog(_components[i]; userlocations=_ul, collision_integral)
            _prop = SelfDiffusionCoefficient(i)
            _eos = only(CL.split_model(eos, [i]))
        else
            _ul = (; sigma=[mean(σ[[i,j]])], epsilon=[geomean(ε[[i,j]])], Mw=[calc_M_CE(Mw[[j,i]])])
            _ce = ChapmanEnskog("$(_components[i]) in $(_components[j])"; userlocations=_ul, collision_integral)
            _prop = InfDiffusionCoefficient(i => j)
            _eos = eos
        end
        optf(x) = property_CE_plus(_prop, _ce, _eos, x[1]; i=j)
        sol = optimize(optf, [2*Tc[j]], LBFGS(), Optim.Options(f_reltol=1e-8); autodiff=AutoForwardDiff())
        ce[i,j] = _ce
        Y₀⁺min.values[i,j] = Optim.minimum(sol)[1]
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

function scaling(param::ESFrameworkParam{P,TT}, eos, Y, T, ϱ, s, z=Z1; inverse=false) where {P,TT}
    k = !inverse ? 1 : -1

    prop   = transport_property(param)
    Y₀⁺    = property_CE_plus(prop, param.ce, eos, T, z)
    Y₀⁺min = mix_CE(prop, param.ce, param.Y₀⁺min.values, z)

    sˢ  = scaling_variable(param, s, z)
    Ws  = W(sˢ)
    base = BaseParam(P(),param.ce.Mw)
    Yˢ  = (Ws/Y₀⁺ + (1-Ws)/Y₀⁺min)^k * plus_scaling(base, Y, T, ϱ, s, z; inverse)
    return Yˢ
end

function scaling_variable(param::Union{ESFrameworkParam, ESFrameworkDiffParam}, s, z=Z1)
    return -s / R / _dot(param.m, z)
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

function _set_ij_diff_param!(param::ESFrameworkParam, param_all::ESFrameworkDiffParam, i, j)
    for k in (:α0, :α1, :α2, :α3, :αln, :Y₀⁺min)
        getproperty(param, k).values[1] = getproperty(param_all, k).values[i,j]
    end
    for k in (:Mw, :sigma, :epsilon)
        getproperty(param.ce, k).values[1] = first(getproperty(getproperty(param_all.ce[i,j], k), :values))
    end
    return nothing
end
