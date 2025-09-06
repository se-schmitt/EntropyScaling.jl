module EntropyScalingPlotsExt

using EntropyScaling
import Plots

const ES = EntropyScaling

function Plots.plot(model::ES.AbstractEntropyScalingModel, dat::ES.TransportPropertyData; 
                    slims=nothing, cprop=nothing, kwargs...)
    # Preprocessing 
    param = model[dat.prop]

    ϱdat = deepcopy(dat.ϱ)
    what_ϱ_nan = isnan.(ϱdat)
    ϱdat[what_ϱ_nan] = [ES.molar_density(model.eos, dat.p[i], dat.T[i]; phase=dat.phase[i]) 
                        for i in findall(what_ϱ_nan)]
    sdat = ES.entropy_conf.(model.eos, ϱdat, dat.T)
    sˢdat = ES.scaling_variable.(param, sdat)
    Yˢdat = ES.scaling.(param, model.eos, dat.Y, dat.T, ϱdat, sdat)

    slims = isnothing(slims) ? extrema(sˢdat) : slims
    sˢx = [range(slims..., length=100);]
    Yˢx = ES.scaling_model.(param, sˢx)

    # Settings
    cb_title = cprop == :T ? "T / K" : 
               cprop == :p ? "p / Pa" : 
               cprop == :ϱ ? "ϱ / mol/m³" : "$cprop"

    # Create the plot with general settings
    p = Plots.plot(;
        yscale=:log10,
        framestyle=:box,
        xlabel="sˢ",
        ylabel="$(ES.symbol(dat.prop))ˢ",
        xlims=slims,
        kwargs...
    )
    
    if !isnothing(cprop)
        p.attr[:colorbar_framestyle] = :box
        p.attr[:colorbar_title] = cb_title
    end

    # Plot data points
    if isnothing(cprop)
        Plots.scatter!(p, sˢdat, Yˢdat;
            msw=0,
            markersize=5,
            label=nothing,
            markercolor=:dodgerblue,
            z_order=1
        )
    else
        Plots.scatter!(p, sˢdat, Yˢdat;
            msw=0,
            markersize=5,
            label=nothing,
            zcolor=getfield(dat, cprop),
            z_order=1
        )
    end

    # Plot model
    Plots.plot!(p, sˢx, Yˢx;
        linecolor=:black,
        linewidth=1.5,
        label=nothing
    )
    
    return p
end

end # module