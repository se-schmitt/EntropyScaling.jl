export cite_model

cite_model(::Any) = "No citation available for this model."

const SAMPLE_PATH = normpath(Base.pkgdir(EntropyScaling),"test","data")

function load_sample_data(;prop="eta")
    if prop == "eta"
        @info "Experimental data for the viscosity of n-butane.\n" *
              "\tUnits: [T] = K, [ϱ] = mol/m³, [η] = Pa·s"
    elseif prop == "lambda"
        @info "Experimental data for the thermal conductivity of n-butane.\n" *
              "\tUnits: [T] = K, [ϱ] = mol/m³, [λ] = W/m·K"
    elseif prop == "D"
        @info "Experimental data for the self-diffusion coefficient of n-butane.\n" *
              "\tUnits: [T] = K, [ϱ] = mol/m³, [D] = m²/s"
    else
        error("Invalid keyword value: prop = $prop")
    end
    path = joinpath(SAMPLE_PATH,"exp_$prop.csv")
    dat = readdlm(path,',',skipstart=1)
    return [dat[:,i] for i in 1:3]
end

function EStoJD end
export EStoJD
