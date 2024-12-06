export cite_model

cite_model(::Any) = "No citation available for this model."

const SAMPLE_PATH = normpath(Base.pkgdir(EntropyScaling),"test","data")
function load_sample_data(;prop="eta")
    path = joinpath(SAMPLE_PATH,"exp_$prop.csv")
    dat = readdlm(path,',',skipstart=1)
    return [dat[:,i] for i in 1:3]
end