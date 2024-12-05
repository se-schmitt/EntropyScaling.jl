using Documenter, EntropyScaling

makedocs(
    sitename="EntropyScaling.jl",
    format = Documenter.HTML(
        canonical = "https://se-schmitt.github.io/EntropyScaling.jl/",
        # assets = ["assets/logo.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Transport Properties" => "transport_properties.md",
        "Models" => "models.md",
    ]
)

deploydocs(
    repo = "github.com/se-schmitt/EntropyScaling.jl.git",
)