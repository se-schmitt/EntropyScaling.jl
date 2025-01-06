using Documenter, DocumenterCitations, EntropyScaling

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(
    sitename="EntropyScaling.jl",
    format = Documenter.HTML(
        canonical = "https://se-schmitt.github.io/EntropyScaling.jl/",
        assets = ["assets/citations.css"],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Transport Properties" => "transport_properties.md",
        "Models" => [
            "models/ES_models.md",
            "models/CE_models.md"
        ],
        "References" => "references.md"
    ],
    plugins=[bib]
)

deploydocs(
    repo = "github.com/se-schmitt/EntropyScaling.jl.git",
)