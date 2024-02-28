using InvariantManifolds
using Documenter

DocMeta.setdocmeta!(InvariantManifolds, :DocTestSetup, :(using InvariantManifolds); recursive=true)

makedocs(;
    modules=[InvariantManifolds],
    authors="XiaomingZhang",
    sitename="InvariantManifolds.jl",
    format=Documenter.HTML(;
        canonical="https://Xiaomingzzhang.github.io/InvariantManifolds.jl",
        edit_link="master",
        assets=String[],
    ),
    pages = [
        "index.md",
        "state.md",
        "smooth.md",
        "nonsmooth.md",
        "api.md"
    ]
)

deploydocs(;
    repo="github.com/Xiaomingzzhang/InvariantManifolds.jl",
    devbranch="main",
)
