using InvariantManifolds
using Documenter

DocMeta.setdocmeta!(InvariantManifolds, :DocTestSetup, :(using InvariantManifolds); recursive=true)

makedocs(;
    modules=[InvariantManifolds],
    authors="XiaomingZhang",
    sitename="InvariantManifolds.jl",
    warnonly = true,
    format=Documenter.HTML(;
        canonical="https://Xiaomingzzhang.github.io/InvariantManifolds.jl",
        edit_link="master",
        assets=String[],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1
    ),
    pages = [
        "主页"=>"index.md",
        "开始使用: 一维光滑流形"=>"smooth_one.md",
        "光滑两维流形"=>"smooth_two.md",
        "非光滑一维流形"=>"non_smooth_one.md",
        "非光滑两维流形"=>"non_smooth_two.md",
        "类型与函数"=>"api.md"
    ]
)

deploydocs(;
    repo="github.com/Xiaomingzzhang/InvariantManifolds.jl",
    devbranch = "master"
)
