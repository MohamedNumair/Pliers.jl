using Pliers
using Documenter
using DocumenterVitepress

#DocMeta.setdocmeta!(Pliers, :DocTestSetup, :(using Pliers); recursive=true)  

makedocs(;
    modules = [Pliers, Pliers.PMDUtils, Pliers.PMDSEUtils, Pliers.PMDPlotting],
    repo = "https://github.com/MohamedNumair/Pliers.jl",
    authors = "Mohamed Numair <Mo7amednumair@gmail.com>",
    sitename = "Pliers.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/MohamedNumair/Pliers.jl",
        md_output_path = ".",    # comment when deploying
        build_vitepress = false, # comment when deploying
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    clean = false,
)

deploydocs(;
    repo = "https://github.com/MohamedNumair/Pliers.jl",
    target = "build", # this is where Vitepress stores its output
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)