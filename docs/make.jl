using Pliers
using Documenter
using DocumenterVitepress
using Literate
using Pkg

# 1. Automate Metadata from Project.toml
function get_project_toml()
    project_toml_path = joinpath(dirname(@__DIR__), "Project.toml")
    return Pkg.TOML.parsefile(project_toml_path)
end

project_info = get_project_toml()
PROJECT_VERSION = project_info["version"]
NAME = project_info["name"]
AUTHORS = join(project_info["authors"], ", ")

# 2. Automate Tutorials from examples/
# This converts .jl files in /examples to markdown in /docs/src/tutorials
tutorial_source = joinpath(dirname(@__DIR__), "examples")
tutorial_output = joinpath(@__DIR__, "src", "tutorials")

if isdir(tutorial_output)
    rm(tutorial_output, recursive=true)
end
mkpath(tutorial_output)

tutorial_pages = []
if isdir(tutorial_source)
    for file in readdir(tutorial_source)
        if endswith(file, ".jl")
            Literate.markdown(
                joinpath(tutorial_source, file),
                tutorial_output;
                documenter = true
            )
            # Add to menu
            clean_name = replace(file, ".jl" => "")
            push!(tutorial_pages, titlecase(clean_name) => joinpath("tutorials", "$clean_name.md"))
        end
    end
end

# 3. Automate TODO list
# Copies the root TODO.md to the docs source so it renders in the dev section
todo_src = joinpath(dirname(@__DIR__), "TODO.md")
todo_dest = joinpath(@__DIR__, "src", "TODO.md")
if isfile(todo_src)
    cp(todo_src, todo_dest; force=true)
end

# Define the pages structure
pages = [
    "Home" => "index.md",
    "Tutorials" => tutorial_pages,
    "API Reference" => "api.md",
    "Development" => [
        "To-Do List" => "TODO.md",
    ]
]

makedocs(;
    modules = [Pliers, Pliers.PMDUtils, Pliers.PMDSEUtils, Pliers.PMDPlotting],
    repo = "https://github.com/MohamedNumair/Pliers.jl",
    authors = AUTHORS,
    sitename = "$NAME.jl",
    version = PROJECT_VERSION,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/MohamedNumair/Pliers.jl",
        md_output_path = ".",    # comment when deploying
        build_vitepress = false, # comment when deploying
    ),
    pages = pages,
    clean = false,
    checkdocs = :exports, # Warn if exported functions are not documented
)

deploydocs(;
    repo = "github.com/MohamedNumair/Pliers.jl",
    target = "build",
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)