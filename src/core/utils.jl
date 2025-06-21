
#extract_keys (compares missing keys and shows a warning)
"""
    extra_keys(eng_data::Dict{String, Any}, expected_keys)
Checks if there are extra keys in the KeySet of the dictionary `eng_data` that are not in the list `keys` and prints them.

# Arguments
- `eng_data::Dict{String, Any}`: A dictionary containing the data to be checked.
- `expected_keys`: A list of keys that are expected to be in the dictionary.
- `show_keys=false`: A boolean indicating whether to print the expected keys, the keys in the dictionary, and the extra keys.

# Description
This function compares the keys in the dictionary `eng_data` with the list of `expected_keys` and prints a warning message if there are extra keys in the dictionary that are not in the list. If `show_keys` is set to `true`, the expected keys, the keys in the dictionary, and the extra keys are printed.

# Example
```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
extra_keys(eng, ["bus", "line", "linecode", "load", "voltage_source", "time_series", "conductor_ids", "name", "settings", "files", "data_model"])
```
"""
function extra_keys(eng_data::Dict{String, Any}, expected_keys; show_keys=false)
    eng_data_keys = keys(first(eng_data).second)
    extra_keys = setdiff(eng_data_keys,expected_keys)
    
    if length(extra_keys) > 0
        warning_text("Extra coulmns exist in the data model but are not displayed: $(extra_keys)")
    end

    if show_keys
        @show expected_keys
        @show eng_data_keys
        @show extra_keys
    end
end


"""
    diff_vectors(vec1::Vector{Float64}, vec2::Vector{Float64})

Prints the difference between two vectors element-wise.

# Arguments
- `vec1::Vector{Float64}`: The first vector.
- `vec2::Vector{Float64}`: The second vector.

# Example

diff_vectors([1.0, 2.0, 3.0], [1.0, 2.0, 4.0])


"""
function diff_vectors(vec1::Vector{Float64}, vec2::Vector{Float64})
    green = _CRN.Crayon(foreground = :green)
    red = _CRN.Crayon(foreground = :red)
    display(vec1)
    display(vec2)
    println("diff :")
    for (num1, num2) in zip(vec1, vec2)
        str1 = string(num1)
        str2 = string(num2)
    
        # Ensure the strings have the same length for comparison
        len = max(length(str1), length(str2))
        str1 = rpad(str1, len)
        str2 = rpad(str2, len)
    
        for (d1, d2) in zip(str1, str2)
            if d1 == d2
                print(green(string(d1)))
            else
                print(red(string(d2)))
            end
        end
        println()  # Newline after each comparison    
    end
end

# Exporting
export diff_vectors


function convert_keys_to_symbols(data)
    if isa(data, Dict)
        new_data = Dict{Symbol, Any}()
        for (key, value) in data
            new_data[Symbol(key)] = convert_keys_to_symbols(value)
        end
        return new_data
    elseif isa(data, Vector)
        return [convert_keys_to_symbols(item) for item in data]
    else
        return data
    end
end


# fix MetaGraphs set_indexing_prop! function


"""
    set_indexing_prop!(g, prop)
    set_indexing_prop!(g, v, prop, val)

Make property `prop` into an indexing property. If any values for this property
are already set, each vertex must have unique values. Optionally, set the index
`val` for vertex `v`. Any vertices without values will be set to a default
("(prop)(v)").
"""
function set_indexing_prop!(g::AbstractMetaGraph, prop::Symbol; exclude=nothing)
    in(prop, g.indices) && return g.indices
    index_values = [g.vprops[v][prop] for v in keys(g.vprops) if haskey(g.vprops[v], prop)]
    length(index_values) != length(union(index_values)) && error("Cannot make $prop an index, duplicate values detected")
    index_values = Set(index_values)

    g.metaindex[prop] = Dict{Any,Integer}()
    for v in vertices(g)
        if !haskey(g.vprops, v) || !haskey(g.vprops[v], prop)
            val = default_index_value(v, prop, index_values, exclude=exclude)
            set_prop!(g, v, prop, val)
        end
        g.metaindex[prop][g.vprops[v][prop]] = v
    end
    push!(g.indices, prop)
    return g.indices
end


"""
    default_index_value(v, prop, index_values; exclude=nothing)

Provides a default index value for a vertex if no value currently exists. The default is a string: "\$prop\$i" where `prop` is the property name and `i` is the vertex number. If some other vertex already has this name, a randomized string is generated (though the way it is generated is deterministic).
"""
function default_index_value(v::Integer, prop::Symbol, index_values::Set{}; exclude=nothing)
    val = string(prop) * string(v)
    if in(val, index_values) || val == exclude
        seed!(v + hash(prop))
        val = randstring()
        @warn("'$(string(prop))$v' is already in index, setting ':$prop' for vertex $v to $val")
    end
    return val
end



"""Print a matrix as a pretty table with highlighted diagonal and off-diagonal elements"""
function _pretty_diag_matrix(mat::Matrix)
    isa(mat, Array) || error("The $key is not a matrix")
            mx = deepcopy(mat)
            mx = hcat([range(1, size(mx,1))...], mx)
            pretty_table(mx,
                         header=[range(0, size(mx,2)-1)...],
                         highlighters=(highlight_row_label, highlight_diagonal, highlight_off_diagonal),
                         formatters    = ft_printf("%1.0f", 1:1),
                         )
end


function show_example(dict::Dict)
    return first(dict).second
end

polarize(complex::Complex; scale=1.0::Float64, rounding_digits = 4) = "$(round(abs(complex)*scale, digits=rounding_digits)) ∠ $(round(rad2deg(angle(complex)), digits=rounding_digits))"



function _phase_letter(phase::Int) 
    if phase == 0
        return "g"
    elseif phase == 1
        return "a"
    elseif phase == 2
        return "b"
    elseif phase == 3
        return "c"
    elseif phase == 4
        return "n"
    else
        error("Invalid phase value: $phase")
    end
end


function search_files(directory, file_name)
    files = []
    for (root, dirs, file) in walkdir(directory)
        for f in file
            if occursin(file_name, f)
                push!(files, joinpath(root, f))
            end
        end
    end
    return files
end

function list_files(directory)
    files = []
    for (root, dirs, file) in walkdir(directory)
        for f in file
            push!(files, joinpath(root, f))
        end
    end
    return files
end

# write a search for directory

function search_directories(directory, directory_name)
    directories = []
    for (root, dirs, file) in walkdir(directory)
        for d in dirs
            if occursin(directory_name, d)
                push!(directories, joinpath(root, d))
            end
        end
    end
    return directories
end
