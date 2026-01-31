
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

"""
    convert_keys_to_symbols(data)

Recursively convert all dictionary keys from strings to symbols.

# Arguments
- `data`: Input data structure (Dict, Vector, or other types).

# Returns
- If `data` is a Dict: A new Dict with Symbol keys and recursively converted values.
- If `data` is a Vector: A new Vector with recursively converted elements.
- Otherwise: The original value unchanged.

# Examples
```julia
data = Dict("key1" => Dict("nested" => 1), "key2" => [Dict("a" => 2)])
result = convert_keys_to_symbols(data)
# result[:key1][:nested] == 1
```
"""
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





"""
    _pretty_diag_matrix(mat::Matrix)

Print a matrix as a pretty table with highlighted diagonal and off-diagonal elements.

# Arguments
- `mat::Matrix`: The matrix to display.

# Description
Displays the matrix with row indices prepended and highlights:
- Diagonal elements in cyan
- Off-diagonal elements in white
- Row labels in bold
"""
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

"""
    show_example(dict::Dict)

Return the value of the first entry in a dictionary.

# Arguments
- `dict::Dict`: A dictionary to inspect.

# Returns
The value associated with the first key in the dictionary.

# Examples
```julia
d = Dict("a" => 1, "b" => 2)
show_example(d)  # Returns 1 or 2 depending on iteration order
```
"""
function show_example(dict::Dict)
    return first(dict).second
end

"""
    polarize(complex::Complex; scale=1.0, rounding_digits=4)

Convert a complex number to polar notation string format.

# Arguments
- `complex::Complex`: The complex number to convert.

# Keyword Arguments
- `scale::Float64`: Scaling factor for the magnitude (default: 1.0).
- `rounding_digits::Int`: Number of decimal places for rounding (default: 4).

# Returns
A string in the format "magnitude ∠ angle_degrees".

# Examples
```julia
polarize(1 + 1im)  # "1.4142 ∠ 45.0"
polarize(1 + 1im; scale=1000)  # "1414.2136 ∠ 45.0"
```
"""
polarize(complex::Complex; scale=1.0::Float64, rounding_digits = 4) = "$(round(abs(complex)*scale, digits=rounding_digits)) ∠ $(round(rad2deg(angle(complex)), digits=rounding_digits))"


"""
    _phase_letter(phase::Int)

Convert a phase index to its corresponding letter designation.

# Arguments
- `phase::Int`: Phase index (0=ground, 1=a, 2=b, 3=c, 4=neutral).

# Returns
A single character string representing the phase.

# Throws
- `ErrorException`: If phase is not in range 0-4.
"""
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

"""
    search_files(directory, file_name)

Search recursively for files containing a specific string in their name.

# Arguments
- `directory`: Root directory to start the search.
- `file_name`: String to search for in file names.

# Returns
A vector of full paths to matching files.

# Examples
```julia
files = search_files("/path/to/dir", "test")
```
"""
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

"""
    list_files(directory)

List all files recursively in a directory.

# Arguments
- `directory`: Root directory to list files from.

# Returns
A vector of full paths to all files in the directory tree.

# Examples
```julia
all_files = list_files("/path/to/dir")
```
"""
function list_files(directory)
    files = []
    for (root, dirs, file) in walkdir(directory)
        for f in file
            push!(files, joinpath(root, f))
        end
    end
    return files
end

"""
    search_directories(directory, directory_name)

Search recursively for directories containing a specific string in their name.

# Arguments
- `directory`: Root directory to start the search.
- `directory_name`: String to search for in directory names.

# Returns
A vector of full paths to matching directories.

# Examples
```julia
dirs = search_directories("/path/to/dir", "test")
```
"""
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
