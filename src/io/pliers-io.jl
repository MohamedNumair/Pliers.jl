using Serialization

"""
    save_data(data, path::String)

Serializes and saves `data` to the specified `path`.
"""
function save_data(data, path::String)
    open(path, "w") do io
        serialize(io, data)
    end
end

"""
    read_data(path::String)

Reads and deserializes data from the specified `path`.
"""
function read_data(path::String)
    data = open(path, "r") do io
        deserialize(io)
    end
    return data
end


"""
    write_dat(x, filename; append=false)

Write an arbitrary Julia value `x` (commonly a dictionary) to `filename` (a `.dat` file).
- If `filename` does not end with `.dat` the extension is appended.
- Dictionaries and nested containers are flattened into key-path lines.
- Arrays/tuples are written with index notation (e.g. `arr[1] = 42`).
- Scalars and strings are written using `repr` so strings are quoted.

Examples:
    write_dat(mydict, "out.dat")
    write_dat(somearray, "array")  # writes to "array.dat"
"""
function write_dat(x, filename::AbstractString; append::Bool=false)
    fname = endswith(lowercase(filename), ".dat") ? filename : string(filename, ".dat")
    open_mode = append ? "a" : "w"
    open(fname, open_mode) do io
        write(io, "# write_dat output; root_type = $(typeof(x))\n")
        _write_entry(io, "", x)
    end
    return nothing
end

"""
    read_dat(filename)

Read a `.dat` file (formatted by `write_dat`) into a nested Dictionary.
- Lines starting with `#` are ignored.
- Keys like `a.b[1]` are parsed into nested dictionaries (e.g. `d["a"]["b"][1]`).
- Array indices in the file become integer keys in the dictionary.
- Values are parsed using `Meta.parse`.
"""
function read_dat(filename::AbstractString)
    data = Dict{Any,Any}()
    fname = endswith(lowercase(filename), ".dat") ? filename : string(filename, ".dat")
    
    isfile(fname) || error("File not found: $fname")

    for line in eachline(fname)
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end
        
        parts = split(line, "=", limit=2)
        length(parts) != 2 && continue
        
        key_str = strip(parts[1])
        val_str = strip(parts[2])
        
        val = try
            Meta.parse(val_str)
        catch
            val_str 
        end
        
        _insert_entry!(data, key_str, val)
    end
    return data
end

function _insert_entry!(root, key_path, value)
    cursor = root
    tokens = Any[]
    # Match property names or array indices
    for m in eachmatch(r"([^\.\[\]]+)|\[(\d+)\]", key_path)
        if m.captures[1] !== nothing
            push!(tokens, String(m.captures[1]))
        elseif m.captures[2] !== nothing
            push!(tokens, parse(Int, m.captures[2]))
        end
    end
    
    for (i, token) in enumerate(tokens)
        if i == length(tokens)
            cursor[token] = value
        else
            if !haskey(cursor, token)
                cursor[token] = Dict{Any,Any}()
            end
            cursor = cursor[token]
        end
    end
end

# internal helper: flatten nested structures into lines
function _write_entry(io::IO, prefix::AbstractString, x)
    if x isa AbstractDict
        for (k, v) in x
            keystr = _safe_key(k)
            newpref = isempty(prefix) ? keystr : string(prefix, ".", keystr)
            _write_entry(io, newpref, v)
        end
    elseif x isa AbstractArray || x isa Tuple
        for (i, v) in enumerate(x)
            newpref = isempty(prefix) ? string("[", i, "]") : string(prefix, "[", i, "]")
            _write_entry(io, newpref, v)
        end
    else
        name = isempty(prefix) ? "value" : prefix
        # use repr to produce valid Julia-like literal for strings/numbers
        write(io, string(name, " = ", repr(x), "\n"))
    end
    return nothing
end

# convert arbitrary dict key to a safe string for a key path
_safe_key(k) = k isa Symbol ? string(k) : k isa AbstractString ? k : replace(string(k), r"[.\[\]]" => "_")