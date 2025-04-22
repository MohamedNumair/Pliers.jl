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
