
function extract_network_and_feeder(file_path::String; pattern = r"(?i)\\Network_(\d+)\\Feeder_(\d+)\\Master.dss")
    
    # Match the pattern against the file path
    matchedPattern = match(pattern, file_path)
    
    if matchedPattern !== nothing
        # Extract the network and feeder numbers
        Ntw = parse(Int, matchedPattern.captures[1])
        Fdr = parse(Int, matchedPattern.captures[2])
        return Ntw, Fdr
    else
        error("Pattern not found in the given file path")
    end
end


# extract network and feeder from jld file in the form joinpath(BASE_DIR,"test/data/MENWL/Ntw_$(ntw)_Fdr_$(fdr).jld2")

function extract_ntw_fdr_jld(file_path::String; pattern = r"(?i)\\Ntw_(\d+)_Fdr_(\d+).jld2")
    
    # Match the pattern against the file path
    matchedPattern = match(pattern, file_path)
    
    if matchedPattern !== nothing
        # Extract the network and feeder numbers
        Ntw = parse(Int, matchedPattern.captures[1])
        Fdr = parse(Int, matchedPattern.captures[2])
        return Ntw, Fdr
    else
        error("Pattern not found in the given file path")
    end
end


transBritishto4326 = Proj.Transformation("EPSG:27700", "EPSG:4326", always_xy=true) # British National Grid OSGB36 to WGS84
