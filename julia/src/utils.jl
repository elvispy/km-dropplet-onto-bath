function prev_index_for(numl::AbstractVector{<:Real}, idx::Int, value::Real)
    for offset in 0:(idx - 1)
        if numl[idx - offset] != value
            return idx - offset
        end
    end
    return idx
end

function normalize_indexes(indexes::Vector{Int})
    if length(indexes) >= 2 && indexes[2] == 1
        return indexes[2:end]
    end
    return indexes
end

function save_results(path::AbstractString, results::AbstractDict{String, <:Any})
    JLD2.jldopen(path, "w") do file
        for (k, v) in results
            file[k] = v
        end
    end
end
