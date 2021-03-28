"""
---
Absolute distance classifier
---
"""
function absClassifier(dictionary,constellation,signal)

    (_,ids) = findmin(
            abs.(transpose(constellation[:]) .- signal[:]),
            dims=2);

    chunks = map(p -> dictionary[constellation[p[2]]],ids);
    return chunks
end