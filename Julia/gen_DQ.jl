using Distributions
using LatinHypercubeSampling

function total_degree_indices(n::Int, d::Int)
    # Helper function to generate all combinations recursively
    function generate_combinations(sum::Int, depth::Int)
        if depth == 0
            return sum == 0 ? [[]] : []
        else
            combinations = []
            for i = 0:sum
                for tail in generate_combinations(sum - i, depth - 1)
                    push!(combinations, [i; tail])
                end
            end
            return combinations
        end
    end

    # Generate combinations for each degree up to n and flatten the result
    indices = []
    for degree in 0:n
        append!(indices, generate_combinations(degree, d))
    end

    return hcat(indices...)'
end

d, p = 2, 2
aind = total_degree_indices(d, p)

# Assuming `aind` is previously defined or obtained from a function like `hyperbolic_cross_indices(d, p)`
# aind = hyperbolic_cross_indices(d, p)

# Number of terms for total order
n_terms = size(aind, 1)

# Initialization for Gaussian
n_s, d = size(aind) # Assuming `n_s` and `d` are defined or can be inferred from `aind`
# Generate Latin Hypercube Samples
lhs_samples = sample(LatinHypercube(n_s), d)
# Apply the inverse normal distribution function to each sample
b = [quantile(Normal(), lhs_samples[i, j]) for i in 1:size(lhs_samples, 1), j in 1:size(lhs_samples, 2)]

rad = vecnorm(b .- 1, 2, 2) # Adjusted for Julia's 1-based indexing and function names
# For high d centralize the distances since the norm of nodes are high i.e., cent=1 otherwise cent=0;
cent = 1
w = exp.(-(rad .- cent .* mean(rad)).^2 ./ 2)
w = (w .* n_terms ./ sum(w))'

delta = 1; count = 1
xnew = b[:, 1]
for i = 2:d
    xnew = vcat(xnew, b[:, i])
end
xnew = vcat(xnew, w)

RHS = zeros(n_terms + n_s * (d + 1), 1); RHS[1] = 1