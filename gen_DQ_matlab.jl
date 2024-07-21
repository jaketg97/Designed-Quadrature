using MAT

function DQ_rules(d, p, n_s, mat)

    if isnothing(mat)
        # Run the MATLAB command from Julia
        run(`'/Applications/MATLAB_R2024a.app/bin/maca64/MATLAB' -nodisplay -nosplash -nodesktop -r "addpath(genpath('./matlab_code')); [XW,deltamain]=generator($d,$p,$n_s); save('matlab_code/DQ_rules.mat', 'XW'); exit;"`)

        # Load the output from MATLAB
        mat = matread("matlab_code/DQ_rules.mat")
    end
    
    XW = mat["XW"]
    # nodes are first d columns, weights are last column
    nodes = XW[:,1:d]
    weights = XW[:,d+1]

    return nodes, weights
end

# define DQ_integrate, which takes a function f, a dimension d, a polynomial degree p, and a number of samples n_s, and runs DQ_rules, then for each node and weight sum f(node) * weight
function DQ_integrate(f, d, p, n_s, mat)
    nodes, weights = DQ_rules(d, p, n_s, mat)
    return sum(f(nodes[i, axes(nodes, 2)]) * weights[i] for i in 1:size(nodes, 1))
end

# define f(x1,x2,x3,x4,x5) = x1^2 + x2^2 + x3^2 + x4^2 + x5^2
f(x::Vector{Float64}) = sum(xi for xi in x)
mat = matread("matlab_code/DQ_rules.mat")
DQ_integrate(f, 5, 5, 300, nothing)