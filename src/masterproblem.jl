mutable struct MasterProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
    V::Vector{Vector{Float64}}
    R::Vector{Vector{Float64}}
end

mutable struct MasterProblem
    model::Model
    λ::Vector{ConstraintRef}
    λ_0::ConstraintRef
end

function build_masterproblem(data::MasterProblemData)::MasterProblem
    model = Model(CPLEX.Optimizer)
    
    p = length(data.V)
    q = length(data.R)

    c_0 = data.c[1:2]
    c = data.c[3:end]

    A_0 = data.A[:,1:2]
    A = data.A[:,3:end]

    @variables(model, begin
        z >= 0
        0 <= w <= 1
        α[i=1:p] >= 0
        β[i=1:q] >= 0
    end)

    constraint_exp = A_0 * [z, w]
    objective_exp = c_0'*[z, w]

    @inbounds for i = 1:p
        constraint_exp += (A * data.V[i]) * α[i]
        objective_exp += (c' * data.V[i]) * α[i]
    end

    @inbounds for j = 1:q
        constraint_exp += (A * data.R[j]) * β[j]
        objective_exp += (c' * data.R[i]) * β[j]
    end

    @constraint(model, λ, constraint_exp .>= data.b)
    @constraint(model, λ_0, sum(α) == 1)

    @objective(model, Min, objective_exp)

    return MasterProblem(model, λ, λ_0)
end