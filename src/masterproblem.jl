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
    model = Model(CPLEX.Optimizer, bridge_constraints=false)
    set_silent(model)

    p = length(data.V)
    q = length(data.R)

    c_0 = data.c[1:2]
    c = data.c[3:end]

    A_0 = data.A[:,1:2]
    A = data.A[:,3:end]
    m, _ = size(A)

    @variables(model, begin
        z >= 0
        0 <= w <= 1
        α[i=1:p] >= 0
        β[i=1:q] >= 0
    end)

    constraint_vector = Vector{AffExpr}()
    objective_exp = c_0' * [z, w]

    @inbounds for j = 1:p
        add_to_expression!(objective_exp, (c' * data.V[j]) * α[j])
    end

    @inbounds for j = 1:q
        add_to_expression!(objective_exp, (c' * data.R[j]) * β[j])
    end

    @inbounds for i = 1:m
        constraint_exp = A_0[i,:]' * [z, w]

        @inbounds for j = 1:p
            add_to_expression!(constraint_exp, (A[i,:]' * data.V[j]) * α[j])
            add_to_expression!(objective_exp, (c' * data.V[j]) * α[j])
        end

        @inbounds for j = 1:q
            add_to_expression!(constraint_exp, (A[i,:]' * data.R[j]) * β[j])
            add_to_expression!(objective_exp, (c' * data.R[j]) * β[j])
        end

        push!(constraint_vector, constraint_exp)
    end

    @constraint(model, λ_0, sum(α) == 1)

    @objective(model, Min, objective_exp)

    return MasterProblem(model, λ, λ_0)
end