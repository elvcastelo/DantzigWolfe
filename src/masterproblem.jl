mutable struct MasterProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
    V::Vector{Vector{Float64}}
end

mutable struct MasterProblem
    model::Model
    λ::Vector{ConstraintRef}
    λ_0::ConstraintRef
end

mutable struct MasterProblemDual
    model::Model
    λ::Vector{VariableRef}
    λ_0::VariableRef
end

function build_dualmasterproblem(data::MasterProblemData)::MasterProblemDual
    model = Model(CPLEX.Optimizer)
    # set_silent(model)

    m, _ = size(data.A)
    p = length(data.V)

    @variables(model, begin
        λ[i=1:m-1] >= 0
        λ_w <= 0
        λ_0
    end)

    variable_vector = [λ..., λ_w]

    @objective(model, Max, λ' * data.b + λ_0 + λ_w)

    @constraint(model, sum(λ) <= 2)
    @constraint(model, variable_vector' * data.A[:,2] <= 1)

    @inbounds for i = 1:p
        @constraint(model, variable_vector' * (data.A[:,3:end] * data.V[i]) + λ_0 <= 0)
    end

    return MasterProblemDual(model, variable_vector, λ_0)

end

function build_masterproblem(data::MasterProblemData)::MasterProblem
    model = Model(CPLEX.Optimizer, bridge_constraints=false)
    set_silent(model)

    p = length(data.V)

    c_0 = data.c[1:2]

    A_0 = data.A[:,1:2]
    A = data.A[:,3:end]
    m, _ = size(A)

    @variables(model, begin
        z >= 0
        0 <= w <= 1
        α[i=1:p] >= 0
    end)
    
    constraint_vector = Vector{AffExpr}()
    objective_exp = c_0' * [z, w]

    @inbounds for i = 1:m
        constraint_exp = A_0[i,:]' * [z, w]

        @inbounds for j = 1:p
            add_to_expression!(constraint_exp, (A[i,:]' * data.V[j]) * α[j])
        end

        push!(constraint_vector, constraint_exp)
    end

    @constraint(model, λ, constraint_vector .>= data.b)
    @constraint(model, λ_0, sum(α) == 1)

    @objective(model, Min, objective_exp)

    return MasterProblem(model, λ, λ_0)
end