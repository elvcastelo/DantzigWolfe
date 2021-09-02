mutable struct SubProblem
    model::Model
    x::Vector{VariableRef}
end

mutable struct SubProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
end

"""
    build_subproblem(A::Matrix{Int}, A_line::Matrix{Int}, b::Vector{Int}, c::Vector{Int}, λ::Vector{Float64}, λ_0::Float64)

Contrói o subproblema do problema original utilizando uma repartição fornecida pelo usuário.
"""
function build_subproblem(data::SubProblemData, A_line::Matrix{Int}, λ::Vector, λ_0::Number, instance::ModelInstance, initial_solution=false)::SubProblem
    model = Model(CPLEX.Optimizer, bridge_constraints=false)
    set_silent(model)

    set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_NET)

    _, n = size(data.A)

    @variable(model, x[i=1:n] >= 0)
    
    if !initial_solution
        @objective(model, Min, (data.c' - λ' * A_line) * x - λ_0)
        @constraint(model, data.A[1:instance.vertices,:] * x .== data.b[1:instance.vertices])
        @constraint(model, data.A[instance.vertices+1:end,:] * x .<= data.b[instance.vertices+1:end])
    else
        @constraint(model, data.A * x .== data.b)
        @objective(model, Min, data.c' * x)
    end

    return SubProblem(model, x)
end