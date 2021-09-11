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
function build_subproblem(data::SubProblemData, instance::ModelInstance)::SubProblem
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)

    set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_NET)

    _, n = size(data.A)

    @variable(model, x[i=1:n] >= 0)

    @constraint(model, data.A[1:instance.vertices,:] * x .== data.b[1:instance.vertices])
    @constraint(model, data.A[instance.vertices+1:end,:] * x .<= data.b[instance.vertices+1:end])

    return SubProblem(model, x)
end