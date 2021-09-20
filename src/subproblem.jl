"""
    SubProblemData

Representação dos dados necessários para a construção do subproblema.
"""
mutable struct SubProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
end

"""
    SubProblem

Representação do modelo do subproblema, onde `x` é uma referência para as suas variáveis.
"""
mutable struct SubProblem
    model::Model
    x::Vector{VariableRef}
end

"""
    build_subproblem(data::SubProblemData, instance::ModelInstance)::SubProblem

Contrói o subproblema do problema original utilizando uma repartição fornecida pelo usuário. O subproblema trata-se de um problema de fluxo.

# Argumentos

- `data`: Dados do subproblema na representação de `SubProblemData`.
- `instance`: Instância do problema, onde estão os dados gerais do problema linear, como a matriz de capacidade. 
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