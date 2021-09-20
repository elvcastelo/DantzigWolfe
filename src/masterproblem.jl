"""
    MasterProblemData

Representação dos dados do subproblema.
"""
mutable struct MasterProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
    V::Vector{Vector{Float64}}
end

"""
    MasterProblem

Modelo e referência para as variáveis do problema mestre. A variável λ_0 é referente à restrição de convexidade de α.
"""
mutable struct MasterProblem
    model::Model
    λ::Vector{VariableRef}
    λ_0::VariableRef
end

"""
    build_dualmaster(data::MasterProblemData)::MasterProblem

Constrói o problema mestre com os dados fornecidos em `data`.

# Argumentos

- `data`: Dados do problema mestre representados por `MasterProblemData`.
"""
function build_dualmaster(data::MasterProblemData)::MasterProblem
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)

    m, _ = size(data.A)
    p = length(data.V)

    @variables(model, begin
        λ[i=1:m-1] >= 0, Int
        λ_w <= 0, Int
        λ_0, Int
    end)

    variable_vector = [λ..., λ_w]

    @objective(model, Max, λ' * data.b + λ_0 + λ_w)

    @constraint(model, sum(λ) <= 2)
    @constraint(model, variable_vector' * data.A[:,2] <= 1)

    @inbounds for i = 1:p
        @constraint(model, variable_vector' * (data.A[:,3:end] * data.V[i]) + λ_0 <= 0)
    end

    return MasterProblem(model, variable_vector, λ_0)
end