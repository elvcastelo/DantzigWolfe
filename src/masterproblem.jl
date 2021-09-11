mutable struct MasterProblemData
    A::Matrix{Int}
    b::Vector{Int}
    c::Vector{Int}
    V::Vector{Vector{Float64}}
end

mutable struct MasterProblemDual
    model::Model
    λ::Vector{VariableRef}
    λ_0::VariableRef
end

function build_dualmasterproblem(data::MasterProblemData)::MasterProblemDual
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

    return MasterProblemDual(model, variable_vector, λ_0)

end