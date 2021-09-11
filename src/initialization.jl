function get_initial_point(instance::ModelInstance, verbose::Bool)::Tuple{SubProblem, SubProblemData}
    restrictions_size = instance.vertices + instance.arcs
    A = zeros(restrictions_size, instance.arcs)
    b = zeros(restrictions_size)
    c = zeros(instance.arcs)

    b[1:instance.vertices] = instance.demands

    # Cria-se as restrições de fluxo
    @inbounds for i = 1:instance.vertices
        arcs_out = instance.arcs_indexes[i,:][findall(x->(x>0), instance.capacity_matrix[i,:])]
        arcs_in = instance.arcs_indexes[:,i][findall(x->(x>0), instance.capacity_matrix[:,i])]

        A[i, arcs_out] .= 1
        A[i, arcs_in] .= -1
    end

    # Cria-se restrições de capacidade
    @inbounds for i = instance.vertices+1:restrictions_size
        arc = i-instance.vertices
        A[i, arc] = 1

        u, v = instance.arcs_set[arc]
        b[i] = instance.capacity_matrix[u, v]
    end

    sub_data = SubProblemData(A, b, c)
    sub = build_subproblem(sub_data, instance)
    @objective(sub.model, Min, 0)

    if verbose; printstyled("[get_initial_point] Resolvendo subproblema para obter ponto extremo inicial.\n", color=:blue, bold=true) end
    optimize!(sub.model)

    return sub, sub_data
end

function initialize(instance::ModelInstance, verbose::Bool)
    # Obtêm um ponto viável inicial que satisfaz as restrições (2) e (5).
    sub, sub_data = get_initial_point(instance, verbose)

    if verbose; printstyled("[initialize] Obtendo dados para a formulação do problema mestre.\n", color=:blue, bold=true) end

    V = Vector{Vector{Float64}}()

    # Coloca o ponto adquirido em V
    push!(V, value.(sub.x)[1:instance.arcs])

    # Obtêm as demandas pares e ímpares
    even_demands = findall(instance.demands .% 2 .== 0)
    odd_demands = setdiff(1:instance.vertices, even_demands)

    A = zeros(Int, instance.vertices, instance.arcs+2)
    b = zeros(Int, instance.vertices)
    c = zeros(Int, instance.arcs+2)
    # Designa os coeficientes da função objetivo para z e w respectivamente
    c[1] = 2 ; c[2] = 1

    # Cria as restrições para as demandas pares e ímpares utilizando as definições do modelo
    @inbounds for vertex in even_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(instance.capacity_matrix[:,vertex] .> 0)]
        A[vertex, arcs.+2] .= -1
        A[vertex, 1] = 1
        b[vertex] = demand / 2
    end

    @inbounds for vertex in odd_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(instance.capacity_matrix[:,vertex] .> 0)]
        A[vertex, arcs.+2] .= -1
        A[vertex, [1, 2]] .= 1
        b[vertex] = cld(demand, 2)
    end

    last_constraint = zeros(Int, 1, instance.arcs+2)
    last_constraint[2] = 1
    A = vcat(A, last_constraint)

    master_data = MasterProblemData(A, b, c, V)
    
    # initialize_decomposition(masterproblem_data, subproblem_data, instance)
    return DantzigWolfe(master_data, sub, sub_data, verbose)
end

function originalModel(instance::ModelInstance)::Model
    model = Model(CPLEX.Optimizer, bridge_constraints=false)
    set_silent(model)

    @variables(model, begin
        z >= 0
        0<= w <= 1 
        ϕ[i=1:instance.arcs] >= 0
    end)

    even_demands = findall(x->(x%2==0), instance.demands)
    odd_demands = setdiff(1:instance.vertices, even_demands)

    # Cria as restrições para as demandas pares e ímpares utilizando as definições do modelo
    @inbounds for vertex in even_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(x->(x>0), instance.adjacency_matrix[:,vertex])]
        ϕ_neighborhood = ϕ[arcs]

        @constraint(model, z - sum(ϕ_neighborhood) >= demand / 2)
    end

    @inbounds for vertex in odd_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(x->(x>0), instance.adjacency_matrix[:,vertex])]
        ϕ_neighborhood = ϕ[arcs]

        @constraint(model, z + w - sum(ϕ_neighborhood) >= ceil(demand / 2))
    end

    @inbounds for i = 1:instance.vertices
        arcs_out = instance.arcs_indexes[i,:][findall(x->(x>0), instance.adjacency_matrix[i,:])]
        arcs_in = instance.arcs_indexes[:,i][findall(x->(x>0), instance.adjacency_matrix[:,i])]
        ϕ_out = ϕ[arcs_out]
        ϕ_in = ϕ[arcs_in]

        @constraint(model, sum(ϕ_out) - sum(ϕ_in) == instance.demands[i])
    end

    # Cria-se restrições de capacidade
    @inbounds for i = 1:instance.arcs
        arc = i
        u, v = instance.arcs_set[arc]
        
        @constraint(model, ϕ[i] <= instance.capacity_matrix[u, v])
    end

    @objective(model, Min, 2z + w)

    optimize!(model)

    return model
end