function get_initial_point(instance::ModelInstance)
    restrictions_size = instance.vertices + instance.arcs
    A = zeros(restrictions_size, instance.arcs)
    b = zeros(restrictions_size)

    # Criamos um vetor de custo com 1s nas posições das variáveis artificiais.
    c = zeros(2instance.arcs + instance.vertices)
    c[2instance.arcs+1:end] .= 1

    # Adicionamos as variáveis de folga para as restrições de capacidade. Pelo modelo são |A| restrições de capacidade.
    slack_matrix = Matrix(1.0I, instance.arcs, instance.arcs)
    slack_matrix = vcat(zeros(instance.vertices, instance.arcs), slack_matrix)

    # Adicionamos as variáveis artificiais para as restrições de fluxo. Pelo modelo são |V| restrições de fluxo.
    artificial_matrix = Matrix(1.0I, instance.vertices, instance.vertices)
    artificial_matrix = vcat(artificial_matrix, zeros(instance.arcs, instance.vertices))

    A = hcat(A, slack_matrix, artificial_matrix)
    b[1:instance.vertices] = instance.demands

    # Cria-se as restrições de fluxo
    @inbounds for i = 1:instance.vertices
        arcs_out = instance.arcs_indexes[i,:][findall(x->(x>0), instance.adjacency_matrix[i,:])]
        arcs_in = instance.arcs_indexes[:,i][findall(x->(x>0), instance.adjacency_matrix[:,i])]

        A[i, arcs_out] .= 1
        A[i, arcs_in] .= -1
    end

    # Cria-se restrições de capacidade
    for i = instance.vertices+1:restrictions_size
        arc = i-instance.vertices
        A[i, arc] = 1

        u, v = instance.arcs_set[arc]
        b[i] = instance.capacity_matrix[u, v]
    end

    subproblem_data = SubProblemData(A, b, c)
    subproblem = build_subproblem(subproblem_data, [0 0], [0], 0, instance, true)

    subproblem_data.A = A[:,1:instance.arcs]
    subproblem_data.c = c[1:instance.arcs]

    printstyled("[get_initial_point] Resolvendo subproblema para obter ponto extremo inicial.\n", color=:blue, bold=true)
    optimize!(subproblem.model)

    return subproblem, subproblem_data
end

function initialize(instance::ModelInstance)
    # Obtêm um ponto viável inicial que satisfaz as restrições (2) e (5).
    subproblem, subproblem_data = get_initial_point(instance)
    # subproblem, subproblem_data = get_initial_point(instance)
    V = Vector{Vector{Float64}}()
    R = Vector{Vector{Float64}}()

    # Coloca o ponto adquirido em V
    push!(V, value.(subproblem.x)[1:instance.arcs])

    # Obtêm as demandas pares e ímpares
    even_demands = findall(x->(x%2==0), instance.demands)
    odd_demands = setdiff(1:instance.vertices, even_demands)

    A = zeros(Int, instance.vertices, instance.arcs+2)
    b = zeros(Int, instance.vertices)
    c = zeros(Int, instance.arcs+2)
    c[1] = 2 ; c[2] = 1

    # Cria as restrições para as demandas pares e ímpares utilizando as definições do modelo
    @inbounds for vertex in even_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(x->(x>0), instance.adjacency_matrix[:,vertex])]
        A[vertex, arcs.+2] .= -1
        A[vertex, 1] = 1
        b[vertex] = demand / 2
    end

    @inbounds for vertex in odd_demands
        demand = instance.demands[vertex]
        arcs = instance.arcs_indexes[:,vertex][findall(x->(x>0), instance.adjacency_matrix[:,vertex])]
        A[vertex, arcs.+2] .= -1
        A[vertex, [1, 2]] .= 1
        b[vertex] = ceil(Int, demand / 2)
    end

    masterproblem_data = MasterProblemData(A, b, c, V, R)

    return DantzigWolfe(masterproblem_data, subproblem_data, instance)
end

function originalModel(instance::ModelInstance)
    model = Model(CPLEX.Optimizer)

    @variables(model, begin
        z >= 0
        0<= w <= 1 
        ϕ[i=1:instance.arcs] >= 0
    end)

    even_demands = findall(x->(x%2==0), instance.demands)
    odd_demands = setdiff(1:instance.vertices, even_demands)

    # # Cria as restrições para as demandas pares e ímpares utilizando as definições do modelo
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