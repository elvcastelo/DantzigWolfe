using CPLEX
using JuMP

"""
    ModelInstance

Instância de um problema do modelo onde os campos representam as suas propriedades úteis.
"""
struct ModelInstance
    arcs_indexes::Matrix{Int}
    adjacency_matrix::Matrix{Int}
    capacity_matrix::Matrix{Int}
    arcs_set::Vector{Tuple}
    demands::Vector{Int}
    vertices::Int
    arcs::Int
end

"""
    model_from_file(path::String)

Lê os arquivos no ficheiro `path` e os transforma em uma matriz de adjacência, uma matriz de capacidade e um vetor de demandas.

# Argumentos
- `path`: Caminho do ficheiro.
"""
function model_from_file(path::String, verbose=false)
    files = readdir(path, join=false)
    for file in files
        if file == "capacities"
            continue
        end

        if verbose; printstyled("[model_from_file] Trabalhando com o arquivo \"$file\"\n", color=:blue, bold=true) end
        # Realiza uma concatenação do caminho absoluto do arquivo
        join_file = path * "\\" * file
        join_capacities_file = path * "\\capacities\\" * file

        # Lê os arquivos por linhas, retornando um vetor de strings
        file_buffer = readlines(join_file, keep=false)
        file_capacities_buffer = readlines(join_capacities_file, keep=false)

        # Número de vetores e metade do número de arcos respectivamente.
        n, m = split(file_buffer[1])
        # Conversão dos valores para inteiro
        n, m = parse(Int, n), parse(Int, m)

        # Cria um vetor de demanda nos vértices
        demands = zeros(Int, n)
        if verbose; printstyled("[model_from_file] Obtendo vetor de demandas \n", color=:blue, bold=true) end
        # Visita cada uma das próximas n linhas para obter a demanda de cada vértice
        for line = 2:n+1
            if file_buffer[line] != "0"
                demand = parse(Int, file_buffer[line])
                # A linha é substraída por um pois embora inicie a contagem na segunda linha em diante essas se referem ao primeiro vetor em diante.
                demands[line-1] = demand
            end
        end

        # Matriz de adjacência inicial
        adjacency_matrix = zeros(Int, n, n)
        # Matriz capacidade inicial
        capacity_matrix = zeros(Int, n, n)
        # Vetor de arcos
        arc_set = Vector{Tuple}()
        # Matriz dos índices dos arcos
        arc_indexes = zeros(Int, n, n)

        # Valor de soma para normalização dos índices
        i = 1

        if verbose; printstyled("[model_from_file] Criando estrutura da instância \n", color=:blue, bold=true) end
        # Visita cada uma das próximas m linhas para obter os arcos e criar a matriz de adjacência
        @inbounds for line = n+2:length(file_buffer)
            # Obtenção dos vértices que formam o arco
            u, v = split(file_buffer[line])
            # Obtenção das capacidades dos arcos
            capacityIn, capacityOut = split(file_capacities_buffer[line-n-1])
            # Conversão das capacidades para inteiro
            capacityIn, capacityOut = parse(Int, capacityIn), parse(Int, capacityOut)
            # Conversão dos valores para inteiro e incrementa em 1 uma vez que o índice começa em 1 ao invés de 0
            u, v = parse(Int, u)+1, parse(Int, v)+1

            push!(arc_set, (u, v))
            push!(arc_set, (v, u))

            adjacency_matrix[u, v] = 1
            adjacency_matrix[v, u] = 1

            capacity_matrix[u, v] = capacityIn
            capacity_matrix[v, u] = capacityOut

            # Gambiarra para normalizar os índices
            if line - n - 1 == 1
                arc_indexes[u, v] = line - n - 1
                arc_indexes[v, u] = line - n
            else
                arc_indexes[u, v] = line - n - 1 + i
                arc_indexes[v, u] = line - n + i
                i += 1
            end
        end

        if verbose; printstyled("[model_from_file] Inicializando modelo.\n", color=:blue, bold=true) end
        model_instance = ModelInstance(arc_indexes, adjacency_matrix, capacity_matrix, arc_set, demands, n, 2m)

        # Inicializa a decomposição de Dantzig-Wolfe
        initialize(model_instance, verbose)
    end
end