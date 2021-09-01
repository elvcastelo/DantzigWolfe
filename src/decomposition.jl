function DantzigWolfe(masterproblem_data::MasterProblemData, subproblem_data::SubProblemData, instance::ModelInstance)
    printstyled("[DantzigWolfe] Iniciando decomposição de Dantzig-Wolfe.\n", color=:blue, bold=true)
    max_iter = 2000
    iter = 1
    while true
        masterproblem = build_masterproblem(masterproblem_data)
        optimize!(masterproblem.model)
        
        if termination_status(masterproblem.model) == MOI.OPTIMAL
            printstyled("[DantzigWolfe] Iteração $iter: O problema mestre possui uma solução ótima.\n", color=:blue, bold=true)
            
            A_line = masterproblem_data.A[:,3:end]
            λ = masterproblem.λ
            λ_0 = masterproblem.λ_0

            subproblem = build_subproblem(subproblem_data, A_line, dual.(λ), dual(λ_0), instance)
            optimize!(subproblem.model)

            if termination_status(subproblem.model) == MOI.OPTIMAL
                printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui uma solução ótima.\n", color=:blue, bold=true)

                if objective_value(subproblem.model) >= 0
                    printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ >= 0, portanto temos uma solução ótima.\n", color=:blue, bold=true)
                    return masterproblem
                else
                    printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ < 0, portanto iremos acrescentar uma nova coluna. \n", color=:blue, bold=true)
                    push!(masterproblem_data.V, value.(subproblem.x)) 
                end
            end
        else
            printstyled("[DantzigWolfe] Iteração $iter: O problema é ilimitado ou inviável. Retornando.\n", color=:blue, bold=true)
            return masterproblem
        end

        iter += 1
        # Evita que o computador pegue fogo
        if iter == max_iter
            printstyled("[DantzigWolfe] Iteração $iter: Quantidade máxima de iterações excedida. Encerrando algoritmo.\n", color=:red, bold=true)
            return masterproblem
        end

    end
end