function DantzigWolfe(masterproblem_data::MasterProblemData, subproblem_data::SubProblemData, instance::ModelInstance, verbose::Bool)::MasterProblemDual
    if verbose; printstyled("[DantzigWolfe] Iniciando decomposição de Dantzig-Wolfe.\n", color=:blue, bold=true) end

    iter = 1
    max_iter = 10000
    same_objective_iter = 0

    z = 10^6
    bounds = ""
    masterproblem = build_dualmasterproblem(masterproblem_data)
    while true
        optimize!(masterproblem.model)

        if iter == max_iter
            printstyled("[DantzigWolfe] Iteração $iter: O limite de iterações foi atingido.\n", color=:red, bold=true)
            return masterproblem
        end

        if iter == 1
            z = objective_value(masterproblem.model)        
        else
            if z - objective_value(masterproblem.model) == 0
                if same_objective_iter == 10
                    if verbose; printstyled("[DantzigWolfe] O problema não obteve nenhuma melhora.", color=:red, bold=true) end
                    println(bounds)
                    return masterproblem
                else
                    same_objective_iter += 1
                end
            else
                z = objective_value(masterproblem.model)
                same_objective_iter = 0
            end
        end

        
        if termination_status(masterproblem.model) == MOI.OPTIMAL
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema mestre possui uma solução ótima.\n", color=:blue, bold=true) end
            
            A_line = masterproblem_data.A[:,3:end]
            λ = masterproblem.λ
            λ_0 = masterproblem.λ_0

            subproblem = build_subproblem(subproblem_data, A_line, value.(λ), value(λ_0), instance)
            optimize!(subproblem.model)

            if termination_status(subproblem.model) == MOI.OPTIMAL
                if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui uma solução ótima.\n", color=:blue, bold=true) end
                δ = objective_value(subproblem.model)
                l_bound = z + δ
                u_bound = z
                bounds *= "Iteração $iter: $l_bound <= V[PM] <= $u_bound\n"

                if δ >= 0
                    if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ >= 0, portanto temos uma solução ótima.\n", color=:blue, bold=true) end
                    println(bounds)
                    return masterproblem
                else
                    if verbose printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ < 0, portanto iremos acrescentar uma nova coluna. \n", color=:blue, bold=true) end
                    # Obtêm o valor do
                    extreme_point = value.(subproblem.x)
                    # Adiciona uma nova restrição ao modelo
                    @constraint(masterproblem.model, λ' * (A_line * extreme_point) + λ_0 <= 0)
                end
            end
        else
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema é ilimitado ou inviável. Retornando.\n", color=:blue, bold=true) end
            println(bounds)
            return masterproblem
        end

        iter += 1

    end
end