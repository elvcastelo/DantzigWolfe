function DantzigWolfe(masterproblem_data::MasterProblemData, subproblem_data::SubProblemData, instance::ModelInstance, verbose::Bool)::MasterProblem
    if verbose; printstyled("[DantzigWolfe] Iniciando decomposição de Dantzig-Wolfe.\n", color=:blue, bold=true) end
    iter = 1
    z = 10^6
    bounds = ""
    while true
        masterproblem = build_masterproblem(masterproblem_data)
        optimize!(masterproblem.model)

        if iter == 1
            z = objective_value(masterproblem.model)
        else
            if z - objective_value(masterproblem.model) == 0
                if verbose; printstyled("[DantzigWolfe] O problema não obteve nenhuma melhora.", color=:red, bold=true) end
                return masterproblem
            else
                z = objective_value(masterproblem.model)
            end
        end

        
        if termination_status(masterproblem.model) == MOI.OPTIMAL
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema mestre possui uma solução ótima.\n", color=:blue, bold=true) end
            
            A_line = masterproblem_data.A[:,3:end]
            λ = masterproblem.λ
            λ_0 = masterproblem.λ_0

            subproblem = build_subproblem(subproblem_data, A_line, dual.(λ), dual(λ_0), instance)
            optimize!(subproblem.model)

            if termination_status(subproblem.model) == MOI.OPTIMAL
                if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui uma solução ótima.\n", color=:blue, bold=true) end
                δ = objective_value(subproblem.model)
                l_bound = z + δ
                u_bound = z
                bounds *= "Iteração $iter: $l_bound <= V[PM] <= $u_bound\n"

                if δ >= 0
                    if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ >= 0, portanto temos uma solução ótima.\n", color=:blue, bold=true) end
                    return masterproblem
                else
                    if verbose printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ < 0, portanto iremos acrescentar uma nova coluna. \n", color=:blue, bold=true) end
                    push!(masterproblem_data.V, value.(subproblem.x)) 
                end
            end
        else
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema é ilimitado ou inviável. Retornando.\n", color=:blue, bold=true) end
            return masterproblem
        end

        iter += 1

    end
end