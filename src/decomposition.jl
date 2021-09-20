"""
    DantzigWolfe(master_data::MasterProblemData, sub::SubProblem, sub_data::SubProblemData, verbose::Bool)::MasterProblem

Realiza a decomposição de Dantzig-Wolfe. O algoritmo cria o modelo mestre inicial, uma vez que o problema mestre é o dual, ao invés de adicionarmos variáveis estaremos adicionando restrições, o que possibilita uma maior facilitadade para efetuar essas manipulações e na possibilidade de melhorias usando Lazy Constraints. O subproblema é resolvido a cada iteração para se obter um ponto extremo, de maneira que a sua função objetiva é atualizada a cada iteração.

Dado a natureza do problema fornecido, o algoritmo não obtêm as direções extremas pelo fato de o poliedro do subproblema ser limitado, dado as suas restrições de capacidade e não negatividade.

# Argumentos
- `master_data`: Dados do problema mestre que serão utilizados para a construção do modelo e para a obtenção de `A'` (`A_line`).
- `sub`: Subproblema de fluxo criado anteriormente para obter o primeiro ponto inicial viável para o problema mestre.
- `sub_data`: Dados do subproblema.
- `verbose`: Booleano indicando se deseja receber saídas conforme o algoritmo executa.
"""
function DantzigWolfe(master_data::MasterProblemData, sub::SubProblem, sub_data::SubProblemData, verbose::Bool)::MasterProblem
    if verbose; printstyled("[DantzigWolfe] Iniciando decomposição de Dantzig-Wolfe.\n", color=:blue, bold=true) end

    iter = 1
    same_objective_iter = 0

    z = 10^6
    # bounds = ""
    master = build_dualmaster(master_data)
    A_line = master_data.A[:,3:end]

    while true
        optimize!(master.model)

        if iter == 1
            z = objective_value(master.model)        
        else
            if z - objective_value(master.model) == 0
                if same_objective_iter == 10
                    if verbose; printstyled("[DantzigWolfe] O problema não obteve nenhuma melhora.", color=:red, bold=true) end
                    # println(bounds)
                    return master
                else
                    same_objective_iter += 1
                end
            else
                z = objective_value(master.model)
                same_objective_iter = 0
            end
        end

        
        if termination_status(master.model) == MOI.OPTIMAL
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema mestre possui uma solução ótima.\n", color=:blue, bold=true) end
            
            @objective(sub.model, Min, (sub_data.c' - value.(master.λ)' * A_line) * sub.x - value(master.λ_0))
            optimize!(sub.model)

            if termination_status(sub.model) == MOI.OPTIMAL
                if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui uma solução ótima.\n", color=:blue, bold=true) end
                δ = objective_value(sub.model)
                # l_bound = z + δ
                # u_bound = z
                # bounds *= "Iteração $iter: $l_bound <= V[PM] <= $u_bound\n"

                if δ >= 0
                    if verbose; printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ >= 0, portanto temos uma solução ótima.\n", color=:blue, bold=true) end
                    # println(bounds)
                    return master
                else
                    if verbose printstyled("[DantzigWolfe] Iteração $iter: O subproblema possui δ < 0, portanto iremos acrescentar uma nova coluna. \n", color=:blue, bold=true) end
                    # Obtêm o valor do ponto extremo
                    extreme_point = value.(sub.x)
                    # Adiciona uma nova restrição ao modelo
                    @constraint(master.model, master.λ' * (A_line * extreme_point) + master.λ_0 <= 0)
                end
            end
        else
            if verbose; printstyled("[DantzigWolfe] Iteração $iter: O problema é ilimitado ou inviável. Retornando.\n", color=:blue, bold=true) end
            # println(bounds)
            return master
        end

        iter += 1

    end
end