function initialize_decomposition(masterproblem_data::MasterProblemData, subproblem_data::SubProblemData, instance::ModelInstance)

    master = build_dualmasterproblem(masterproblem_data)
    subproblem = build_subproblem(subproblem_data, [0 0], [0], 0, instance, false, true)

    iter = 0
    A_line = masterproblem_data.A[:,3:end]

    function dantzigwolfe_lazy_constraint_callback(cb_data)
        iter += 1

        λ_current = callback_value.(Ref(cb_data), master.λ)
        λ_0_current = callback_value(cb_data, master.λ_0)

        @objective(subproblem.model, Min, (subproblem_data.c' - λ_current' * A_line) * subproblem.x - λ_0_current)
        optimize!(subproblem.model)

        δ = objective_value(subproblem.model)

        println("Iteração $iter: $δ")
        if δ < 0
            printstyled("Adicionando coluna\n", color=:green, bold=true)
            extreme_point = value.(subproblem.x)
            new_constraint = @build_constraint(λ_current' * (A_line * extreme_point) + λ_0_current <= 0)

            MOI.submit(
                master.model,
                MOI.LazyConstraint(cb_data),
                new_constraint
            )
        end
    end

    MOI.set(master.model, MOI.NumberOfThreads(), 1)

    MOI.set(
        master.model,
        MOI.LazyConstraintCallback(),
        dantzigwolfe_lazy_constraint_callback
    )

    optimize!(master.model)
end