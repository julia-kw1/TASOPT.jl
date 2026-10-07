module OptimizationCommon

using Printf

struct Constraint
    name::String
    value::Float64
    limit::Float64
    penalty::Float64
end

function clamp_initial!(initial, lower, upper)
    for index in eachindex(initial)
        initial[index] = clamp(initial[index], lower[index], upper[index])
    end
    return initial
end

function print_setup(optimizer, initial, lower, upper, tolerance, maxeval)
    println("="^60)
    println("TASOPT optimization setup")
    println("="^60)
    println("Optimizer: $(optimizer.algorithm)")
    println("Relative tolerance: $tolerance")
    println("Maximum evaluations: $maxeval")
    for index in eachindex(initial)
        println("  x[$index] = $(round(initial[index]; digits = 3)) ∈ [$(lower[index]), $(upper[index])]")
    end
    println("="^60)
end

function print_progress(iteration::Int, diagnostics, violations)
    if iteration == 1 || iteration % 10 == 0
        @printf("%-5s", "Iter")
        for (name, _) in diagnostics
            @printf("│ %-10s", name)
        end
        println()
    end

    @printf("%-5d", iteration)
    for (_, value) in diagnostics
        @printf("│ %10.3f", value)
    end
    for violation in violations
        @printf("│ %10s", violation.name * "!")
    end
    println()
end

end
