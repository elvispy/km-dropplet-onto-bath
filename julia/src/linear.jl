function solve_linear_system(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real})
    rcond_val = 0.0
    try
        rcond_val = 1 / LinearAlgebra.cond(A)
    catch
        rcond_val = 0.0
    end
    NEAR_SINGULAR[] = NEAR_SINGULAR[] || (rcond_val < NEAR_SINGULAR_TOL)
    return A \ b
end

