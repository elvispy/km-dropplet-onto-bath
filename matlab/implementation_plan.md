# Implementation Plan - Refactor solve_motion.m

## Analysis
The current `solve_motion.m` contains a `attempt_step_wrapper` function that encapsulates the contact point search logic. However, the implementation is verbose, with specific `if/elseif` blocks for `k=0`, `k=1`, `k=2`, and a generic block for `k>2`. This violates the DRY (Don't Repeat Yourself) principle and makes the code harder to maintain.

The reference file `matlab/.tmp/solve_motion_reference.m` uses a more generic approach in its main loop, testing neighbors `k-1` and `k+1`, and extending the search to `k-2` or `k+2` if the gradient suggests it.

## Refactoring Goal
Refactor `attempt_step_wrapper` in `solve_motion.m` to use a generic, data-driven approach for selecting candidate contact point counts (`k`). The new logic should:
1.  Remove special handling for `k=0, 1, 2`.
2.  Implement a generic search strategy that checks `k`, `k-1`, and `k+1`.
3.  Based on the error gradient, optionally check `k-2` or `k+2`.
4.  Handle boundary conditions (`k >= 0`) naturally.

## Proposed Logic (Pseudocode)
```matlab
function res = attempt_step_wrapper(prev_k, k, dt, st, PC, N, tolP)
    % 1. Try current k
    r_curr = try_k(k);
    if is_acceptable(r_curr); return r_curr; end

    % 2. Try neighbors k-1 and k+1
    r_prev = try_k(k-1); % Returns inf error if k-1 < 0
    r_next = try_k(k+1);

    % 3. Determine direction
    if r_prev.err < r_next.err
        best_neighbor = r_prev;
        direction = -1;
    else
        best_neighbor = r_next;
        direction = 1;
    end

    % 4. Compare with current (if we didn't return early)
    % Note: The reference logic compares neighbors against each other mostly.
    % We should follow the reference priority.

    % Reference logic adaptation:
    % Compare k, k-1, k+1.
    % If k-1 is best, check k-2.
    % If k+1 is best, check k+2.
    % If k is best, keep k.

    % ... implementation details ...
end
```

## Testing Plan
1.  **Regression Test Script**: Create `matlab/test_solve_motion.m`.
    - Setup a small domain (D5Q20) as requested.
    - Run `solve_motion` for ~10 timesteps.
    - Assert that execution finishes without error.
    - (Ideally) Compare output with a baseline, but without a running environment, we focus on code correctness and logical equivalence.

## Files to Modify
- `matlab/1_code/simulation_code/solve_motion.m`: Refactor `attempt_step_wrapper`.
- `matlab/test_solve_motion.m`: Create new file.
