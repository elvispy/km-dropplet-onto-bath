# Simulation Logic Analysis

## Divergence at `t=0.105`
The refactored `solve_motion.m` initially diverged from the original `solve_motion_original.m` around `t=0.105`. The original code successfully navigated a transition in contact geometry (likely advancing or receding contact line), while the refactored code entered an infinite loop of time-step reductions.

## Root Cause: Contact Angle Hysteresis Physics
The critical missing piece was the calculation of `nlprev` (previous number of contact points).

- **Original Code**: Explicitly searches the history `numl(tentative_index:-1:1)` to find the last state *different* from the current target state. This determines whether the contact line is advancing or receding, which changes the boundary condition (contact angle) used in `solvenDDCusp`.
- **Refactored Code (Initial)**: Simply passed `k` (current target) as `prev_k`, or used the immediately preceding steps without deep history search. This effectively treated every step as "static" or lost the directionality of the contact line movement, leading to incorrect pressure/deformation solutions that failed consistency checks (tangency error > 0.5).

## Fix Implementation
Implemented `find_nlprev` helper in `solve_motion.m` that replicates the `find` logic from the original code:
```matlab
% Search backwards from tentative_index for a value different from target_k
hist_vec = numl(tentative_index:-1:1);
co = find(hist_vec ~= target_k, 1);
np = numl(tentative_index - co + 1);
```
This ensures the solver uses the correct physical regime (advancing/receding) when calling the core `solvenDDCusp` function.

## Verification Status
Currently running `test_solve_motion.m` with the fix. Debug logs enabled to verify `nlprev` values and error convergence.
