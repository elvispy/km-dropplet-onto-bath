# Checkpoint: Refactoring `solve_motion.m`

## Current Status
- **Objective**: Refactor the monolithic `solve_motion.m` into a modular, testable structure while maintaining bit-for-bit (or within tolerance) output compatibility with `solve_motion_original.m`.
- **Progress**: 
    - Created `test_solve_motion.m` for automated regression testing.
    - Simplified `solve_motion.m` main loop and extracted `attempt_step`, `update_state`, and `save_result_safe`.
    - Fixed a bug in `solveDD0` call where arguments were misaligned.
    - **Crucial Discovery**: The original code depends on history-dependent contact state tracking (contact angle hysteresis). Implemented `find_nlprev` in `solve_motion.m` to replicate this behavior.
- **Last Action**: Started a long-running test via `test_solve_motion.m` (using `U0=10, N=15`) to verify if the hysteresis fix prevents the divergence previously seen at `t=0.105`.

## Next Steps for the Incoming Agent
1. **Check Test Result**: Read the output of the background command/terminal from the last turn. If successful, both `original` and `new` versions should finish without crashing.
2. **Differential Debugging**: If divergence persists, use the `debug_solve_motion.log` (generated when `debug_flag=true`) to compare `try_k` results between the new version and any instrumentation added to the original.
3. **Ma Hardcoding**: `solvenDDCusp.m` has `Ma = 4*pi/3` hardcoded, ignoring the input parameter. This was temporarily changed but reverted to ensure the "new" version behaves *exactly* like the "original" first. Consider fixing this once parity is achieved.
4. **Final Validation**: Ensure all `.mat` outputs (`z`, `numl`, `etas`, etc.) match between versions.

## Key Files
- `matlab/1_code/simulation_code/solve_motion.m`: The refactored version.
- `matlab/1_code/simulation_code/solve_motion_original.m`: The legacy version.
- `matlab/test_solve_motion.m`: The regression test script.
- `.context/LOGIC_FINDINGS.md`: Detailed notes on the physics/logic divergence found.
