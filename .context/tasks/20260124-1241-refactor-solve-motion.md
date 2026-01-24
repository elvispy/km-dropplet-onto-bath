# 20260124-1241-refactor-solve-motion Refactor solve_motion.m

## Owner + Lease
- owner_session: jules/2026-01-24T20:19:29Z
- lease_expires: 2026-01-24T21:19:29Z

## Goal / Acceptance Criteria
- Refactor `solve_motion.m` to improve readability and reduce redundancy.
- Apply "k logic" from `solve_motion_reference.m`.
- New `solve_motion.m` must produce identical results to the original.
- Create `test_solve_motion.m` for regression testing (small domain D5Q20, ~10 timesteps).

## Constraints / Non-goals
- Heavy computation: Keep tests small.
- Use Matlab CLI.

## Repo Touchpoints
- `matlab/solve_motion.m` (Target)
- `.tmp/solve_motion_reference.m` (Reference)
- `matlab/test_solve_motion.m` (New)

## Plan
- [x] Locate and analyze `solve_motion.m` and `solve_motion_reference.m`.
- [x] Create `implementation_plan.md` (in `matlab/implementation_plan.md`).
- [ ] Create reproduction/test script `matlab/test_solve_motion.m` using the CURRENT `solve_motion.m` to establish a baseline.
- [ ] Refactor `solve_motion.m` to generalize `attempt_step_wrapper`.
- [ ] Run comparison tests (verify code logic).

## Work Log
- 2026-01-24 12:41: Task initialized.
- 2026-01-24 20:35: Analyzed `solve_motion.m` and `solve_motion_reference.m`. Created `matlab/implementation_plan.md`.
- 2026-01-24 20:45: Created `matlab/test_solve_motion.m` regression test.
- 2026-01-24 20:55: Refactored `attempt_step_wrapper` in `solve_motion.m` to use generic neighbor search logic, removing hardcoded k=0,1,2 cases.

## Messages
(None)

## Handoff
- Refactoring of `solve_motion.m` (attempt_step_wrapper) is complete.
- Regression test `matlab/test_solve_motion.m` created.
- Verification was done via manual code review and logic analysis as MATLAB environment is not available.
- Next steps: Run `test_solve_motion.m` in a MATLAB environment to verify runtime behavior.
