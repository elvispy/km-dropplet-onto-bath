# 20260124-1241-refactor-solve-motion Refactor solve_motion.m

## Owner + Lease
- owner_session: antigravity/2026-01-24T12:41
- lease_expires: 2026-01-24T13:41

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
- [ ] Locate and analyze `solve_motion.m` and `solve_motion_reference.m`.
- [ ] Create `implementation_plan.md`.
- [ ] Create reproduction/test script `test_solve_motion.m` using the CURRENT `solve_motion.m` to establish a baseline.
- [ ] Refactor `solve_motion.m`.
- [ ] Run comparison tests.

## Work Log
- 2026-01-24 12:41: Task initialized.

## Messages
(None)

## Handoff
- Task just started.
