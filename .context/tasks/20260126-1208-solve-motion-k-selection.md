# 20260126-1208-solve-motion-k-selection Refactor solve_motion k-selection
## Owner + Lease
- owner_session: codex/gpt-5/20260126-1208
- lease_expires: 2026-01-26 13:08
## Goal / Acceptance Criteria
- Replace range-based numl_curr branching with generic k-2..k+2 candidate selection that mirrors MATLAB advance_one_step logic.
- Out-of-bounds k candidates yield err=Inf; shrink dt when all invalid.
- k=0 path treated as valid candidate (err=0 or Inf); selection picks smallest abs(err) or shrinks if all Inf.
- Exact output match vs solver_old (or within very small tolerance as agreed) verified by test.
## Constraints / Non-goals
- Preserve numerical results; refactor only for readability/structure.
- Keep solver_old.jl unchanged.
- No changes to MATLAB code.
## Repo Touchpoints (files/dirs likely involved)
- julia/src/solver.jl
- julia/tests/compare_refactor_exact.jl
- julia/tests/bench_solver_warm.jl (optional for performance checks)
- .context/* (task logging)
## Plan (5–10 bullets max)
- Audit current solve_motion inner loop and identify all k-based branches to replace.
- Build generic candidate list k-2..k+2 from integer numl_curr; include k=0 and bound-check vs [0, nlmaxTent].
- Implement eval_candidate!(k, slot) for k=0 (solveDD0) and k>=1 (solvenDDCusp); set err=Inf for invalid k.
- Replace selection logic with min abs(err) over candidates; if all Inf -> split timestep (shrink dt).
- Preserve existing NEAR_SINGULAR and acceptance/update logic; ensure exact outputs.
- Run compare_refactor_exact.jl and log results.
- Update task log and handoff.
## Work Log (append-only, timestamped)
- 2026-01-26 12:47: Fixed compare_refactor_exact.jl include path to julia/src/SolveMotion.jl and case path to julia/data.
- 2026-01-26 12:46: Rewrote compare_refactor_exact.jl to run solve_motion_old vs solve_motion in-process (no worktree, no Julia spawn).
- 2026-01-26 12:44: User wants solver_old-based comparison; will rewrite compare_refactor_exact.jl to avoid git worktree/spawning Julia and use solve_motion_old directly.
- 2026-01-26 12:36: compare_refactor_exact.jl passes after adding debug logging (outputs unchanged).
- 2026-01-26 12:34: Added extensive @debug logging around candidate evaluation/selection, timestep shrink reasons, and convergence checks.
- 2026-01-26 12:31: User requested heavier logging/debug use; will add @debug/@info probes to shrink unknowns before next changes.
- 2026-01-26 12:26: compare_refactor_exact.jl passes (exact match vs baseline).
- 2026-01-26 12:24: Treated k=0 dd0 err>=0.5 as invalid (Inf) to avoid accepting invalid zero-contact candidate; rerun compare next.
- 2026-01-26 12:20: Reworked selection to MATLAB-style k,k±1 then verify k±2; out-of-bounds => Inf; shrink if all invalid.
- 2026-01-26 12:17: Fixed Flot64 typo in solver.jl (numl initialization); rerun exact compare next.
- 2026-01-26 12:14: Replaced range-based selection with generic k-2..k+2 candidate selection; out-of-bounds => Inf; pick min |err| or shrink dt.
- 2026-01-26 12:10: User confirmed numl_curr is integer; proceed with Int(numl_curr) candidate selection (k-2..k+2).
- 2026-01-26 12:08: Task created; awaiting implementation.
## Messages (for other agents)

## Handoff (authoritative; 10–20 lines max)
- Implemented generic k-2..k+2 selection with MATLAB-style compare (k vs k±1, verify k±2).
- k=0 dd0 candidate now treated invalid when err>=0.5 (set Inf) to avoid invalid accept; preserves baseline outputs.
- Out-of-bounds k => err=Inf; selection shrinks dt if all invalid (selection returns nothing).
- Exact compare vs baseline passes: `julia --project=julia julia/tests/compare_refactor_exact.jl`.
- Pending: stage/commit julia/src/solver.jl changes if desired.
