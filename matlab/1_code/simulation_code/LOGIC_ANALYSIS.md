### Logic Analysis of Original `solve_motion`

The original `solve_motion_original.m` uses a nested `if/elseif` structure based on `numl` (number of contact points) from the *previous* accepted step (`tentative_index`). Here is the breakdown of the decision tree for `numlTent` (the candidate for current step).

#### 1. Case `numl == 0` (In flight)
- **Primary Attempt**: Solve with `numlTent = 0` (`solveDD0`).
- **Check**: Is `abs(errortan) < 0.5`?
    - **Yes**: Accept `numlTent = 0`.
    - **No**: 
        - Try `numlTent = 1` (`solvenDDCusp` with `nl=1`). Get `err1`.
        - Try `numlTent = 2` (`solvenDDCusp` with `nl=2`). Get `err2`.
        - **Selection**:
            - If `abs(err1) < abs(err2)`: Accept `numlTent = 1`.
            - Else: Reduce time step (`reduc = 1`). (Wait, it doesn't accept 2? Line 430 accepts 1. Line 431 `else` -> `reduc`. It seems it only accepts 1 or fails to reduc. It *does not* accept 2 directly from 0? Careful reading needed.)
            - *Correction*: Line 423 `if abs(err1) < abs(err2)`. If true, accept 1. If false (err2 better or equal), it reduces dt. This implies a transition 0 -> 2 is forbidden or considered unstable, enforcing continuity (0->1).

#### 2. Case `numl == 1` (Point contact)
- **Primary Attempt**: Solve with `numlTent = 0`.
- **Check**: Is `abs(err0) < 0.5`?
    - **Yes**: Accept `numlTent = 0`. (Release)
    - **No**:
        - Try `numlTent = 1`. Get `err1`.
        - Try `numlTent = 2`. Get `err2`.
        - **Selection**:
            - If `abs(err1) < abs(err2)`: Accept `numlTent = 1` (Maintain).
            - Else:
                - Try `numlTent = 3`. Get `err3`.
                - If `abs(err2) < abs(err3)`: Accept `numlTent = 2` (Spread).
                - Else: Reduce time step. (Again, enforcing 1 -> 2, blocking 1 -> 3).

#### 3. Case `numl == 2` (Two points)
- **Primary Attempt**: `solveDD0` (Test flight/release? No, just check error?). Actually Line 482 calls `solvenDDCusp`? No, calls `solveDD0` -> `errortan(1)`.
- **Check**: Is `abs(err0) < 0.5`?
    - **Yes**: Reduce time step. (Wait, if `numl=2` and `err0` (flight) is good, it reduces dt? Why? Maybe rapid detachment requires finer resolution? Or it considers jump 2->0 invalid and forces refinement to capture 2->1->0?)
    - **No**: (Flight is bad, continue contact logic)
        - Try `numlTent = 2`. Get `err2`.
        - Try `numlTent = 1`. Get `err1`.
        - **Selection**:
            - If `abs(err1) < abs(err2)`: Accept `numlTent = 1` (Retract).
            - Else (2 is better than 1):
                - Try `numlTent = 3`. Get `err3`.
                - If `abs(err2) < abs(err3)`: Accept `numlTent = 2` (Maintain).
                - Else (3 is better than 2):
                    - Try `numlTent = 4`. Get `err4`.
                    - If `abs(err3) < abs(err4)`: Accept `numlTent = 3` (Spread).
                    - Else: Reduce dt.

#### 4. Case `numl > 2` (General)
- General pattern seems to be:
    - Compare `k` vs `k-1`.
    - If `k-1` better:
        - Check `k-2`. 
        - If `k-1` better than `k-2`: Accept `k-1`.
        - Else: Reduce dt.
    - Else (`k` better than `k-1`):
        - Check `k+1`.
        - If `k` better than `k+1`: Accept `k`.
        - Else (`k+1` better):
            - Check `k+2`.
            - If `k+1` better than `k+2`: Accept `k+1`.
            - Else: Reduce dt.

### Core Philosophy Inference
1.  **Laziness**: The code prioritizes the *current* state (or 0) if valid.
2.  **Continuity**: It strictly enforces `|k_new - k_old| <= 1`. It rarely jumps more than 1 (except 0->0). Transitions like 0->2 or 1->3 trigger a time-step reduction (`reduc`), effectively forcing the simulation to take smaller steps to resolve the intermediate state (0->1->2).
3.  **Local Optimization**: It looks at neighbors `k-1` and `k+1`. If a neighbor is better, it considers moving. But it also checks `k+2` or `k-2` to validate the move (making sure it's a local minimum, not a slope towards a jump).

### Flaw in My Refactor (`solve_motion.m`)
- **Greedy Global Search**: My code checked `k, k-1, k+1` (and `k-2, k+2`) and simply picked the *best* error.
- **Violation of Continuity**: If `k=0` and `err(2) < err(0)` and `err(2) < err(1)`, my code might have jumped to 2 (if I enabled that search). Even if restricted to `k+/-1`, my code didn't enforce the "Reduce dt if best is `k+2`" logic.
- **Missing "Stickiness"**: My code didn't implement "If `err(k)` is good, stop". It always checked neighbors. This allows "phantom" solutions with lower error (but unphysical) to hijack the state.

### Mitigation Plan
1.  **Refactor `attempt_step`** to implement the specific decision tree logic of the original.
2.  Instead of a generic `k` optimizer, implement a `select_next_contact_state(current_k, ...)` function that encapsulates the `if/elseif` logic.
3.  This function will return the `next_k` OR a signal to `reduce_dt`.

