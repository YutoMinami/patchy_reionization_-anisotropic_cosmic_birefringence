# AGENTS.md

## Purpose

This file gives working guidance for coding agents contributing to `PatchyReionization`.

The project studies whether patchy reionization can induce an effective anisotropic cosmic birefringence signal that competes with the standard ALP fluctuation contribution.

## Read order

Before making changes, read these in order:

1. `HANDOFF.md`
2. `RESEARCH_LOG.md`
3. `README.md`
4. `patchy_reionization.py`
5. `scripts/README.md`
6. any script or notebook cell you plan to modify

`HANDOFF.md` is the physics brief.
`RESEARCH_LOG.md` is the running project memory.
`patchy_reionization.py` is the shared implementation that notebook and scripts should follow.

## Current source-of-truth policy

- Keep the notebook because it is useful for plots and interactive inspection.
- Treat `patchy_reionization.py` plus the numbered scripts as the reproducible implementation.
- Avoid duplicating core numerical logic inside the notebook when a shared function can be reused.
- If notebook behavior changes, update the shared module first unless the change is purely presentation.
- Prefer the split workflow `04a` then `04b` over the older combined `04` script when runs are heavy or WSL2 stability matters.
- For heavy scans such as `04a` or `07a`, prefer splitting the mass range into smaller chunks and resuming from saved CSV outputs instead of running the full range in one shot.
- Once a numbered script or notebook has been used for a recorded result, prefer preserving it as-is for reproducibility.
- If a workflow needs a correction or a better-matched rerun, add a new continued number or a clearly versioned suffix instead of silently rewriting the old numbered artifact.

## Current priorities

Focus on the following unless the user redirects:

1. use a stable unit-amplitude response `A_unit(m_a)` with `phi_ini = 1`,
2. rescale linearly for external amplitudes instead of re-solving at tiny amplitudes,
3. keep dense-output evaluation at reionization for stable extraction of `A`,
4. add convergence diagnostics for representative masses,
5. compute `R_tau,max^unit(m_a)`,
6. compute `phi_needed(m_a)` for patchy dominance,
7. compare that requirement to physically allowed ALP amplitudes,
8. use `scripts/05-visualize_phi_needed.ipynb` when a human-friendly 2D view of the current results is helpful.

Do not prioritize detailed modeling of `C_L^{phi tau}` unless explicitly requested.

## Numerical guidance

- Prefer transparent and auditable numerical code.
- Preserve unit conventions in names or nearby comments.
- Separate toy inputs from physically normalized conclusions.
- If a result changes when `phi_ini` changes, assume a numerical issue first because the ODE is linear in `phi`.
- Check solver settings, tolerance scaling, interpolation method, and sampling near reionization before changing the physics interpretation.

## Notebook guidance

- The notebook should remain readable for a human researcher.
- Favor a few clear cells over many fragmented cells.
- Keep plots in the notebook, but keep core calculations in shared Python functions.
- If a long-running cell exists, label it clearly in markdown or comments.
- `scripts/05-visualize_phi_needed.ipynb` is currently the main visualization notebook for the saved `04a/04b` outputs.

## Validation expectations

Whenever you modify the pipeline, try to verify:

1. unit-response rescaling matches direct solves,
2. dense-output extraction is stable,
3. representative masses are converged across reasonable solver settings,
4. heavy runs write reusable intermediate outputs to disk when feasible,
5. any claim based on toy spectra is clearly labeled as toy-only.

For long or memory-hungry runs:

- write partial CSV outputs incrementally,
- add `--resume` support when practical,
- prefer mass-range split runs over single all-range runs on WSL2.

## Communication expectations

When summarizing results:

- separate formal derivation, toy-model behavior, and physically normalized conclusions,
- explicitly state caveats,
- do not overclaim from large toy `R_tau` values,
- emphasize the comparison `phi_needed` vs `phi_amp_max` as the likely main science result.
