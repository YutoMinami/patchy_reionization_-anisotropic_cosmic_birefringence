# HANDOFF.md

## Project
Investigate whether the patchy-reionization-induced contribution to anisotropic cosmic birefringence (CB) can compete with or dominate the usual ALP fluctuation contribution.

Core question:
```math
C_L^{\alpha\alpha}
=
A_\phi^2 C_L^{\phi\phi}
+
A_\tau^2 C_L^{\tau\tau}
+
2A_\phi A_\tau C_L^{\phi\tau},
```
and in particular whether
```math
A_\tau^2 C_L^{\tau\tau}
```
or the cross term can be comparable to or larger than
```math
A_\phi^2 C_L^{\phi\phi}.
```

---

## Main formulas

Define
```math
R_\tau(L)
\equiv
\frac{A_\tau^2 C_L^{\tau\tau}}{A_\phi^2 C_L^{\phi\phi}},
\qquad
R_\times(L)
\equiv
\frac{2A_\phi A_\tau C_L^{\phi\tau}}{A_\phi^2 C_L^{\phi\phi}}.
```

Using
```math
A_\phi = -\frac{g_{a\gamma}}{2},
\qquad
A_\tau = -\frac{g_{a\gamma}}{2}\,\dot\phi\,\frac{d\eta}{d\tau},
```
we get
```math
\frac{A_\tau}{A_\phi}
=
\dot\phi\,\frac{d\eta}{d\tau}
\equiv \mathcal A.
```

Therefore
```math
R_\tau(L)
=
\mathcal A^2 \frac{C_L^{\tau\tau}}{C_L^{\phi\phi}},
```
and
```math
R_\times(L)
=
2\mathcal A\,\frac{C_L^{\phi\tau}}{C_L^{\phi\phi}}.
```

A useful alternative parameterization for the cross term is
```math
\rho_L \equiv \frac{C_L^{\phi\tau}}{\sqrt{C_L^{\phi\phi}C_L^{\tau\tau}}},
```
so that
```math
R_\times(L)
=
2\mathcal A\,\rho_L\,\sqrt{\frac{C_L^{\tau\tau}}{C_L^{\phi\phi}}}.
```

---

## Interpretation

- `A_phi^2 C_L^{phi phi}`: genuine ALP-induced anisotropic CB.
- `A_tau^2 C_L^{tau tau}`: effective birefringence-like term from patchy reionization, through time-dependent ALP background.
- `cross term`: only relevant if ALP fluctuations and optical-depth fluctuations are correlated.

The project focus shifted from “derive delta alpha_tau” to:
- quantify whether the patchy term can dominate,
- determine where in parameter space this can happen,
- or show a no-go result if it never happens after physical normalization.

---

## Effective patchy term derivation summary

The useful first-order relation is
```math
\delta\alpha_\tau(\hat n)
\simeq
-\frac{g_{a\gamma}}{2}
\dot\phi(\bar\eta_{\rm rei})
\frac{d\eta}{d\tau}
\delta\tau(\hat n).
```

Physical interpretation:
patchy reionization does not directly rotate polarization; it shifts the effective emission/scattering time, and a time-dependent ALP background turns that shift into an effective CB fluctuation.

This derivation is not itself the paper’s main result. It is background/formalism only.

---

## Numerical strategy used so far

We set up a toy numerical pipeline to compute
```math
\mathcal A(m_a) = \dot\phi_{\rm conf}(\eta_{\rm rei}) \frac{d\eta}{d\tau},
```
then evaluate toy
```math
R_\tau(L)=\mathcal A(m_a)^2 \frac{C_L^{\tau\tau}}{C_L^{\phi\phi}}.
```

### Background ALP equation
Using physical time:
```math
\ddot\phi + 3H\dot\phi + m_a^2 \phi = 0.
```

Using `x = ln a`:
```math
\phi'' + \left(3 + \frac{d\ln H}{dx}\right)\phi' + \frac{m_a^2}{H^2}\phi = 0,
```
where prime is derivative w.r.t. `x`.

Then
```math
\dot\phi_{\rm phys}=H\,\frac{d\phi}{dx},
\qquad
\dot\phi_{\rm conf}=a\,\dot\phi_{\rm phys}=aH\frac{d\phi}{dx}.
```

For optical depth:
```math
\frac{d\eta}{d\tau}=\frac{1}{a n_e \sigma_T}.
```

### Toy spectra
Used toy models
```math
C_L^{\tau\tau} = C_{\tau0}\exp\left[-\frac{(L-L_c)^2}{2\sigma_L^2}\right]
```
and
```math
C_L^{\phi\phi}=C_{\phi0}(L/L_*)^{-n_\phi}
```
(optionally with cutoff).

This was only to test whether patchy dominance is even plausible in principle.

---

## Numerical results obtained so far

A heavy scan was run. One example reported:

- Picked mass:
  `m_pick = 3.919406774847229e-28 eV`
- At that point:
  `A_dimless = 7.653902236886374e9`
- Toy result:
  `max R_tau = 7.052221316946602e19`
- Maximum near:
  `L = 356`

Also observed:
- around `L ~ 356`, the toy `C_L^{tau tau}` exceeded toy `C_L^{phi phi}`.

This is important qualitatively but **must not yet be interpreted as a physical claim**.

---

## Important caveat: current absolute normalization is not physical yet

The code used `phi_ini = 1.0` with no physical normalization.
Therefore the current `A_dimless` and `R_tau` absolute values are not yet physically meaningful.

What **is** meaningful at this stage:
- the shape of `A(m_a)` as a function of mass,
- the fact that there is a preferred mass region,
- the fact that patchy dominance can occur in toy setups near the bubble-scale multipole window.

What is **not yet** safe to claim:
- physically realistic `R_tau ~ 10^19`,
- patchy dominance in real parameter space,
- reinterpretation of observational limits.

---

## Scaling test already performed

A scaling test was run by varying `phi_ini`:

- `phi_ini=1.0e+00, A=3.815e+09, Rmax=1.752e+19`
- `phi_ini=1.0e-03, A=3.813e+06, Rmax=1.751e+13`
- `phi_ini=1.0e-06, A=3.330e+03, Rmax=1.335e+07`
- `phi_ini=1.0e-09, A=-1.972e+02, Rmax=4.679e+04`

Interpretation:
- first two points are close to expected linear/quadratic scaling,
- last two points break scaling and even flip sign,
- this indicates numerical instability / interpolation error / solver tolerance issues,
- not physical behavior.

Since the ODE is linear in `phi`, the robust approach is:

1. solve once with `phi_ini = 1`,
2. define `A_unit(m)` from that solution,
3. rescale later by physical amplitude:
```math
   A(m;\phi_{\rm amp}) = \phi_{\rm amp} A_{\rm unit}(m).
```

Do **not** repeatedly solve for tiny `phi_ini`.

---

## Current assessment

The present situation is actually promising:

- In toy spectra, patchy term can easily dominate near specific `L`.
- Therefore the real question is no longer whether the effect exists formally.
- The real question is whether it survives **physical normalization** of the ALP background.

This is exactly where the project should go next.

---

## Updated status after the latest numerical round

The project has now moved beyond the original notebook-only feasibility stage.

### What is now in place

- `scripts/02-check_linearity.py`
  confirms that the old small-`phi_ini` pathology was numerical, not physical
- `scripts/03-check_convergence.py`
  shows that `DOP853` is the safer default, especially toward the heavier/oscillatory mass range
- `scripts/04a-scan_a_unit.py`
  computes `A_unit(m_a)` and writes progress incrementally to `results/04a-a-unit/A_unit_scan.csv`
- `scripts/04b-scan_rtau_from_aunit.py`
  reads a saved `A_unit` scan and computes `R_{\tau,\max}^{unit}(m_a)` and `phi_needed(m_a)`
- `scripts/05-visualize_phi_needed.ipynb`
  provides human-readable 2D plots for `A_unit`, `R_{\tau,\max}^{unit}`, and `phi_needed`

### Current best mass from the toy analysis

Using the completed `A_unit` scan out to `m_a = 10^{-26} eV`, the current best mass on the present coarse grid is

```math
m_{\rm best} = 5.878016 \times 10^{-27}\ {\rm eV}.
```

Representative toy results at that mass include:

- `nphi_2_no_cut`
  - `Rtau_max_unit = 2.337756e+20`
  - `L_peak = 356`
- `nphi_3_Lcut_800`
  - `Rtau_max_unit = 3.609159e+20`
  - `L_peak = 389`

Corresponding required amplitudes are already down to the $10^{-11}$ to $10^{-10}$ range in the toy analysis.

The project has now also gone through a matched rerun:

- `scripts/11-recompute_matched_scan.py`
  recomputes $A_{\rm unit}$ and $\phi_{\rm amp,max}$ from the same high-precision background solution
- `scripts/12-ratio_with_matched_scan.py`
  recomputes $\phi_{\rm needed} / \phi_{\rm amp,max}$ from the matched scan, without overwriting the earlier `08` sanity-check products

Current matched output:

- `results/11-matched-scan-global-split/matched_scan.csv`
- `results/12-ratio-with-matched-scan/matched_ratio_scan.csv`

The important point is that the matched rerun did **not** qualitatively change the conclusion:

- $\phi_{\rm needed} / \phi_{\rm amp,max}$ remains far below 1 over the scanned mass range
- for $R_{\rm target} = 1$, the matched-rerun ratio is still at the level of roughly $10^{-15}$ to $10^{-11}$

So the current state of the project is:

- toy patchy dominance is easy to realize in the spectrum-level model,
- the required amplitude is tiny compared with the current phenomenological amplitude bound,
- therefore the bottleneck is no longer basic feasibility, but rather how to formulate the paper-quality claim and its caveats.

---

## High-priority next tasks for Codex

### Task 1: Treat `11/12` as the current best numerical result
Do not overwrite older numbered products such as `08`, `09`, `10`.
Those are useful as historical sanity checks.

The current numerically preferred chain is:

1. `04a` / `04b` for toy $\phi_{\rm needed}$
2. `11` for matched high-precision $A_{\rm unit}$ and $\phi_{\rm amp,max}$
3. `12` for the matched ratio $\phi_{\rm needed} / \phi_{\rm amp,max}$

### Task 2: Turn the current result into a paper claim
The numerics now strongly suggest:

- patchy dominance is toy-feasible,
- the required amplitude is far below the current phenomenological bound,
- and this survives the matched rerun.

So the next job is to phrase the result carefully:

- what is already robust,
- what is still toy/model-dependent,
- what caveats remain about normalization and $C_L^{\tau\tau}$ modeling.

### Task 3: Revisit observational reinterpretation with the matched result in hand
The natural follow-up is to reconnect this to:

- anisotropic CB limits such as $A_{\rm CB}$,
- possible reinterpretation as a bound on $A_\tau^2 C_L^{\tau\tau}$,
- benchmark scenarios tied to isotropic CB.

This should now be done using the matched-rerun understanding, not the earlier mixed `04a/07a` combination.

### Task 4: Keep cross term secondary for now
Do not prioritize $C_L^{\phi\tau}$ modeling yet.
For now:

- set $\rho_L = 0$,
- maybe provide an upper-envelope estimate using $|\rho_L| \le 1$,
- leave realistic cross modeling for later.

### Task 5: Keep the finite-width reionization issue as a TODO, not a blocker
The current treatment uses the effective emit-time-shift approximation for `\delta\alpha_\tau`.
This is good enough for the present feasibility and amplitude-budget stage.

If the patchy term remains interesting after the current matched analysis, then the next physics-level refinement is to test the validity of this approximation against a finite-width reionization kernel treatment.

---

## Suggested near-term paper framing

Possible title:
**Can patchy reionization dominate anisotropic cosmic birefringence in axion-photon coupling models?**

Suggested structure:
1. Introduction
2. Formalism
3. Scaling argument for patchy vs genuine term
4. ALP background dynamics and $A_{\rm unit}(m_a)$
5. Required amplitude for patchy dominance
6. Physical normalization from ALP energy-density bounds
7. Results / viability or no-go
8. Discussion

Most likely central figure:

- matched $\phi_{\rm needed}(m_a)$ for $R_{\tau,\max} = 1$
vs
- matched $\phi_{\rm amp,max}(m_a)$ from the same background solutions

with the old non-matched version retained only as a sanity-check appendix or internal log.

---

## Numerical implementation notes

### Current solver approach
- `solve_ivp`
- variable: `x = ln a`
- background cosmology: flat LCDM with radiation
- `z_rei ~ 7.7`
- `deta/dtau = 1/(a n_e sigma_T)`

### Potential numerical issues
- interpolation at `x_rei` may be too crude if oscillations are fast,
- `RK45` may struggle if oscillations become rapid,
- absolute tolerance dominates for tiny `phi_ini`,
- therefore unit-amplitude response is the safer path.

### Code-level recommendation
If the oscillatory regime becomes important, consider:
- denser sampling near reionization,
- event-free fixed dense output evaluation,
- maybe stiff solver or oscillation-aware treatment if needed.

But first, keep it simple and stabilize the unit-amplitude scan.

---

## Short summary for Codex

The project has moved past the initial feasibility phase.

What is already in hand:

- a stable unit-response workflow,
- toy $\phi_{\rm needed}(m_a)$,
- a phenomenological $\phi_{\rm amp,max}(m_a)$,
- and a matched rerun showing that the very small ratio $\phi_{\rm needed} / \phi_{\rm amp,max}$ is not an artifact of mixing two different solver setups.

So the next job is no longer “make the scan work”.
The next job is to articulate the result:

- what exactly is robust already,
- what is still phenomenological,
- and how this feeds into observational reinterpretation and a paper narrative.
