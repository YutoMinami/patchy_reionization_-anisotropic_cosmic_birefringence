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

- $A_\phi^2 C_L^{\phi\phi}$: genuine ALP-induced anisotropic CB.
- $A_\tau^2 C_L^{\tau\tau}$: effective birefringence-like term from patchy reionization, through time-dependent ALP background.
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

After runs `21` and `22` (natural-unit normalization), the current picture is:

- With physical ALP-photon coupling $g_{a\gamma} = 1.4$–$4.0 \times 10^{-12}\ {\rm GeV}^{-1}$ (Chandra benchmark range), the physical response coefficient is
  ```math
  A_\tau^{\rm physical} \simeq 0.16\text{–}0.45,
  ```
  which puts the patchy term at roughly **1–10% of the anisotropic CB observational limit** (using the Dvorkin-Smith-inspired surrogate template).
- The patchy term is therefore **subdominant but non-negligible** for realistic couplings: it is not excluded, and it is not dominant.
- The earlier toy results (runs `04b`, `11`, `12`) showing $\phi_{\rm needed}/\phi_{\rm amp,max} \sim 10^{-15}$–$10^{-11}$ were computed in mixed code units without the $g_{a\gamma}$ factor. They established that patchy dominance is not forbidden by the amplitude budget in principle, but those dimensionless ratios are **not physical $A_\tau^{\rm physical}$ values** and should not be quoted as the final result.
- The paper-level claim is no longer "patchy can dominate" but rather "patchy induces a percent-to-ten-percent level contribution to the anisotropic CB power spectrum that existing constraints already partially probe, and that should be accounted for in future interpretations of CB data."

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

The matched rerun did **not** qualitatively change the conclusion in the mixed-unit framework:

- $\phi_{\rm needed} / \phi_{\rm amp,max}$ (in code units) remains far below 1 over the scanned mass range
- for $R_{\rm target} = 1$, the matched-rerun ratio was at the level of roughly $10^{-15}$ to $10^{-11}$ in those code units

**However, this ratio is not the same as a physical $A_\tau^{\rm physical}$.** Runs `21` and `22` subsequently identified a unit-system mismatch in `11`/`12`: the code-unit $A_{\rm unit}$ and $\phi_{\rm amp,max}$ each carry implicit unit factors that cancel differently from what a direct comparison implies. Once converted to natural units with explicit $g_{a\gamma}$, the physical response is $A_\tau^{\rm physical} \simeq 0.16$–$0.45$ (Chandra benchmark range), placing the patchy contribution at the percent-to-ten-percent level of the anisotropic CB observational limit.

The `11`/`12` runs remain useful as evidence that the amplitude budget is not violated and that solver-mixing was not an artifact. Their quantitative ratios should not be quoted as the final physical result; `21`/`22` are the correct reference for that.

---

## High-priority next tasks for Codex

### Task 1: Treat `21`/`22` as the current canonical result
Do not overwrite older numbered products such as `08`, `09`, `10`, `11`, `12`.
Those are useful as historical sanity checks and intermediate steps.

The current numerically preferred chain is:

1. `04a` / `04b` for toy $\phi_{\rm needed}$ (code-unit feasibility scan)
2. `11` for matched high-precision $A_{\rm unit}$ and $\phi_{\rm amp,max}$ (matched solver, still code units)
3. `12` for the matched ratio $\phi_{\rm needed} / \phi_{\rm amp,max}$ (code-unit amplitude budget check)
4. `21` for natural-unit conversion of the matched background quantities
5. `22` for the $D_L^{\alpha\alpha}$ budget in natural units with explicit $g_{a\gamma}$

Physical claims should be drawn from `21`/`22`, not from `11`/`12` alone.

### Task 2: Turn the current result into a paper claim
The numerics now indicate:

- with Chandra-benchmark coupling, the patchy term contributes at the **percent-to-ten-percent level** of the anisotropic CB observational limit,
- this is a non-negligible, potentially observable contribution,
- it is subdominant in the current benchmark but non-trivially so,
- and it is robust against solver-mixing artifacts (established by `11`/`12`).

The paper claim should be framed as:

- existing anisotropic CB constraints probe the total birefringence power and therefore already partially constrain the patchy term,
- the patchy contribution is not negligible at the percent level for Chandra-range couplings,
- a complete interpretation of CB data requires accounting for this term.

What is robust, what is still toy/model-dependent, and what caveats remain about $C_L^{\tau\tau}$ modeling should all be stated explicitly.

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
The current treatment uses the effective emit-time-shift approximation for $\delta\alpha_\tau$.
This is good enough for the present feasibility and amplitude-budget stage.

If the patchy term remains interesting after the current matched analysis, then the next physics-level refinement is to test the validity of this approximation against a finite-width reionization kernel treatment.

---

## Open caveats and known limitations

These are issues that must be stated explicitly in any paper-level claim. They are not blockers for the current stage, but they bound what can and cannot be concluded.

### Caveat 1: $C_L^{\tau\tau}$ is still a toy spectrum

The entire quantitative result rests on a log-normal surrogate for $C_L^{\tau\tau}$, calibrated loosely to Dvorkin-Smith (2009). The amplitude normalization of $C_L^{\tau\tau}$ has an uncertainty of at least a factor of several from:

- bubble size distribution ($\bar R$, $\sigma_{\ln R}$),
- reionization duration ($\Delta y$),
- radiative transfer corrections.

Any numerical claim about the size of the patchy contribution should be accompanied by a range over representative template parameters, not a single fiducial value.

### Caveat 2: Emit-time-shift approximation may miss a resonant amplification

The current effective patchy term

```math
\delta\alpha_\tau \simeq -\frac{g_{a\gamma}}{2}\,\dot\phi(\bar\eta_{\rm rei})\,\frac{d\eta}{d\tau}\,\delta\tau(\hat n)
```

is a thin-shell approximation that evaluates $\dot\phi$ at a single epoch. The full expression beyond this approximation involves integrating over the visibility function:

```math
\delta\alpha_\tau(\hat n) = -\frac{g_{a\gamma}}{2} \int d\eta\, g(\eta)\,\dot\phi(\eta)\,\frac{d\eta}{d\tau}\,\delta\tau(\hat n, \eta).
```

The key issue is not just that this approximation might be slightly wrong. The preferred mass $m_{\rm best} \sim 6\times10^{-27}$ eV satisfies $m_a \sim H(z_{\rm rei})$, meaning the ALP is just beginning to oscillate at reionization. At this mass, the oscillation period of $\dot\phi(\eta)$ is $\sim 1/m_a$, which corresponds to a comoving distance scale

```math
\lambda_{\rm osc} \sim \frac{2\pi\hbar c}{m_a c^2\, a_{\rm rei}} \simeq 0.060\ {\rm Mpc}\ (= 59.5\ {\rm kpc})
\quad\text{at }z_{\rm rei} \simeq 7.7.
```

**This is ~570 times smaller than the bubble scale $R_{\rm eff} \sim 34$ Mpc.** There is therefore no spatial resonance between ALP oscillations and the bubble distribution.

However, $N_{\rm osc} \equiv \Delta\chi_{\rm patchy}/\lambda_{\rm osc} \simeq 15{,}000$ oscillation periods fit within the patchy reionization epoch ($\Delta z \sim 3$, $\Delta\chi_{\rm patchy} \simeq 924$ Mpc). This means the thin-shell approximation picks a single oscillation phase of $\dot\phi$ while the full visibility-weighted integral averages over ~15,000 periods. The two can differ drastically depending on whether a partial phase cancellation occurs.

**The dominant concern is therefore temporal, not spatial**: the thin-shell $A_{\rm unit}$ may significantly overestimate (or in principle underestimate) $A_{\rm eff}$ because of phase averaging across the reionization epoch. The resonant-amplification scenario originally hypothesized does not apply at this mass.

*(Numerical result from `scripts/23-check_osc_scale.py` / `results/23-osc-scale/`.)*

#### Verification plan for Caveat 2

Two physical layers need to be distinguished.

**Layer 1 — temporal matching** (effect of finite reionization width on $A_{\rm eff}$):

The thin-shell approximation evaluates $\dot\phi$ at a single epoch $\bar\eta_{\rm rei}$. The corrected amplitude is

```math
A_{\rm eff}(m_a) = \int d\eta\; g(\eta)\,\dot\phi_{\rm conf}(\eta)\,\frac{d\eta}{d\tau}\bigg|_\eta,
```

where $g(\eta)$ is the visibility function. The ratio $A_{\rm eff}/A_{\rm unit}$ depends on how many oscillation periods of $\dot\phi$ fit inside the visibility window:

- $m_a \sigma_\eta \ll 1$: thin-shell is accurate
- $m_a \sigma_\eta \gg 1$: oscillations average out, suppression
- $m_a \sigma_\eta \sim 1$: phase-dependent, amplification or suppression possible

**Layer 2 — spatial matching** (resonant amplification from bubble-scale coherence):

The oscillation wavelength of $\dot\phi$ at $z_{\rm rei}$ is

```math
\lambda_{\rm osc} \sim \frac{2\pi c}{m_a (1+z_{\rm rei})}.
```

According to `23-check_osc_scale.py`, this spatial matching does **not** occur at the current preferred mass $m_{\rm best} \simeq 5.9\times10^{-27}$ eV. There, $\lambda_{\rm osc} \simeq 0.060$ Mpc, which is far smaller than the bubble scale $R_{\rm eff} \sim 34$ Mpc.

Spatial resonance would instead occur near the mass range

```math
m_{\rm res} \simeq 7\times10^{-30}\text{–}7\times10^{-29}\ {\rm eV},
```

where $\lambda_{\rm osc}$ becomes comparable to the representative bubble scale. If the bubble clustering power $P_{\delta\tau}(k)$ has significant support near $k \sim 2\pi/\lambda_{\rm osc}$ in this range, the line-of-sight window function $W(\chi) \propto g(\chi)\dot\phi(\chi)$ can oscillate coherently with the bubble distribution, and $C_L^{\alpha_\tau\alpha_\tau}$ may be enhanced relative to the thin-shell estimate.

**Verification steps** (in order of difficulty):

| Step | Content | Effort |
|------|---------|--------|
| 1 | Plot $\lambda_{\rm osc}(m_a)$ at $z_{\rm rei}$ alongside $R_{\rm eff}$ as a function of $m_a$ | ~half a day |
| 2 | Compute $A_{\rm eff}(m_a)$ by convolving dense-output $\dot\phi(\eta)$ with $g(\eta)$; plot $A_{\rm eff}/A_{\rm unit}$ | ~1–2 days |
| 3 | Compute the full Limber-integrated $C_L^{\alpha_\tau\alpha_\tau}$ using the oscillating $W(\chi)$ window; compare to thin-shell result | ~1 week |

Step 2 is a lightweight extension of the existing solver pipeline and should be done first. It should be evaluated both at the current preferred mass and near $m_{\rm res}$. If $A_{\rm eff}/A_{\rm unit}$ deviates significantly from 1, or if the $m_{\rm res}$ region looks promising, Step 3 becomes necessary.

### Caveat 3: Cross term $C_L^{\phi\tau}$ is set to zero without justification

Setting $\rho_L = 0$ (no correlation between ALP fluctuations and $\delta\tau$) is a simplification. Patchy reionization is sourced by matter perturbations, which also contribute to ALP clustering. The cross correlation is unlikely to be exactly zero. At minimum, an order-of-magnitude bound on the cross term using $|\rho_L| \le 1$ should be provided to show it does not qualitatively change the budget.

### Caveat 4: $\phi_{\rm amp,max}$ is a phenomenological proxy, not a rigorous bound

The maximum ALP amplitude is estimated from the energy-density condition $\rho_\phi \le \rho_{\rm DM}$. For $m_a \sim 10^{-27}$ eV, independent constraints exist from:

- Lyman-$\alpha$ forest,
- galaxy formation (substructure suppression),
- CMB small-scale power.

These can be significantly more restrictive than the raw energy-density bound. The current $\phi_{\rm amp,max}$ should be treated as an upper envelope, and the paper should flag that tighter astrophysical constraints would reduce the allowed amplitude.

---

## Suggested near-term paper framing

Possible title:
**Patchy reionization as a non-negligible source of anisotropic cosmic birefringence in axion-photon coupling models**

(The earlier "Can patchy reionization dominate...?" framing is less appropriate now that the natural-unit budget places the effect at the percent-to-ten-percent level rather than dominance.)

Suggested structure:
1. Introduction
2. Formalism: decomposition of $\alpha = \alpha_\phi + \alpha_\tau$
3. Effective patchy term and emit-time-shift approximation
4. ALP background dynamics and $A_{\rm unit}(m_a)$
5. Natural-unit normalization and physical $A_\tau$
6. Birefringence budget and reinterpretation of observational limits
7. Results: percent-to-ten-percent level contribution for benchmark couplings
8. Discussion: caveats, $C_L^{\tau\tau}$ modeling, cross term

Most likely central figures:

- $A_\tau^{\rm physical}(m_a)$ for benchmark $g_{a\gamma}$ values (from `21`)
- $D_L^{\alpha\alpha}$ budget compared to the anisotropic CB limit (from `22`)

The old $\phi_{\rm needed}$ vs $\phi_{\rm amp,max}$ plot (from `11`/`12`) should be retained as a supplementary figure showing that the code-unit amplitude budget is not violated, but it is not the main result.

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

The project has moved through feasibility, matched rerun, and natural-unit normalization stages.

What is already in hand:

- a stable unit-response workflow (`04a`, `04b`)
- a matched high-precision rerun confirming no solver-mixing artifact (`11`, `12`)
- natural-unit conversion of the matched background quantities (`21`)
- a $D_L^{\alpha\alpha}$ budget in natural units with explicit $g_{a\gamma}$ (`22`)

Current best physical result (from `21`/`22`):

- at matched mass $m_{\rm best} = 5.878016 \times 10^{-27}\ {\rm eV}$,
  $A_\tau^{\rm physical} \simeq 0.16$ for $g = 1.4 \times 10^{-12}\ {\rm GeV}^{-1}$,
  $A_\tau^{\rm physical} \simeq 0.45$ for $g = 4.0 \times 10^{-12}\ {\rm GeV}^{-1}$
- the patchy $D_L^{\alpha\alpha}$ contribution is at the **1–10% level** of the anisotropic CB observational limit for these benchmark couplings

The earlier code-unit ratio $\phi_{\rm needed}/\phi_{\rm amp,max} \sim 10^{-15}$–$10^{-11}$ from `11`/`12` confirmed amplitude feasibility but is **not a physical observable**; `21`/`22` are the correct reference for paper-level claims.

The next job is:

- use `22` as the canonical budget figure,
- reformulate the paper claim around percent-to-ten-percent non-negligibility rather than dominance,
- quantify the $C_L^{\tau\tau}$ template uncertainty and the emit-time-shift approximation error,
- and prioritize follow-up tests near $m_{\rm res} \sim 10^{-29}$ eV, where spatial resonance remains possible.
