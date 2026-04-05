# A Methodological Note on Interpreting Anisotropic Cosmic Birefringence Constraints in the Presence of Patchy Reionization

## Abstract

Current constraints on anisotropic cosmic birefringence are often interpreted as limits on the genuine ALP-induced fluctuation term alone. In this project we revisited that assumption in the presence of an effective patchy-reionization contribution. The basic point is simple: the observable anisotropic birefringence power spectrum is a sum of contributions, not a single component. Therefore, observational limits should in general be read as constraints on the total birefringence power, rather than on $C_L^{\phi\phi}$ alone [1,2].

We constructed a sequence of toy, matched, natural-unit, and visibility-weighted estimates to test whether the patchy term could become phenomenologically important. A thin-shell treatment suggested that the patchy contribution might reach the percent-to-ten-percent level of current anisotropic-CB limits for benchmark couplings. However, once finite-width visibility averaging is included, the old preferred mass point is strongly suppressed. The physically relevant mass window moves toward the resonance scale $m_{\rm res} \sim 10^{-29}\,{\rm eV}$, where oscillatory structure remains, but in the current surrogate setup the observable patchy signal stays below the present anisotropic-CB limit by a substantial margin. Even optimistic rescalings of the surrogate template normalization only lift the signal to at most the few-percent level for the larger benchmark coupling.

The main conclusion is therefore methodological rather than discovery-driven: anisotropic-CB constraints should be interpreted as applying to the total birefringence power, and any patchy contribution must be evaluated with a finite-width treatment before drawing physical conclusions.

## 1. Motivation

The starting question of this project was whether patchy reionization could induce an effective anisotropic cosmic birefringence term large enough to compete with, or even dominate over, the genuine ALP fluctuation term. In a thin-shell picture, the answer can easily look encouraging, because one is effectively probing the ALP time derivative at a single emission time. This naturally raises the possibility that a suitable matching between the ALP oscillation scale and the reionization scale could enhance the signal.

However, the observable quantity is not the thin-shell response itself. What matters is the visibility-weighted response after integrating across the finite thickness of reionization. The project therefore evolved from a feasibility study into a question of interpretation:

- What is actually constrained by current anisotropic-CB limits?
- How much of the apparent thin-shell signal survives finite-width averaging?
- Can one meaningfully separate a constraint on the genuine ALP term from a constraint on the patchy-induced effective term?

## 2. Decomposition of the Birefringence Signal

The anisotropic birefringence field is decomposed as

```math
\delta\alpha(\hat n)=\delta\alpha_\phi(\hat n)+\delta\alpha_\tau(\hat n),
```

where

```math
\delta\alpha_\phi(\hat n)
=
-\frac{g_{a\gamma}}{2}\,
\delta\phi(\bar\eta_{\rm emit}, \chi_{\rm emit}\hat n)
```

is the genuine ALP fluctuation term, and

```math
\delta\alpha_\tau(\hat n)
=
-\frac{g_{a\gamma}}{2}\,
\dot{\bar\phi}(\bar\eta_{\rm emit})
\frac{d\eta}{d\tau}\,
\delta\tau(\hat n)
```

is the effective patchy term induced by fluctuations in the optical-depth field.

The angular power spectrum therefore becomes

```math
C_L^{\alpha\alpha}
=
A_\phi^2 C_L^{\phi\phi}
+
A_\tau^2 C_L^{\tau\tau}
+
2 A_\phi A_\tau C_L^{\phi\tau}.
```

If $C_L^{\phi\tau}$ is neglected, observational constraints apply to

```math
A_\phi^2 C_L^{\phi\phi} + A_\tau^2 C_L^{\tau\tau},
```

not to $A_\phi^2 C_L^{\phi\phi}$ by itself. This is the central interpretational point.

## 3. What This Project Established

### 3.1 Thin-shell feasibility

The early toy analysis indicated that a patchy contribution can look large in a thin-shell approximation. In particular, the old preferred mass around

```math
m_{\rm best} \simeq 5.878016\times10^{-27}\,{\rm eV}
```

gave a strong unit-response signal. In natural units, this translated into a patchy contribution at roughly the 1–10% level of current anisotropic-CB limits for the benchmark couplings

```math
g_{a\gamma} = 1.4\times10^{-12}\ {\rm GeV}^{-1},
\qquad
4.0\times10^{-12}\ {\rm GeV}^{-1}.
```

These results are still useful, but they should be understood as thin-shell upper bounds. The coupling benchmarks are motivated by cluster/X-ray bounds in the ultra-light mass regime [4].

### 3.2 Finite-width visibility averaging

The crucial correction came from directly computing a visibility-weighted effective response $A_{\rm eff}$. Once this was done, the old thin-shell preferred point was found to be strongly suppressed:

```math
\left.\frac{A_{\rm eff}}{A_{\rm unit}}\right|_{m_{\rm best}}
\simeq 0.0148.
```

This means that the old thin-shell story at $m_{\rm best}$ is not physically robust.

At the resonance-scale mass

```math
m_{\rm res} = 1.023963\times10^{-29}\,{\rm eV},
```

the suppression is much weaker,

```math
\left.\frac{A_{\rm eff}}{A_{\rm unit}}\right|_{m_{\rm res}}
\simeq 0.734,
```

and the physically interesting window shifts toward the $m_{\rm res}$ region. A finer scan (201 points) showed that the response there is oscillatory and sign-changing as a function of mass, because the ALP oscillation phase at reionization shifts continuously with $m_a$. The current best point appears near $m/m_{\rm res} \simeq 0.895$.

![A_eff/A_unit vs m/m_res](figures/fig1_aeff_over_aunit.png)

*Figure 1. Visibility-averaged suppression factor $A_{\rm eff}/A_{\rm unit}$ as a function of $m/m_{\rm res}$ from the 201-point fine scan. The response is oscillatory and sign-changing. The vertical dotted line marks $m = m_{\rm res}$. The thin-shell limit corresponds to $A_{\rm eff}/A_{\rm unit} = 1$.*

### 3.3 Current observable impact

Using the visibility-weighted natural-unit pipeline, the current best patchy contribution remains small:

- about $2.15\times10^{-4}$ of the anisotropic-CB limit for $g_{a\gamma}=1.4\times10^{-12}\,{\rm GeV}^{-1}$
- about $1.75\times10^{-3}$ of the anisotropic-CB limit for $g_{a\gamma}=4.0\times10^{-12}\,{\rm GeV}^{-1}$

This already indicates that the finite-width effect is not a minor correction: it qualitatively changes the phenomenological conclusion.

![Budget vs limit](figures/fig2_budget_vs_limit.png)

*Figure 2. Ratio $D_L^{\alpha\alpha,\tau,{\rm eff}} / D_L^{\rm limit}$ as a function of multipole $L$ at the best visibility-weighted mass ($m \simeq 9.16\times10^{-30}$ eV, $m/m_{\rm res} \simeq 0.895$), for both benchmark couplings. Even at the peak multipole $L \simeq 289$, the patchy contribution is below the current anisotropic-CB limit by roughly three orders of magnitude.*

## 4. Normalization Sensitivity

The remaining major uncertainty is the normalization of the surrogate $C_L^{\tau\tau}$ template. To quantify this, we asked how much the template normalization would need to be boosted in order to raise the patchy contribution to selected fractions of the current anisotropic-CB limit. The surrogate family itself was inspired by the Dvorkin-Smith patchy-reionization framework [3].

For the current best visibility-weighted point, the required boosts are:

- low-$g$ benchmark:
  - about $4.66\times10^1$ to reach 1%
  - about $4.66\times10^2$ to reach 10%
- high-$g$ benchmark:
  - about $5.70$ to reach 1%
  - about $5.70\times10^1$ to reach 10%

We also tested steeper Dvorkin-Smith-inspired amplitude scalings of the form

```math
D_{\rm peak} \propto (A/A_{\rm fid})^{p_{\rm amp}} (b/b_{\rm fid})^{p_b}
```

with $p_{\rm amp}, p_b \le 3$. Within the currently explored parameter range, the maximum available boost was only

```math
{\rm max\ boost} \simeq 3.27\times10^1.
```

Combining this with the current best visibility-weighted budget gives an optimistic best-case estimate of

- about $7.0\times10^{-3}$ for the low-$g$ benchmark
- about $5.7\times10^{-2}$ for the high-$g$ benchmark

as fractions of the current anisotropic-CB limit.

Thus, even after fairly aggressive internal rescaling of the surrogate family, the patchy term does not automatically become a large observable signal.

## 5. Methodological Implication

The strongest claim supported by this project is not that patchy reionization generates a currently detectable anisotropic-CB signal. Rather, it is the following.

1. Constraints on anisotropic cosmic birefringence should be interpreted as constraints on the total birefringence power.
2. A patchy-induced effective term can in principle contribute to that total, and therefore reduces the room available for the genuine ALP term.
3. Thin-shell estimates can substantially overstate the relevance of the patchy term.
4. Finite-width visibility averaging is essential for any physically meaningful assessment.

This is enough to justify a methodological caution: one should not read current anisotropic-CB limits as direct limits on $C_L^{\phi\phi}$ alone unless the patchy contribution has first been shown to be negligible in a finite-width treatment [1,2].

## 6. Why This May Still Be Worth Keeping

Even if the project does not presently motivate a standalone discovery paper, it still provides a useful internal note for future work. It documents:

- which parts of the optimistic signal story were artifacts of the thin-shell treatment,
- where the physically relevant mass window actually moves once finite-width effects are included,
- how strongly the result depends on surrogate-template normalization,
- and how existing anisotropic-CB constraints should be reinterpreted when multiple birefringence sources contribute.

In that sense, the most durable output of this project is a clean cautionary framework for future anisotropic-CB analyses.

## 7. Bottom Line

The current project does **not** support a strong claim that patchy reionization generates a large observable anisotropic-CB signal under the present surrogate assumptions. What it *does* support is a more limited but solid methodological statement:

current anisotropic-CB constraints should be read as constraints on the total birefringence power, and any patchy contribution must be assessed with a finite-width visibility treatment before attributing the bound solely to the genuine ALP fluctuation term.

## References

1. T. Namikawa, "Reionization-era cosmic birefringence with finite last-scattering thickness," *Phys. Rev. D* **109**, 123521 (2024) [arXiv:2404.13771].
2. A. Lonappan, B. Keating, and K. Arnold, "Constraints on Anisotropic Cosmic Birefringence from CMB B-mode Polarization," arXiv:2504.13154.
3. V. Dvorkin and K. M. Smith, "Reconstructing Patchy Reionization from the Cosmic Microwave Background," *Phys. Rev. D* **79**, 043003 (2009).
4. M. Berg, J. P. Conlon, F. Day, N. Jennings, S. Krippendorf, A. J. Powell, and M. Rummel, "Constraints on Axion-Like Particles from X-ray Observations of NGC1275," *Astrophys. J.* **847**, 101 (2017) [arXiv:1605.01043].
