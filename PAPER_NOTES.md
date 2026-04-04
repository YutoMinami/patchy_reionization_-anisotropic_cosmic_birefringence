# PAPER_NOTES.md

## いまの主張候補

このプロジェクトの現時点での主張候補は、次の3本です。

1. patchy reionization による effective birefringence term は、toy spectrum では genuine ALP term を十分に上回りうる
2. required amplitude $\phi_{\rm needed}(m_a)$ は、現在の phenomenological $\phi_{\rm amp,max}(m_a)$ と比べて非常に小さい
3. この「非常に小さい ratio」は、mixed solver setup の artifact ではなく、matched rerun (`11/12`) でも維持される

## 何が robust で、何がまだ toy か

### かなり robust な点

- $\delta\alpha_\tau$ という effective patchy term の定式化
- $\phi_{\rm ini}$ 小振幅での異常が数値 artifact だったこと
- unit-response rescaling を使うべきこと
- matched rerun でも $\phi_{\rm needed} / \phi_{\rm amp,max} \ll 1$ が維持されること

### まだ toy / phenomenological な点

- $C_L^{\tau\tau}$ の形と正規化
- $C_L^{\phi\phi}$ の toy power-law choice
- $\phi_{\rm amp,max}$ の energy-density interpretation
- finite-width reionization kernel を省いた emit-time-shift 近似

## 一番自然な paper storyline

1. 異方的 CB を $\alpha_\phi + \alpha_\tau$ に分解する
2. patchy term は effective emission-time shift として現れる
3. toy spectra で $R_\tau(L)$ を見ると、bubble-scale multipole 付近で patchy dominance がありうる
4. required amplitude $\phi_{\rm needed}(m_a)$ を作る
5. phenomenological amplitude bound $\phi_{\rm amp,max}(m_a)$ を作る
6. matched rerun で solver-mixing の疑いを潰す
7. その上でも $\phi_{\rm needed} / \phi_{\rm amp,max}$ は十分小さい
8. よって patchy contribution は current proxy では strongly viable

## figure 候補

### Main figures

1. $R_{\tau,\max}^{\rm unit}(m_a)$ または $\phi_{\rm needed}(m_a)$
2. matched $\phi_{\rm needed}(m_a)$ と matched $\phi_{\rm amp,max}(m_a)$ の比較
3. matched $\phi_{\rm needed} / \phi_{\rm amp,max}$

### support / appendix figures

1. `02` の linearity check
2. `03` の convergence check
3. `10` の constant-scaling check
4. old non-matched `08` と matched `12` の比較

## 読者への言い方

強く言ってよいこと:

- current toy + phenomenological framework では no-go ではない
- むしろ amplitude budget の観点ではかなり余裕がある
- solver mismatch を直しても結論は変わらない

慎重に言うべきこと:

- realistic $C_L^{\tau\tau}$ をまだ使っていない
- $\phi_{\rm amp,max}$ はまだ phenomenological proxy
- したがって final physical prediction ではなく “encouraging viability result” の段階

## 次に書くべき短い文章

### abstract 的な一文

We show that a patchy-reionization-induced effective birefringence term can dominate over the genuine ALP fluctuation contribution in toy spectra near the reionization bubble scale, and that the ALP amplitude required for such dominance remains far below a phenomenological energy-density bound even after a matched high-precision rerun.

### introduction 的な一文

The key question is not whether an effective patchy contribution exists formally, but whether it can be large enough to matter after the ALP background amplitude is physically constrained.

### discussion 的な一文

The remaining uncertainty is not the numerical stability of the background solver, but the degree to which the present toy $C_L^{\tau\tau}$ and phenomenological amplitude bound capture the realistic reionization and ALP parameter space.

## 次の実務

1. `README` / `RESEARCH_LOG` / `HANDOFF` を current best result にそろえる
2. matched `12` を見ながら paper figure 候補を固定する
3. $A_{\rm CB}$ 再解釈を matched result に基づいて書き直す
4. 必要なら realistic $C_L^{\tau\tau}$ の benchmark を1本だけ追加する
