# PAPER_NOTES.md

## いまの主張候補

このプロジェクトの現時点での主張候補は、次の3本です。

1. patchy reionization による effective birefringence term は、toy spectrum では genuine ALP term を上回りうるが、physical normalization を入れると現状では percent-to-ten-percent level の寄与として現れる
2. required amplitude $\phi_{\rm needed}(m_a)$ は、matched rerun (`11/12`) でも phenomenological $\phi_{\rm amp,max}(m_a)$ より十分小さい
3. 既存の anisotropic CB constraint は genuine term 単独ではなく total birefringence power への制限として読むべきであり、patchy contribution は subdominant でも無視できない

## 何が robust で、何がまだ toy か

### かなり robust な点

- $\delta\alpha_\tau$ という effective patchy term の定式化
- $\phi_{\rm ini}$ 小振幅での異常が数値 artifact だったこと
- unit-response rescaling を使うべきこと
- matched rerun でも $\phi_{\rm needed} / \phi_{\rm amp,max} \ll 1$ が維持されること
- natural-unit へ直した `21/22` で、huge な $C_L^{\alpha\alpha}$ が主に unit-system mismatch に由来していたと分かったこと

### まだ toy / phenomenological な点

- $C_L^{\tau\tau}$ の形と正規化
- $C_L^{\phi\phi}$ の toy power-law choice
- $\phi_{\rm amp,max}$ の energy-density interpretation
- finite-width reionization kernel を省いた emit-time-shift 近似
- $g_{a\gamma}$ を抜いた raw response と physical observable の切り分け

補足:

- $A_{\rm unit}$ と $\phi_{\rm amp}$ から作る quantity は、まだ $g_{a\gamma}$ を抜いた raw response である
- $R_\tau$ のような ratio ではそれでよいが、
  absolute な $C_L^{\alpha\alpha}$ budget と比べるときは
  $A_\tau = -(g_{a\gamma}/2)\phi_{\rm amp} A_{\rm unit}$ を使う必要がある
- したがって、anisotropic CB limit との比較では
  coupling を入れた physical normalization を明示しないと overclaim になりやすい

## anisotropic CB constraint をどう読むか

ここは論文の主張としてかなり大事な部分で、まず一般論と benchmark 依存部分を分けて書くのがよい。

### 一般論として言えること

観測が直接制限しているのは genuine ALP term 単独ではなく、全 birefringence power

```math
C_L^{\alpha\alpha}
=
A_\phi^2 C_L^{\phi\phi}
+
A_\tau^2 C_L^{\tau\tau}
+
2 A_\phi A_\tau C_L^{\phi\tau}
```

である。

したがって、もし patchy-induced effective term が無視できないなら、
「既存の anisotropic CB constraint は $C_L^{\phi\phi}$ にそのままかかっている」
という読みは一般には正しくない。

特に $C_L^{\phi\tau}=0$ を置いても、観測が制限するのは

```math
A_\phi^2 C_L^{\phi\phi} + A_\tau^2 C_L^{\tau\tau}
```

の和であり、$C_L^{\tau\tau}$ が大きいほど genuine term に許される余地は小さくなる。

### まだ benchmark が必要なこと

一方で、数値的に「どれだけ制限が強くなるか」は $C_L^{\tau\tau}$ の形や正規化に依る。

つまり、

- 一般論としては「和に制限がかかる」と言える
- しかし数値的な bound の強さを言うには representative な $C_L^{\tau\tau}$ template が必要

という二段構えになる。

この意味で、論文の主張として最も安全なのは

1. 既存 anisotropic CB constraint は total birefringence power への制限として読むべき
2. patchy term があれば、genuine ALP term に使える budget は減る
3. representative templates では、この再解釈は数値的にも重要になりうる

という形である。

### 何を避けるべきか

現段階では、realistic $C_L^{\tau\tau}$ をまだ確定していないので、

- model-independent を装って specific $C_L^{\tau\tau}$ bound を強く主張すること
- 既存の $A_{\rm CB}$ constraint を無条件に $C_L^{\tau\tau}$ 単独へ写すこと

は避けた方がよい。

むしろ、

- 本体では general reinterpretation
- 後半では benchmark/template-based illustration

という構成が自然である。

## patchy $C_L^{\tau\tau}$ template をどう置くか

ここは次の具体作業として重要だが、現時点では「独立な複数モデルを並べる」より、
まず Dvorkin-Smith (2009) 系の patchy-$\tau$ family を 1 つ採用し、その中で
代表パラメータを振る方が自然である。

Roy などの後続研究も、少なくとも現段階の理解では Dvorkin-Smith と無関係な
全く別 family というより、同系統の bubble-based prescription を更新したものと
見なすのがよい。

したがって、paper 上の safest framing は

1. Dvorkin-Smith 系を representative patchy template family として採用する
2. bubble radius, duration, normalization などを振って幅を見る
3. Roy などの後続文献は「同系統 family に対する updated astrophysical input / 振幅感」
   の参照として使う

という形である。

## 実装方針: いきなり厳密積分へ行かない

見たところ、realistic な $C_L^{\tau\tau}$ は line-of-sight 的な積分や bubble-model の
積分を通して出てくるため、そのまま厳密に実装するとかなり重い可能性が高い。

このため、次の段階ではいきなり full integral を実装するのではなく、

- 文献の代表図や代表スケールに合わせた lightweight template family
- peak multipole, width, normalization を持つ phenomenological surrogate
- 必要なら低L・高L の傾きを文献風に寄せた broken-power-law / log-normal 近似

をまず試すのがよい。

つまり、

- 第1段階: literature-inspired but lightweight template family
- 第2段階: 必要なら Dvorkin-Smith 系の積分実装
- 第3段階: その後に visibility Gaussian 近似で emit-time-shift の限界を点検

の順で進めるのが、安全で再現性も高い。

`14` ではこの方針に沿って、

- $D_L^{\tau\tau}$ を unit-peak の log-normal bump family で近似し
- $L_{\rm peak}$ と $\sigma_{\ln L}$ を振り
- 観測 window にどれだけ重なるかで allowed normalization がどう変わるか

を見るところまで進めた。

ここでの主張は、

- exact Dvorkin-Smith template を再現した
ではなく、
- toy Gaussian 1本だけに依らず、文献系 patchy family を意識しても
  anisotropic CB reinterpretation は依然かなり強い

という、family-level の傾向である。

ただし `14` の unit-peak normalization は、shape 依存性を見る補助結果としては useful でも、
振幅の任意性が大きすぎて本命の physical claim にはしにくい。

したがって次の本命は、

- Dvorkin-Smith の `B. Reionization model parameters` の整理に従い
- observable amplitude $A$ で振幅を固定し
- Eq.(78) の $R_{\rm eff} = \bar R \exp(4\sigma_{\ln R}^2)$ で代表スケールを与え
- その上で shape だけを log-normal surrogate で軽量化する

という $A + R_{\rm eff}$ surrogate である。

この方針なら、

- `14` の「shape family を見る」利点を残しつつ
- 振幅は free normalization ではなく文献の observable combination に縛る

ことができる。

## 一番自然な paper storyline

1. 異方的 CB を $\alpha_\phi + \alpha_\tau$ に分解する
2. patchy term は effective emission-time shift として現れる
3. toy spectra で $R_\tau(L)$ を見ると、bubble-scale multipole 付近で patchy dominance がありうることを示す
4. その後、required amplitude $\phi_{\rm needed}(m_a)$ と phenomenological amplitude bound $\phi_{\rm amp,max}(m_a)$ を matched rerun で比較する
5. natural-unit への再解釈で unit-system mismatch を取り除く
6. 既存 anisotropic CB constraint は total birefringence power にかかることを整理する
7. Chandra benchmark coupling では patchy contribution は percent-to-ten-percent level で、subdominant だが無視できないと結論する

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
- amplitude budget の観点では no-go ではなく、natural-unit でも現実的 coupling で patchy contribution は残る
- solver mismatch を直しても結論は変わらない

慎重に言うべきこと:

- realistic $C_L^{\tau\tau}$ をまだ使っていない
- $\phi_{\rm amp,max}$ はまだ phenomenological proxy
- したがって final physical prediction ではなく “subdominant but non-negligible benchmark result” の段階

## 次に書くべき短い文章

### abstract 的な一文

We show that a patchy-reionization-induced effective birefringence term can provide a subdominant but non-negligible contribution to anisotropic cosmic birefringence, and that existing constraints should be interpreted as limits on the total birefringence power rather than on the genuine ALP fluctuation term alone.

### introduction 的な一文

The key question is not whether an effective patchy contribution exists formally, but whether it can be large enough to matter after the ALP background amplitude is physically constrained.

### introduction 的な別案

Existing anisotropic cosmic-birefringence constraints should be interpreted as limits on the total birefringence power, not automatically as limits on the genuine ALP fluctuation term alone, if patchy reionization induces an additional effective contribution.

### discussion 的な一文

The remaining uncertainty is not the numerical stability of the background solver, but the degree to which the present toy $C_L^{\tau\tau}$ and phenomenological amplitude bound capture the realistic reionization and ALP parameter space.

## 次の実務

1. `README` / `RESEARCH_LOG` / `HANDOFF` を current best result にそろえる
2. matched `12` を見ながら paper figure 候補を固定する
3. $A_{\rm CB}$ 再解釈を matched result に基づいて書き直す
4. Dvorkin-Smith 系 inspired の lightweight $C_L^{\tau\tau}$ family を1つ用意する
5. unit-peak ではなく $A + R_{\rm eff}$ surrogate に置き換える
6. その surrogate で anisotropic CB reinterpretation をやり直す
7. 必要ならその後に full integral か visibility Gaussian approximation に進む

## $21/22$ 以後の framing

`19` までの raw / naive physical budget は、結局のところ unit-system mismatch をかなり含んでいた。
したがって、これ以後の paper claim は $21/22$ を正本にして組み立てるべきである。

具体的には、matched mass 付近で natural-unit に写すと

- $A_\tau \simeq 0.157$ for $g = 1.4 \times 10^{-12}\,{\rm GeV}^{-1}$
- $A_\tau \simeq 0.448$ for $g = 4.0 \times 10^{-12}\,{\rm GeV}^{-1}$

となり、`22` の $D_L^{\alpha\alpha}$ budget では patchy term の最大は anisotropic-CB limit の
約 $1.2\%$ から $10\%$ に収まる。

したがって、少なくとも Chandra benchmark coupling の範囲では、

- patchy term は absurdly excluded ではない
- むしろ percent-to-ten-percent level の budget を使う可能性がある
- 以後の主張は「巨大すぎる patchy term」ではなく「見落とされやすい subdominant だが non-negligible な patchy contribution」

として組み立てる方が自然である。

## 次の実務の更新

1. `22` を現在の正本 budget plot として採用する
2. `23` を踏まえて、preferred mass では temporal averaging が主要 caveat であり、spatial resonance の本命は $m_{\rm res} \sim 10^{-29}\,{\rm eV}$ 付近だと明記する
3. 次の検証は $m_{\rm res}$ 近傍を中心に進める
4. anisotropic-CB reinterpretation は raw quantity ではなく natural-unit quantity でやり直す
5. Dvorkin-Smith 系 template の normalisation をもっと物理的に寄せる
6. 必要なら isotropic CB benchmark と組み合わせて $g_{a\gamma}$ のベンチマークを追加する
