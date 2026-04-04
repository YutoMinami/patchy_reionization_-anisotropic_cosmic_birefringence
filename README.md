# PatchyReionization

patchy reionization が、標準的な ALP fluctuation 起源の anisotropic cosmic birefringence (CB) と競合する、あるいはそれを上回る有効寄与を生みうるかを調べる数値研究プロジェクトです。

## 目的

このプロジェクトでは、anisotropic CB の角度パワースペクトルを

```math
C_L^{\alpha\alpha}
=
A_\phi^2 C_L^{\phi\phi}
+
A_\tau^2 C_L^{\tau\tau}
+
2 A_\phi A_\tau C_L^{\phi\tau}
```

と分解して考えます。特に注目しているのは、patchy reionization に対応する

```math
A_\tau^2 C_L^{\tau\tau}
```

が、本来の ALP fluctuation 項

```math
A_\phi^2 C_L^{\phi\phi}
```

と同程度以上になりうるかどうかです。

実際の研究上の問いは、「ALP 背景場を物理的に正規化した後でも patchy 項が優勢になりうるか」です。

## $\alpha_\phi$ と $\alpha_\tau$ の分解

このプロジェクトでは、観測される anisotropic CB の揺らぎ $\alpha(\hat n)$ を

```math
\alpha(\hat n) = \alpha_\phi(\hat n) + \alpha_\tau(\hat n)
```

と分けて考えます。

- $\alpha_\phi$
  ALP fluctuation そのものに由来する、本来の anisotropic CB 成分
- $\alpha_\tau$
  patchy reionization による有効的な scattering time のずれを通して生じる成分

このとき角度パワースペクトルは

```math
C_L^{\alpha\alpha}
=
C_L^{\alpha_\phi \alpha_\phi}
+
C_L^{\alpha_\tau \alpha_\tau}
+
2 C_L^{\alpha_\phi \alpha_\tau}
```

と書けます。

さらに係数を明示すると、

```math
C_L^{\alpha\alpha}
=
A_\phi^2 C_L^{\phi\phi}
+
A_\tau^2 C_L^{\tau\tau}
+
2 A_\phi A_\tau C_L^{\phi\tau}
```

です。

それぞれの意味は次の通りです。

- $A_\phi^2 C_L^{\phi\phi}$
  genuine な ALP 起源の anisotropic CB
- $A_\tau^2 C_L^{\tau\tau}$
  patchy reionization による有効 birefringence 項
- $2 A_\phi A_\tau C_L^{\phi\tau}$
  両者が相関しているときに現れる cross term

現段階では、主な関心は $\alpha_\tau$ の自己相関項が $\alpha_\phi$ の自己相関項に勝てるかどうかにあります。

## 現在の理解

- 定式化自体は `HANDOFF.md` に整理されています。
- toy spectrum を使った初期検証では、bubble scale に対応する $L \sim 300$ 付近で patchy 項が優勢になりうることが示唆されています。
- ただし、これまでに出ていた非常に大きな $R_\tau$ は ALP 振幅の物理正規化前の値であり、そのまま物理的主張には使えません。
- 元の feasibility notebook では、$phi_ini$ を小さくすると応答係数 $A$ が線形に縮まず、場合によっては符号まで反転する問題がありました。
- この問題は現在、物理ではなく数値設定に由来するものだと分かっています。

## 安定な計算方針

応答係数は

```math
A(m_a) = \dot{\phi}_{\rm conf}(\eta_{\rm rei}) \frac{d\eta}{d\tau}
```

として評価します。

背景 ODE は $\phi$ に対して線形なので、推奨するワークフローは次の通りです。

1. $phi_ini = 1$ で一度だけ解く
2. その解から $A_unit(m_a)$ を定義する
3. 任意の振幅については後から線形に rescale する

```math
A(m_a; \phi_{\rm amp}) = \phi_{\rm amp} A_{\rm unit}(m_a)
```

これにより、極小振幅で何度も ODE を解くことによる不安定性を避けられます。

また、reionization での値の評価には sampled interpolation ではなく dense output を使う方針にしています。

## リポジトリ構成

- [HANDOFF.md](./HANDOFF.md)
  物理的背景、主要式、注意点、今後の優先課題。
- [RESEARCH_LOG.md](./RESEARCH_LOG.md)
  人間向けの研究ログ。何を試し、何が分かり、何が未解決かを追えるようにした記録。
- [patchy_reionization.py](./patchy_reionization.py)
  notebook と script で共有する数値計算ロジック。
- [scripts/README.md](./scripts/README.md)
  `scripts/` 配下の役割説明と推奨実行順。
- [scripts/01-check_feasibility.ipynb](./scripts/01-check_feasibility.ipynb)
  図を見ながら確認するための notebook。
- [scripts/01-check_feasibility.py](./scripts/01-check_feasibility.py)
  初期 feasibility scan の再現用スクリプト。
- [scripts/02-check_linearity.py](./scripts/02-check_linearity.py)
  `phi_ini` に対する線形性の検証。
- [scripts/03-check_convergence.py](./scripts/03-check_convergence.py)
  solver と tolerance に対する収束性の確認。

## 典型的な使い方

環境を用意:

```bash
uv sync
```

再現性のある script 実行:

```bash
uv run python scripts/02-check_linearity.py
uv run python scripts/03-check_convergence.py
uv run python scripts/01-check_feasibility.py --no-show
```

図を見ながら進める場合:

```bash
jupyter lab scripts/01-check_feasibility.ipynb
```

補足:

- `scripts/01-check_feasibility.py` はフルスキャンだとかなり時間がかかります。
- notebook は図を見るために残しています。
- 数値実装の正本は `patchy_reionization.py` と script 群です。

## 数値問題の現状

元の問題は、

- 小さい `phi_ini` に対して `A` が線形に縮まらない
- 場合によっては `A` の符号が反転する

というものでした。

現時点の診断では、

- tiny amplitude での繰り返し solve
- 固定絶対誤差が効きすぎる設定
- reionization 点での sampled interpolation

の組み合わせが原因でした。

dense output と unit-response rescaling を使う現在の方法では、線形性は回復しています。

## 現時点での `C_L^{\tau\tau}` の扱い

重要な点として、現段階では `\tau` fluctuation そのものの物理モデルはまだ入れていません。

つまり今の解析では、

- `\delta\alpha_\tau \propto \dot{\bar\phi}\,\delta\tau`
  という有効記述を使い
- `C_L^{\tau\tau}` は toy spectrum として与えています

したがって、現在の結果は

- patchy reionization の具体的 bubble model や radiative transfer model を仮定した予言

ではなく、

- 「もし `C_L^{\tau\tau}` がこれくらい大きければ、patchy 項はどれくらい効きうるか」

を見る feasibility study です。

この意味で、現状の主張は `\tau` の詳細モデルよりも、`C_L^{\tau\tau}` に対する必要条件や感度評価に近いものです。

## 次の科学的ステップ

1. `R_tau,max^unit(m_a)` を質量全域で評価する
2. `R_tau,max = 1` などを達成するために必要な `phi_needed(m_a)` を求める
3. ALP energy-density bound から許容される `phi_amp_max(m_a)` を見積もる
4. `phi_needed(m_a)` と `phi_amp_max(m_a)` を比較する
5. その後で必要なら cross spectrum の扱いを深める

## 重要な注意

toy model で得られる大きな `R_tau` は、feasibility を考える手がかりとしては有用ですが、それ自体が最終的な物理結果ではありません。最終的な中心結果は、「patchy 優勢に必要な振幅」と「物理的に許される振幅」の比較になる想定です。
