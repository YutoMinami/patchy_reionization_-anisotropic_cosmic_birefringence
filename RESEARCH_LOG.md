# RESEARCH_LOG.md

## 目的

このファイルは `PatchyReionization` プロジェクトの人間向け研究ログです。
何を試し、何が分かり、どこに注意が必要で、次に何をやるべきかを、後から読み返しやすい形で残すことを目的にしています。

## 研究課題

中心的な問いは次の通りです。

patchy reionization に由来する anisotropic cosmic birefringence の有効寄与は、通常の ALP fluctuation 寄与と同程度以上になりうるか。

実際には

```math
R_\tau(L) = \frac{A_\tau^2 C_L^{\tau\tau}}{A_\phi^2 C_L^{\phi\phi}},
\qquad
R_\times(L) = \frac{2 A_\phi A_\tau C_L^{\phi\tau}}{A_\phi^2 C_L^{\phi\phi}}
```

を見ていますが、現段階で主に追っているのは `R_tau` です。

## `\alpha_\phi` と `\alpha_\tau` の見方

このプロジェクトでは、異方的な CB の揺らぎを

```math
\alpha(\hat n) = \alpha_\phi(\hat n) + \alpha_\tau(\hat n)
```

と分けて整理しています。

- `\alpha_\phi`
  ALP fluctuation に直接対応する成分
- `\alpha_\tau`
  patchy reionization が有効的な散乱時刻の揺らぎを作り、それが時間依存する ALP 背景場を通して CB に写る成分

そのため、角度パワースペクトルは

```math
C_L^{\alpha\alpha}
=
C_L^{\alpha_\phi \alpha_\phi}
+
C_L^{\alpha_\tau \alpha_\tau}
+
2 C_L^{\alpha_\phi \alpha_\tau}
```

と分解できます。

これを係数つきで書くと

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

この分解のうち、

- 第1項は本来の ALP 起源の anisotropic CB
- 第2項は patchy reionization による有効項
- 第3項は両者の相関があるときの cross term

という解釈になります。

現段階では、まず第2項が第1項に勝てるかどうかを見るのが中心課題で、cross term はまだ二次的な扱いです。

## `\delta\alpha_\phi` と `\delta\alpha_\tau` の式

上の分解を、実際の揺らぎの式として書くと

```math
\delta\alpha(\hat n)=\delta\alpha_\phi(\hat n)+\delta\alpha_\tau(\hat n)
```

です。

ここで、

- `genuine ALP term`

```math
\delta\alpha_\phi(\hat n)
=
-\frac{g_{a\gamma}}{2}\,
\delta\phi(\bar\eta_{\rm emit}, \chi_{\rm emit}\hat n)
```

- `effective patchy term`

```math
\delta\alpha_\tau(\hat n)
=
-\frac{g_{a\gamma}}{2}\,
\dot{\bar\phi}(\bar\eta_{\rm emit})
\frac{d\eta}{d\tau}\,
\delta\tau(\hat n)
```

と書けます。

解釈としては、

- `\delta\alpha_\phi` は ALP fluctuation そのものが作る本来の anisotropic CB
- `\delta\alpha_\tau` は patchy reionization が effective emission time をずらし、その時間ずれを時間依存する背景 ALP 場が CB のゆらぎへ変換したもの

です。

この見方をすると、`C_L^{\alpha\alpha}` の分解は単なる記号操作ではなく、「本来の ALP 項」と「patchy reionization を介して現れる有効項」の足し合わせとして理解できます。

## ここまでで分かっていること

### 定式化の面

一次の近似では、patchy reionization による有効 birefringence fluctuation は

```math
\delta \alpha_\tau(\hat n)
\simeq
-\frac{g_{a\gamma}}{2}
\dot\phi(\eta_{\rm rei})
\frac{d\eta}{d\tau}
\delta\tau(\hat n)
```

と書けます。

解釈としては、

- patchy reionization 自体が直接偏光面を回すわけではない
- しかし有効的な散乱時刻を方向依存でずらす
- 時間依存する ALP 背景場があると、そのずれが有効的な anisotropic CB になる

という構図です。

### toy spectrum の面

toy な `C_L^{\tau\tau}` と `C_L^{\phi\phi}` を使うと、`L ~ 300` 近傍で patchy 項が優勢になる状況は作れます。

これは feasibility の観点では重要です。
ただし、まだ物理正規化前なので、この時点では「そういう可能性が toy setup ではある」としか言えません。

### `\tau` fluctuation 自体のモデルはまだ入れていない

ここで注意すべきなのは、現段階では `\tau` fluctuation そのものに対して具体的な patchy reionization model をまだ入れていないことです。

いまやっているのは、

- `\delta\alpha_\tau \propto \dot{\bar\phi}\,\delta\tau`
  という有効関係を使う
- その上で `C_L^{\tau\tau}` は toy spectrum として与える

という段階です。

したがって、現在の `R_\tau` や `phi_needed` の結果は

- 特定の bubble model の予言
- realistic radiative transfer calculation の出力

ではなく、

- 「もし `C_L^{\tau\tau}` がある程度大きいなら、patchy 項はどれくらい効きうるか」
- 「そのとき必要な ALP 側の応答や振幅はどれくらいか」

を見るための feasibility / sensitivity study です。

これは弱点であると同時に利点でもある。

- 弱点:
  realistic `C_L^{\tau\tau}` の正規化や形状はまだ保証していない
- 利点:
  特定モデルに依存する前に、`C_L^{\tau\tau}` に対してどの程度の大きさが必要かを先に整理できる

今後もし patchy 項が本当に有望だと分かれば、その次の段階で realistic `C_L^{\tau\tau}` model を入れるのが自然である。

## 見つかっていた数値的問題

元の feasibility notebook では、`phi_ini` を変えながら毎回 ODE を解いていました。
すると `phi_ini` を十分小さくしたときに、

- `A` が `phi_ini` に比例して縮まらない
- ときには `A` の符号が変わる

という不自然な挙動が見られました。

背景 ODE は `phi` に対して線形なので、これは物理より数値の問題を疑うべき状況でした。

## 現在の診断

この問題は数値的なものです。

不安定化の主な原因は、

- tiny amplitude に対して何度も solve していたこと
- fixed absolute tolerance の影響が相対的に大きくなったこと
- reionization 点で sampled interpolation を使っていたこと

の組み合わせだと考えています。

## 現在の修正方針

現在採用している安定な流れは次の通りです。

1. `phi_ini = 1` で一度だけ解く
2. そこから `A_unit(m_a)` を作る
3. reionization での評価には dense output を使う
4. 振幅依存性は後から線形に rescale する

```math
A(m_a; \phi_{\rm amp}) = \phi_{\rm amp} A_{\rm unit}(m_a)
```

この方針は、共有モジュールと script 群に反映済みです。

## いまのファイル構成

- `patchy_reionization.py`
  共有の数値計算ロジック。
- `scripts/01-check_feasibility.py`
  初期 mass scan と toy `R_tau(L)` を再現する script。
- `scripts/02-check_linearity.py`
  `phi_ini` に対する線形性の確認。
- `scripts/03-check_convergence.py`
  solver / tolerance に対する収束確認。
- `scripts/01-check_feasibility.ipynb`
  図を見ながら検証する notebook。

## 最近の重要な結果

### 結果 1: 安定版 feasibility scan

安定版の script で feasibility scan を回したところ、定性的には以前の見通しを再現しました。

- ある質量帯で応答が大きくなる
- toy `R_tau` は非常に大きくなりうる
- 現在の toy setup では peak は引き続き `L = 356` 付近に出る

代表的な一例としては、

- `m_pick ~ 5.36e-28 eV`
- `A_unit ~ 8.04e9`
- `max R_tau ~ 7.79e19`
- peak は `L = 356`

という値が出ています。

ただし、これは toy model の結果であり、まだ物理正規化済みの結論ではありません。

### 結果 2: 線形性チェック

新しい線形性チェックでは、

- `legacy_A`: 元の notebook 風の取り方
- `robust_A`: dense output を使った安定版
- `unit_rescaled`: `A_unit * phi_ini`

を比較しました。

その結果、

- `robust_A` は `unit_rescaled` と実質的に一致
- `legacy_A` は系統的にずれる

ことが確認されました。

つまり、以前の small-amplitude pathology は物理ではなく、数値実装側の問題でした。

### 結果 3: 収束チェック

代表質量について、

- より厳しい `DOP853` を reference
- 中程度の `DOP853`
- 中程度の `RK45`
- loose な `RK45`

を比較しました。

現時点の印象としては、

- 軽い質量から中間質量まではかなり安定
- 問題になっていた質量帯でも `DOP853` は十分近い
- `RK45` でも使えなくはないが、高質量側では `DOP853` の方が安心

という状況です。

### 結果 4: `A_unit(m_a)` を `10^{-26} eV` まで延長

`A_unit(m_a)` の mass scan を `10^{-26} eV` まで延長した。

この作業では、高質量端がかなり重く、combined script では WSL2 側の資源制限に引っかかる疑いがあった。そのため workflow を次の 2 段に分割した。

- `scripts/04a-scan_a_unit.py`
  `A_unit(m_a)` だけを計算し、`results/04a-a-unit/A_unit_scan.csv` に逐次保存する
- `scripts/04b-scan_rtau_from_aunit.py`
  保存済み `A_unit` を読み込んで `R_{\tau,\max}^{unit}(m_a)` と `phi_needed(m_a)` を軽く後処理する

この分割は実務上かなり重要で、重い solve と軽い post-processing を切り離せるようになった。

`A_unit_scan.csv` の高質量端の代表値は例えば次の通り。

- `m = 1.1937766417144358e-27 eV`, `A_unit = 1.0127134845424482e10`
- `m = 2.0309176209047306e-27 eV`, `A_unit = 1.10500484397749e10`
- `m = 5.878016072274924e-27 eV`, `A_unit = 1.3935403275476799e10`
- `m = 1e-26 eV`, `A_unit = -2.632126036501822e9`

ここから分かることは、

- 少なくとも今回の粗い mass grid では、`A_unit` は `10^{-27} eV` 付近から `10^{-26} eV` にかけて単純に消えていくわけではない
- むしろ大きな応答が続きつつ、位相の都合で符号が振れる

ということ。

### 結果 5: `R_{\tau,\max}^{unit}(m_a)` と `phi_needed(m_a)` の更新

completed `A_unit` scan を使って `04b` を更新した結果、現時点では全 toy case で best mass は

```math
m_{\rm best} = 5.878016 \times 10^{-27}\ {\rm eV}
```

となった。

代表的な結果は以下。

- `nphi_1_no_cut`
  - `Rtau_max_unit = 2.042154e+20`
  - `L_peak = 330`
- `nphi_2_no_cut`
  - `Rtau_max_unit = 2.337756e+20`
  - `L_peak = 356`
- `nphi_3_no_cut`
  - `Rtau_max_unit = 2.865968e+20`
  - `L_peak = 379`
- `nphi_3_Lcut_800`
  - `Rtau_max_unit = 3.609159e+20`
  - `L_peak = 389`

特に `nphi_2_no_cut` の代表値では、

- `phi_needed(R = 0.1) = 2.068237474586814e-11`
- `phi_needed(R = 1) = 6.540341161808949e-11`
- `phi_needed(R = 10) = 2.0682374745868143e-10`

となった。

この段階で重要なのは、

- `1e-26 eV` まで見ても best mass は `5.878016e-27 eV` に残ったこと
- `R_{\tau,\max}^{unit}` が `10^{20}` 級まで大きくなること
- その結果、必要振幅 `phi_needed` が `10^{-11}` から `10^{-10}` 級まで下がっていること

である。

### 結果 6: 可視化 notebook の追加

保存済み CSV をそのまま 2 次元グラフで見られるように、

- `scripts/05-visualize_phi_needed.ipynb`

を追加した。

この notebook では次を可視化できる。

- `A_unit(m_a)` の絶対値と符号つきの形
- toy case ごとの `R_{\tau,\max}^{unit}(m_a)`
- `R_target = 0.1, 1, 10` に対する `phi_needed(m_a)`
- best mass の一覧と、その位置を示した比較図

人間が状況を理解するための可視化入口として重要。

### メモ: `A_unit(m_a)` の高質量側 plateau について

`A_unit(m_a)` が `10^{-31} eV` 以上でほぼフラットに見えるのは、今のところ不自然というより、むしろかなり自然な可能性があります。

見ている量は

```math
A_{\rm unit}(m_a) \propto \dot\phi_{\rm conf}(\eta_{\rm rei})
```

です。

`m_a \gtrsim H_{\rm rei}` では、reionization の時点で ALP 背景場はすでに振動相に入っていると考えられます。このとき大まかには

- `\phi` の包絡線は `a^{-3/2}` で減衰する
- `\dot\phi_{\rm phys} \sim m_a \phi`
- `\dot\phi_{\rm conf} = a \dot\phi_{\rm phys}`

となります。

さらに、振動開始時刻はおおよそ `H \sim m_a` で決まり、matter-dominated 近似では `a_{\rm osc}` の `m_a` 依存が入ることで、`m_a` の依存性がかなり相殺されます。その結果、固定された `a_{\rm rei}` で見た `\dot\phi_{\rm conf}`、ひいては `A_unit(m_a)` の包絡線が高質量側で平坦に近づいても不思議ではありません。

現時点での解釈は次の通りです。

- `m_a \ll H_{\rm rei}` では凍結に近く、質量依存性が強く出やすい
- `m_a \gtrsim H_{\rm rei}` では振動相に入り、`A_unit(m_a)` の包絡線は平坦化しやすい

今回の `H_{\rm rei}` のオーダーを考えると、「`10^{-31} eV` 以上でほぼフラット」という見え方は、少なくとも第一印象としては理屈に合っています。

ただし、これを物理的 plateau と断定する前に、次の確認はしておきたいです。

1. mass range を `10^{-26} eV` まで広げても包絡線として平坦か
2. `1e-29`, `1e-28`, `1e-27`, `1e-26 eV` 付近で収束性が十分か
3. 高質量端でギザギザ構造が数値 aliasing ではないか

## 何が分かっていて、何がまだ分からないか

### 分かっていること

- patchy 項は定式化として確かに現れる
- toy spectrum では `L ~ 300` 近傍で優勢になりうる
- `phi_ini` の異常な振る舞いは数値起源だった
- 安定な `A_unit(m_a)` ベースのワークフローが用意できた

### まだ分からないこと

- 物理正規化後にも patchy 優勢が残るか
- `phi_needed(m_a)` がどのような形になるか
- `phi_amp_max(m_a)` がどのような制限を与えるか
- その両者に overlap があるか

## 次にやるべきこと

1. `R_tau,max^unit(m_a)` を質量全域で計算する
2. `R_tau,max = 1`, `0.1`, `10` などに必要な `phi_needed(m_a)` を出す
3. ALP energy density から `phi_amp_max(m_a)` を見積もる
4. `phi_needed(m_a)` と `phi_amp_max(m_a)` を比較する
5. その後必要なら cross term の扱いを強化する

## TODO: emit-time shift 近似の有効条件と限界

現在の整理では、patchy 項を

```math
\delta\alpha_\tau(\hat n)
\simeq
-\frac{g_{a\gamma}}{2}
\dot{\bar\phi}(\bar\eta_{\rm emit})
\frac{d\eta}{d\tau}\,
\delta\tau(\hat n)
```

という「実効的に emit time がずれる」という近似で扱っている。

ただし、この近似が本当に十分かどうかは、まだ独立に検証すべき課題として残っている。

特に、このプロジェクトの元の発想には

- ALP 背景の振動周期
- patchy reionization の厚さ、あるいは visibility kernel の幅

がちょうどよいときに信号が強く見えるかもしれない、という見通しがある。

もし ALP の時間変化スケールが reionization の幅と同程度なら、単純な time-shift 近似だけでは不十分で、本当は有限幅の kernel で時間方向に平均した扱いが必要かもしれない。

したがって、将来的には次を確認したい。

1. どの条件で emit-time shift 近似が有効か
2. どの条件で有限幅 kernel の効果が重要になるか
3. 厚さと振動周期のマッチングで増強または平均化による抑制が起こるか

ただし、これは「今すぐ解決しないと何も言えない」という種類の問題ではない。

もし実際に `C_L^{\tau\tau}` が十分大きく、patchy 項が観測的に無視できないことが分かれば、

- patchy reionization 由来の有効項それ自体
- それに基づく `C_L^{\tau\tau}` への制限や感度

だけでも論文として十分価値がある可能性が高い。

なので優先順位としては、

- まず現在の有効理論的な整理で patchy 項の大きさを評価する
- その後、必要に応じて finite-width reionization の扱いを精密化する

という順が自然である。

## 実務上のメモ

- フルの mass scan はかなり時間がかかる
- notebook は図を見るために有用
- script は再現性と比較のしやすさで有用
- コアの数値ロジックは shared module に寄せるのがよい

## ひとことで言うと

このプロジェクトは、patchy 項が形式的に存在するかどうかを確かめる段階は越えつつあります。
いま本当に重要なのは、「物理的に許される ALP 振幅のもとで patchy 優勢が実現可能か」を調べることです。

`A_unit(m_a)` まわりの数値整理が済んだことで、次の段階は solver artifact ではなく physics の問題として進められる状態になりました。
