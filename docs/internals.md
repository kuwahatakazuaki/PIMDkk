# 開発者向け内部仕様

この章は、PIMD8の通常PIMD/CLMDについて、入力から終了までの主要処理をソース上で追うための案内です。各数値積分式の理論導出ではなく、責務とデータの流れを扱います。

## 1. 全体フロー

```text
main.F90
  ├─ MPI初期化、rank/process数取得
  ├─ read_parameter        input.inpのパラメータと原子数
  ├─ Broad1               scalar設定をworkerへ配布
  ├─ Set_Allocate         simulation設定に応じた配列確保
  ├─ read_structure       元素、座標、質量
  ├─ Broad2               構造をworkerへ配布
  ├─ Calc_Constant        単位変換、β、自由度、体積
  ├─ Check_Inp            std.outへ実効入力を表示
  ├─ Isimulationでdispatch
  │    ├─ 0..2 → PI_NEW_MPI
  │    ├─ 3    → PIHMC_normal
  │    └─ 10   → Classical
  ├─ MACE interface終了処理
  ├─ Set_Deallocate
  └─ MPI終了、終了時刻表示
```

入口は [`main.F90`](../main.F90) です。rank 0だけがファイル入力を行い、MPI buildでは設定と構造をbroadcastします。

## 2. 入力と初期化

### `read_parameter`

[`read_input.F90`](../read_input.F90) は `input.inp` を1行ずつ走査し、タグ行の次の行を値として読みます。`$end parameter` がない場合は停止します。その後 `$Coords` ブロック内の行数を数えて `Natom` を決めます。

ここで行う主な正規化・検査は次です。

- `Isimulation=0..2` かつ `Nbead=1` をCLMDへ変更
- simulation名の設定
- `Ndim` が1または3か検査
- PIHMCの `Ndyn` とdual targetを検査
- Result/Scratchディレクトリを作成

### 配列確保

[`Set_Allocate.F90`](../Set_Allocate.F90) は、座標・力・質量・backend出力用配列を確保します。CLMD以外ではnormal-mode変換行列、ring-polymer内力、非centroid thermostat配列も確保します。centroid thermostat配列は `Ncent` に応じて形が変わります。

### 構造と単位

`read_structure` は元素記号から質量を得て、`Langstrom=.T.` の座標をÅからbohrへ変換します。[`Calc_Constant.F90`](../Calc_Constant.F90) はさらに次を行います。

- `dt`: fsからatomic timeへ変換
- 質量: amuからelectron mass基準のatomic unitへ変換
- `beta`, 自由度、ring-polymer周波数を計算
- `lattice` の行列式から体積を計算
- Result/Scratchの内部パスを組み立てる

## 3. PIMD処理

driverは [`PI_NEW_MPI.F90`](../PI_NEW_MPI.F90) です。

### 初期化順序

1. `Setup_time_mass`: `dt_ref`、thermostat mass、Yoshida–Suzuki weightを設定
2. `set_pallarel`: 各MPI rankへ連続したビーズ範囲 `Ista:Iend` を割当
3. `Normal_Mode`: ring-polymer normal-mode変換行列を作成
4. `Init_Mass`: simulation種別に応じて仮想質量を設定
5. `set_Iforce`: backend固有の入力・Scratch・interfaceを準備
6. 新規計算なら初期速度、bath、初期力を計算。restartなら状態を読んで配布

### 1 stepの流れ

```text
centroid thermostat half step
  → physical velocity update
  → Nref回の内側積分
       thermostat / ring-polymer force / position / ring-polymer force
  → normal mode座標を各ビーズの実座標へ変換
  → backendで各ビーズのenergyとforceを評価
  → 実座標forceをnormal-mode forceへ変換
  → physical velocity update
  → centroid thermostat half step
  → Hamiltonian・温度・virial・圧力を評価
  → 条件付き出力とrestart更新
```

ring-polymer spring forceとビーズ平均距離拘束は [`Getforce_Ref.F90`](../Getforce_Ref.F90) で `fur_ref` に構成されます。backend由来の力は別に評価され、normal-modeへ変換されます。

### MPIでの力計算

rank 0が時間積分を進めます。各stepで [`Force_New_MPI.F90`](../Force_New_MPI.F90) が座標を担当rankへscatterし、rankごとに `Ista:Iend` のビーズを評価した後、force・energy・付随量をgatherします。worker rankは同じstep数だけ力計算呼出しを繰り返します。

## 4. CLMD処理

driverは [`Classical.F90`](../Classical.F90) です。`Nbead=1` のためnormal-mode変換とring-polymer内力は使いません。

1 stepは、centroid thermostat half step、velocity update、position update、[`Force_Classical.F90`](../Force_Classical.F90) による力評価、velocity update、thermostat half step、energy・出力の順です。

PIMDの `Nref` による内側分割はCLMDでは使われず、`dt_ref=dt` です。CLMDの力dispatchはPIMD側と完全には同一ではないため、利用可能なbackendは[入力リファレンス](input-reference.md#5-iforce-backend一覧)で確認してください。

## 5. backend接続

### 共通dispatch

初期化時は [`set_Iforce.F90`](../set_Iforce.F90)、PIMDの評価時は [`Force_New_MPI.F90`](../Force_New_MPI.F90)、CLMDの評価時は [`Force_Classical.F90`](../Force_Classical.F90) が `Iforce` をdispatchします。

backendは最終的に次の共有状態を設定します。

| 変数 | 内容 | 内部単位 |
|---|---|---|
| `pot_bead(imode)` | ビーズ別potential energy | Hartree |
| `fr(:,iatom,imode)` | ビーズ別Cartesian force | Hartree/bohr。PIMD側はビーズ平均係数を含む |
| `W_pot_bead(imode)` | 周期系のpotential virial | atomic unit |
| `charge`, `dipoler`, `hfcc` | 条件付き物性 | backend依存 |

PIMDの全potentialは `sum(pot_bead)/Nbead`、CLMDでは `pot_bead(1)` です。

### Gaussian (`Iforce=6`)

[`Set_Gaussian_MPI_tk.F90`](../Set_Gaussian_MPI_tk.F90) が `gauss.tmp` と `g0xrun_p` を確認し、ビーズ別Scratchへheaderを準備します。[`Force_Gaussian.F90`](../Force_Gaussian.F90) は現在座標を追記した `gauss.com` を作り、wrapperを呼び、`gauss.log` から次を検索します。

- `SCF Done` または `EUMP2`: energy
- Gaussian force table: force
- 保存フラグが有効な場合: Mulliken charge、dipole、isotropic Fermi contact

restart直後の最初のGaussian入力では `Guess=Read` を除去します。

### VASP (`Iforce=8`)

[`Set_VASP.F90`](../Set_VASP.F90) が `INCAR`, `POTCAR`, `KPOINTS`, `LATTICE` をビーズ別Scratchへ配置します。[`Force_VASP_MPI.F90`](../Force_VASP_MPI.F90) は座標を `POSCAR` へ追記して `vasp_run.sh` を呼び、`OUTCAR` からenergy、force、external pressureを読みます。圧力は `W_pot_bead` に変換されます。

この接続は内部処理の説明対象ですが、本マニュアルでは実行チュートリアルを提供しません。

### MACE (`Iforce=25`)

MACE接続は `_MACE_` buildで有効です。[`Force_MACE.F90`](../Force_MACE.F90) が環境変数を解決し、埋込みPython interfaceを初期化します。各評価では元素記号を原子番号へ変換し、座標Å、cell Å、PBC、device、model pathをPython helperへ渡します。返されたeV、eV/Å、stressをPIMD8のatomic unitへ変換します。

Python interfaceはsimulation終了時に一度だけfinalizeします。外部プロセスを毎step起動するGaussian/VASP接続とはデータフローが異なります。

## 6. energy、pressure、出力

PIMDのenergy集計は [`Ham_Temp.F90`](../Ham_Temp.F90)、CLMDは `Ham_Temp_Classical` です。ここで物理運動energy、ring-polymer spring energy、thermostat energyを集計し、周期系なら `calc_internal_pressure` を呼びます。

出力は次の3層です。

- `print_ini`: 新規計算で出力ファイルのheaderを作る
- `print_result`: `out_step` ごとの座標、物性、圧力、拘束値
- `print_ham`: `ham.dat` と100 stepごとの `std.out`

列定義は[出力ファイル](outputs.md)を参照してください。

## 7. restartと実行中停止

[`Restart.F90`](../Restart.F90) は、normal-mode座標・運動量・力とthermostat状態を保存します。更新中の破損を避けるため一旦 `restart.tmp` を作り、前世代を `restart1.dat` へ移してから置換します。

[`exit.F90`](../exit.F90) は10 stepごとに `input.inp` の `$Lexit` を読み直します。`.T.` なら共通abort処理を呼びます。通常の完走とは終了メッセージが異なるため、restartファイルを確認してから再開します。

## 8. 実装上の注意

- 入力parserは未知タグを無視します。新規タグ追加時は読込、既定値、MPI broadcast、入力表示、利用箇所を一組で更新します。
- PIMDとCLMDには別のbackend dispatchがあります。backend追加時は両方で必要かを明示します。
- backendの単位変換は接続ごとに行います。共有配列へ格納する時点の単位を維持します。
- MACEの `device` は現在scalar broadcast対象ではありません。MPI利用時にrank間で設定を揃える変更を行うまでは、MACEチュートリアルはserialを基準とします。
- ビーズ平均距離拘束の `cons_val` と `cons_strength` は現在MPI broadcastされません。拘束はserialを基準とします。
- `pressure.dat` は観測量出力であり、通常PIMD/CLMD driverにセル時間発展はありません。