# `input.inp` リファレンス

PIMD8は、タグ行の**直後の行**から値を読みます。タグは行頭から書き、大文字・小文字をソースと同じにしてください。例えば次の形式です。

```text
$Nbead
16
```

パラメータ部は `$end parameter` で閉じ、その後に `$Coords` ブロックを置きます。論理値は `.T.` または `.F.` を使います。

## 1. 最小入力の骨格

```text
$Ndim
3
$Nbead
1
$temperature
300.0
$Nstep
10
$dt
0.2
$Isimulation
10
$out_step
1
$Ncent
3
$Lrestart
.F.
$Iforce
6
$address_result
./Result
$address_scr
./Scr
$seed
101
$Langstrom
.T.
$Lrandom_coor
.F.
$Lperiodic
.F.
$Nref
5
$Nys
5
$Nnhc
4
$gamma
1.0
$Lexit
.F.
$end parameter

$Coords
O  0.0  0.0  0.0
H  1.0  0.0  0.0
H  0.0  1.0  0.0
$end Coords
```

## 2. 基本計算タグ

`既定値なし` の項目は、実用上 `input.inp` に必ず明示してください。

| タグ | 型・単位 | 既定値 | 意味・制約 |
|---|---|---:|---|
| `$The Number of Proccesors` | 整数 | なし | `run_pimd.sh` がMPI process数として読む補助値。Fortran本体は読まない。タグ名の綴りは既存のまま使用する |
| `$Ndim` | 整数 | `3` | 空間次元。`1` または `3` |
| `$Nbead` | 整数 | なし | ビーズ数。CLMDは `1`、PIMDは `2` 以上。既存PIMD例は偶数を使用 |
| `$temperature` | 実数、K | なし | 設定温度。正の値を指定する |
| `$Nstep` | 整数 | なし | 計算を終了するstep番号。restart時はrestart stepより大きくする |
| `$dt` | 実数、fs | なし | 外側の時間刻み。読込後に原子単位へ変換される |
| `$Isimulation` | 整数 | なし | `0`: PIMD、`1`: RPMD、`2`: CMD、`3`: PIHMC、`10`: CLMD。通常利用は `0` と `10` |
| `$out_step` | 整数、step | `1` | 軌跡、各種データ、restartの出力間隔。正の値が必要 |
| `$Ncent` | 整数 | なし | centroid thermostat。`0`: NVE、`1`: centroid全体への1本のNHC、`3`: 原子・方向ごとのNHC |
| `$Lrestart` | 論理 | `.F.` | `.T.` なら `address_result/restart.dat` を読む |
| `$Iforce` | 整数 | なし | 力計算backend。[backend一覧](#5-iforce-backend一覧)を参照 |
| `$address_result` | 文字列 | なし | 結果ディレクトリ。存在しなければ本体が作成する |
| `$address_scr` | 文字列 | なし | Scratchディレクトリ。存在しなければ本体が作成する |
| `$seed` | 整数 | なし | 乱数seed。初期速度などに使う |
| `$Langstrom` | 論理 | `.T.` | `.T.` なら `$Coords` の座標をÅとして読み、内部でbohrへ変換する |
| `$Lrandom_coor` | 論理 | `.F.` | PIMD初期化で `.T.` なら非centroid normal modeを自由粒子分布から乱数生成する。`.F.` なら全ビーズを同じ初期構造から始める |
| `$Lperiodic` | 論理 | `.F.` | 周期境界の有無。`.T.` なら `$lattice` を指定する |
| `$lattice` | 実数3行、Å | 単位行列 | 3本の格子ベクトルを各行に指定する。MACEへはこの行列がそのまま渡る |
| `$Lexit` | 論理 | タグ省略時は停止要求なし | 10 stepごとに本体が `input.inp` を読み直す。`.T.` に書き換えると共通停止処理を呼ぶ |

`Isimulation=0`, `1`, `2` で `Nbead=1` を指定すると、入力処理が `Isimulation=10` へ変更します。

## 3. 積分・thermostatタグ

| タグ | 型・単位 | 既定値 | 意味・制約 |
|---|---|---:|---|
| `$Nref` | 整数 | なし | PIMDの内側時間刻み数。`dt_ref=dt/Nref`。正の値を指定する。CLMDでは `dt_ref=dt` となり分割しない |
| `$Nys` | 整数 | なし | Nose–Hoover chain積分のYoshida–Suzuki分解数。実装値は `1`, `3`, `5` |
| `$Nnhc` | 整数 | なし | Nose–Hoover chainの長さ。`Ncent>0` で使用する |
| `$gamma` | 実数 | `1.0` | PIMD/CMDの仮想質量スケーリング係数。内部では二乗して使用する |
| `$freq1` | 実数、fs | `10.0` | centroid thermostat質量を決める基準周期。`omega_system=2π/freq1` として使用する |

通常のNVT例は `Ncent=3`, `Nref=5`, `Nys=5`, `Nnhc=4`, `gamma=1.0` を使用しています。`Ncent=0` ではcentroid thermostatを適用しません。

## 4. 保存タグ

| タグ | 型 | 既定値 | 生成物・注意 |
|---|---|---:|---|
| `$Lsave_force` | 論理 | `.F.` | `.T.` で `force.dat` |
| `$Lsave_energy` | 論理 | `.F.` | `.T.` でビーズ別ポテンシャルを `ene.dat` |
| `$Lsave_charge` | 論理 | `.F.` | `.T.` で `charge.dat`。値を提供するbackendでのみ有効 |
| `$Lsave_dipole` | 論理 | `.F.` | `.T.` で `dipole.dat`。値を提供するbackendでのみ有効 |
| `$Lsave_hfcc` | 論理 | `.F.` | `.T.` で `hfcc.dat`。Gaussian出力に対応sectionが必要 |

`ham.dat`、`coor.xyz`、`restart.dat` は保存タグに関係なく出力されます。`Lperiodic=.T.` では `pressure.dat`、`Icons>0` では `constraint.dat` も出力されます。

## 5. `Iforce` backend一覧

次は現在の力計算dispatchに接続されている値です。コンパイル条件付きbackendは、対応するプリプロセッサ設定を含む実行ファイルでのみ使えます。

| 値 | backend | PIMD | CLMD | 補足 |
|---:|---|:---:|:---:|---|
| `1` | MOPAC | ○ | ○ | 外部ファイル方式 |
| `6` | Gaussian | ○ | ○ | `gauss.tmp`, `g0xrun_p` が必要 |
| `8` | VASP | ○ | ○ | `INCAR`, `POTCAR`, `KPOINTS`, `LATTICE`, `vasp_run.sh` を使用 |
| `9` | SIESTA | ○ | — | PIMD側dispatchのみ |
| `11` | Morse model | ○ | ○ | 内蔵モデル |
| `15` | Harmonic model | ○ | ○ | 内蔵モデル |
| `16` | Double Morse | ○ | ○ | 内蔵モデル |
| `21` | Araidai NNP | ○ | ○ | 外部ファイル方式 |
| `22` | Matlantis旧接続 | ○ | ○ | サブプロセス方式 |
| `24` | LAMMPS | 条件付き | 条件付き | `_LAMMPS_` 有効時 |
| `25` | MACE | 条件付き | 条件付き | `_MACE_` 有効時。`MACE_PYTHON_DIR`, `MACE_MODEL` が必要 |
| `26` | ænet | ○ | — | PIMD側dispatchのみ |
| `31` | SPC/F2 water | ○ | ○ | 内蔵水モデル |
| `32` | LJ-PBC | ○ | ○ | 周期Lennard-Jonesモデル |

### Gaussian固有タグ

| タグ | 型 | 既定値 | 意味 |
|---|---|---:|---|
| `$version` | 文字列 | `g16` | Gaussian出力のMulliken charge見出しを選ぶ。`g16` または `g09` |

Gaussianは各ビーズについて `gauss.com` を生成し、`g0xrun_p` を呼び、`gauss.log` からenergyとforceを読みます。保存フラグに応じてcharge、dipole、HFCCも読みます。

### MACE固有タグ

| タグ | 型 | 既定値 | 意味 |
|---|---|---:|---|
| `$device` | 文字列 | `cpu` | MACE calculatorへ渡すdevice。通常は `cpu`、対応環境では `cuda` |

モデルとhelperは入力タグではなく環境変数で指定します。

```bash
export MACE_PYTHON_DIR=/path/to/PIMD8/MACE
export MACE_MODEL=/path/to/model.model
```

## 6. PIHMCタグ

通常のPIMD/CLMDでは既定値のままで、指定する必要はありません。

| タグ | 型 | 既定値 | 意味 |
|---|---|---:|---|
| `$Ndyn` | 整数 | `20` | `Isimulation=3` の1 trial trajectory内のMD step数。正の値が必要 |
| `$Ldual` | 論理 | `.F.` | `.T.` ならproposalとMetropolis判定で異なるpotentialを使用する |
| `$dual_Iforce` | 整数 | `0` | `Ldual=.T.` のtarget backend。正のdispatch値が必要 |

ここでは入力の意味だけを示します。通常計算のチュートリアルではPIHMCを扱いません。

## 7. ビーズ平均距離拘束

この拘束はPIMDのnormal-mode内力に加算されます。現在の実装ではMPI集約を行わないため、serial実行で使用してください。CLMDの時間積分経路には拘束力が接続されていません。

| タグ | 型・単位 | 既定値 | 意味 |
|---|---|---:|---|
| `$Icons` | 整数 | `0` | `0`: 無効、`1`: ビーズ平均原子間距離への調和拘束 |
| `$cons_atom1` | 整数 | なし | 第1原子の1始まりindex |
| `$cons_atom2` | 整数 | なし | 第2原子の1始まりindex |
| `$cons_atom3` | 整数 | なし | `Icons=1` では未使用 |
| `$cons_val` | 実数、Å | なし | ビーズ平均距離の目標値 |
| `$cons_strength` | 実数、Hartree/bohr² | なし | 調和拘束の力定数 |

拘束ポテンシャルは `0.5*cons_strength*(平均距離-cons_val)^2` です。`constraint.dat` へ平均距離、微分、拘束エネルギーを出力します。

## 8. 座標ブロック

```text
$Coords
O  0.0  0.0  0.0
H  1.0  0.0  0.0
H  0.0  1.0  0.0
$end Coords
```

各行は `元素記号 x y z` です。原子数は `$Coords` と `$end Coords` の間の行数から決まります。空行も原子行として数えられるため、ブロック内に空行やコメントを入れないでください。

`Langstrom=.T.` なら座標単位はÅ、`.F.` なら内部単位のbohrとして読みます。`D` と `mu` は質量判定後、力計算backendへ渡す元素ラベルとしては `H` に置換されます。

## 9. 入力時の注意

- タグ行の後ろにはコメントを書けますが、タグ自体は行頭から一致する必要があります。
- 未知のタグはエラーにならず無視されます。綴り間違いに注意してください。
- `Ndim`, `Nbead`, `Nstep`, `temperature`, `dt`, `Isimulation`, `Ncent`, `Iforce`, `seed`, `Nref`, `Nys`, `Nnhc`, `address_result`, `address_scr` は明示することを推奨します。
- `Nref`, `out_step`, `Nbead`, `temperature`, `dt` には正の値を使用してください。これらの全条件が入力時に検査されるわけではありません。
- `input.inp` テンプレートにある項目でも、現行実行経路で値が使われないものは、このリファレンスの利用可能タグに含めていません。

実装上の入力元は [`read_input.F90`](../read_input.F90)、既定値は [`Parameter.F90`](../Parameter.F90) です。