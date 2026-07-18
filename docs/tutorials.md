# PIMD8 チュートリアル

この章では、既存のテスト資産を基にGaussianとMACEの最小計算を準備します。外部プログラムと `pimd.exe` は利用できる状態であるとします。インストール方法は扱いません。

既存ファイルには作成者のローカル絶対パスが含まれます。実行前に、この章で示す箇所を自分の環境に合わせて置き換えてください。テスト資産そのものを変更したくない場合は、ディレクトリ全体を作業場所へコピーして使います。

## 1. GaussianによるCLMD

### 1.1 使用する既存例

実際に `Iforce=6` で計算した保存済み例は [`Tests/Gaussian/CL`](../Tests/Gaussian/CL/) です。水1分子を300 K、10 step、`dt=0.2 fs` で計算します。

注意: [`Tests/Gaussian/QM/input.inp`](../Tests/Gaussian/QM/input.inp) はディレクトリ名に `Gaussian/QM` とありますが、現在の入力値は `Iforce=31`（SPC/F2）です。Gaussian-PIMDの完成済みテストとしては使用しません。

### 1.2 作業ディレクトリの準備

```bash
cp -R /path/to/PIMD8/Tests/Gaussian/CL ./gaussian-clmd
cd gaussian-clmd
mkdir -p Result Scr
```

計算ディレクトリには次のファイルが必要です。

| ファイル | 用途 |
|---|---|
| `input.inp` | PIMD8の入力 |
| `gauss.tmp` | Gaussian入力のヘッダー。route section、title、電荷・多重度までを書く |
| `g0xrun_p` | Gaussianを起動するラッパー |

既存ディレクトリには `gauss.tmp` が保存されていないため、例えば次の内容で作成します。

```text
# B3LYP/6-31G force SCF=XQC

PIMD8 force calculation

0 1
```

`g0xrun_p` のGaussian実行ファイルを自分の環境に合わせます。スクリプトには3引数が渡されます。

```bash
#!/bin/bash
export GAUSS_SCRDIR="$3"
/path/to/g16 < "$1" > "$2"
```

第1引数は生成された `gauss.com`、第2引数は `gauss.log`、第3引数はビーズのScratchディレクトリです。

### 1.3 入力の確認

既存の [`input.inp`](../Tests/Gaussian/CL/input.inp) で、少なくとも次を確認します。

```text
$Nbead
1
$Isimulation
10
$Iforce
6
$Lrestart
.F.
$address_result
./Result
$address_scr
./Scr
$version
g16
```

`Lsave_charge=.T.` と `Lsave_dipole=.T.` のため、この例では `charge.dat` と `dipole.dat` も出力されます。Gaussian 09を使う場合は `$version` を `g09` に変更します。

### 1.4 実行

保存済み `run_pimd.sh` には古いローカルパス表記があります。テスト資産を変更せずに実行する場合は、計算ディレクトリから実行ファイルを直接指定します。

```bash
/path/to/PIMD8/pimd.exe
```

### 1.5 成功の確認

```bash
tail -n 8 std.out
ls Result/ham.dat Result/coor.xyz Result/restart.dat
```

`std.out` の入力確認部分で次を確認します。

```text
++++ Simulation type            CLMD
++++ Number of Beads               1
++++ Flag for Force Calc           6
```

末尾に `Simulation End!` があり、`Result/ham.dat` の最終stepが `10` なら、この例は最後まで進んでいます。保存済みの期待結果は [`Tests/Gaussian/CL/std.out`](../Tests/Gaussian/CL/std.out) と [`Tests/Gaussian/CL/Result/ham.dat`](../Tests/Gaussian/CL/Result/ham.dat) で確認できます。数値の完全一致ではなく、同じ列が有限値で出力されることを確認してください。

### 1.6 同じ系をPIMDにする

Gaussian-PIMDでは各stepに全ビーズのGaussian計算が必要です。まず短い計算で確認します。CLMD例のコピーで次の値を変更します。

```text
$Nbead
4
$Nstep
2
$Isimulation
0
```

`Iforce=6`、`Ncent=3`、`Nref=5` はそのまま使用できます。実行すると `Scr/00001/` から `Scr/00004/` が作られ、各ディレクトリに `gauss.com` と `gauss.log` が生成されます。`Result/coor.xyz` の1フレームは `Natom*Nbead` 行、この例では `3*4=12` 原子行になります。

## 2. MACEによるPIMD

### 2.1 使用する既存例

既存例は [`Tests/MACE/Input`](../Tests/MACE/Input/) です。水1分子、16ビーズ、300 K、10 stepのPIMDです。保存済みの [`std.out`](../Tests/MACE/Input/std.out) にはMACE modelとdeviceが記録されています。

### 2.2 作業ディレクトリの準備

```bash
cp -R /path/to/PIMD8/Tests/MACE/Input ./mace-pimd
cd mace-pimd
mkdir -p Result Scr
```

MACE接続では、PIMD8同梱のPython helperディレクトリとモデルファイルを環境変数で指定します。

```bash
export MACE_PYTHON_DIR=/path/to/PIMD8/MACE
export MACE_MODEL=/path/to/model.model
```

`MACE_PYTHON_DIR` には `mace_function.py` が存在する必要があります。既存の `run_pimd.sh` には作成者の絶対パスと `../medium.MACE.model` が書かれているため、そのまま使わず、上記環境変数を設定して実行ファイルを直接起動します。

### 2.3 入力の確認

```text
$Nbead
16
$Nstep
10
$dt
0.1
$Isimulation
0
$Iforce
25
$device
cpu
$Lperiodic
.T.
$lattice
17.0  0.0  0.0
0.0  17.0  0.0
0.0  0.0  17.0
```

既存入力に `$device` がなければ既定値 `cpu` が使われます。GPUを使用する設定では `cuda` を指定できますが、このチュートリアルは `cpu` を基準にします。

### 2.4 実行と確認

```bash
/path/to/PIMD8/pimd.exe
tail -n 12 std.out
head -n 6 Result/ham.dat
head -n 20 Result/coor.xyz
```

`std.out` で次を確認します。

```text
++++ Simulation type            PIMD
++++ Number of Beads              16
++++ Flag for Force Calc          25
++++ MACE force field +++++
++++ MACE model      /path/to/model.model
++++ MACE device     cpu
```

末尾に `Simulation End!` があり、`Result/ham.dat` のstep 0から10までが出力されれば完了です。周期系なので `Result/pressure.dat` も生成されます。

### 2.5 MACEでCLMDを行う

同じ入力のコピーで次を変更します。

```text
$Nbead
1
$Isimulation
10
```

他の設定を変えずに起動できます。CLMDでは `ham.dat` が6列、PIMDでは9列になります。周期系の `pressure.dat` と `coor.xyz` はどちらでも出力されます。

## 3. restart計算

正常終了または `out_step` ごとの出力時に `Result/restart.dat` が更新されます。計算を続けるには、同じ計算ディレクトリで次のようにします。

1. `Nstep` を現在の最終stepより大きくする。
2. `$Lrestart` を `.T.` にする。
3. `Natom`、`Nbead`、`Ndim`、`Ncent`、`Nnhc` をrestart作成時から変えない。
4. 同じ `address_result` を使って実行する。

`restart.dat` の先頭整数が再開stepです。既存の `run_pimd.sh` は `Lrestart=.F.` の場合にResultを削除する処理を持つため、restart時に誤って `.F.` のままにしないでください。

## 4. よくある停止原因

| 症状 | 確認箇所 |
|---|---|
| `There is no input file of "input.inp"` | `input.inp` のあるディレクトリから起動したか |
| `There is no "$end parameter"` | パラメータ部末尾に `$end parameter` があるか |
| `There is no "$Coords"` | `$Coords` と `$end Coords` があるか |
| Gaussianの入力ファイル不足 | `gauss.tmp` と `g0xrun_p` が計算ディレクトリにあるか |
| Gaussian出力の検索エラー | route sectionに `force` があるか、`gauss.log` が正常終了しているか |
| `MACE_PYTHON_DIR is not set` | helperディレクトリの環境変数を設定したか |
| `MACE_MODEL is not set` | モデルファイルの環境変数を設定したか |
| `Invalid $lattice for MACE periodic calculation` | `Lperiodic=.T.` のとき3本の格子ベクトルが有効か |
| restart読込で停止 | restart作成時と配列サイズを決める入力が一致しているか |