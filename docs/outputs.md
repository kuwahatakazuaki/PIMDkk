# 出力ファイル

出力先は `$address_result` で指定します。`std.out` と `std.err` はResultディレクトリではなく、実行時のカレントディレクトリに置かれます。

## 1. 常に確認するファイル

| ファイル | 内容 | 出力間隔 |
|---|---|---|
| `std.out` | 入力確認、100 stepごとのエネルギー概要、開始・終了時刻 | 開始時、100 stepごと、終了時 |
| `Result/ham.dat` | Hamiltonian、温度、エネルギー成分 | `out_step` ごと |
| `Result/coor.xyz` | 全ビーズの実座標、Å | `out_step` ごと |
| `Result/restart.dat` | 再開に必要な座標・運動量・力・thermostat状態 | `out_step` ごと |
| `Result/restart1.dat` | 1世代前のrestart backup | 2回目以降のrestart更新時 |

## 2. `std.out`

実行開始時にファイルを作り直し、次を記録します。

- 読み取ったsimulation種別、原子数、次元、ビーズ数、step数
- 温度、時間刻み、出力時間間隔
- thermostat、backend、乱数seed
- Result/Scratch、restart、周期境界、格子
- MACE使用時のmodelとdevice
- 拘束使用時の拘束パラメータ
- 原子ラベル、質量、初期座標

正常終了の目印は末尾の次の表示です。

```text
***********************
    Simulation End!
***********************
```

この表示がない場合は、実行を捕捉したshellの出力、`std.err`、backendのScratch出力も確認してください。

## 3. `ham.dat`

エネルギーの内部単位はHartree、温度はKです。

### CLMD（`Nbead=1`）

| 列 | 内容 |
|---:|---|
| 1 | step |
| 2 | Hamiltonian |
| 3 | 瞬間温度、K |
| 4 | potential energy |
| 5 | physical kinetic energy (`DKinetic`) |
| 6 | centroid thermostat energy (`EBath_Cent`) |

### PIMD（`Nbead>1`）

| 列 | 内容 |
|---:|---|
| 1 | step |
| 2 | Hamiltonian |
| 3 | 瞬間温度、K |
| 4 | potential energy |
| 5 | physical kinetic energy (`DKinetic`) |
| 6 | ring-polymer spring energy (`QKinetic`) |
| 7 | 非centroid modeのthermostat energy (`EBath`) |
| 8 | centroid thermostat energy (`EBath_Cent`) |
| 9 | virial kinetic-energy estimator (`E_Virial`) |

step 0も `out_step` の条件を満たすため、非restart計算では初期状態が記録されます。Hamiltonianの絶対値だけで良否を判断せず、時間変化に不連続、NaN、Inf、急激な発散がないかを確認します。

## 4. `coor.xyz`

extended XYZ形式で、座標単位はÅです。1フレームの原子数は `Natom*Nbead` です。同じ原子順序をビーズ1から `Nbead` まで繰り返します。

非周期系のコメント行は次の情報を持ちます。

```text
Properties=species:S:1:pos:R:3 pbc="F F F" step=0
```

周期系ではさらに `Lattice="..."` と `pbc="T T T"` を持ちます。格子の単位はÅです。

## 5. 条件付き出力

| ファイル | 生成条件 | 内容・単位 |
|---|---|---|
| `force.dat` | `Lsave_force=.T.` | step見出しの後に、ビーズ・原子順の3成分。内部の原子単位 |
| `ene.dat` | `Lsave_energy=.T.` | step見出しの後にビーズ別potential、Hartree |
| `charge.dat` | `Lsave_charge=.T.` | step、ビーズごとの原子charge。backend依存 |
| `dipole.dat` | `Lsave_dipole=.T.` | step、ビーズごとのx, y, z成分。backend出力の単位 |
| `hfcc.dat` | `Lsave_hfcc=.T.` | step、ビーズごとの原子HFCC。backend出力の単位 |
| `constraint.dat` | `Icons>0` | step、時刻fs、平均距離Å、`dV/dCV`、拘束energy |
| `pressure.dat` | `Lperiodic=.T.` | 圧力・体積・格子長・potential virial |
| `pihmc.out` | `Isimulation=3` | PIHMCのaccept/reject情報 |

保存フラグを有効にしても、backendが該当物性を設定しなければ意味のある値は得られません。Gaussianはcharge、dipole、HFCCの読取経路を持ちます。

## 6. `pressure.dat`

ヘッダーと列は次のとおりです。

| 列 | 内容 | 単位 |
|---:|---|---|
| 1 | step | — |
| 2 | time | fs |
| 3 | centroid-virial pressure (`P_cv`) | bar |
| 4 | primitive pressure (`P_prim`) | bar |
| 5 | instantaneous pressure (`P_inst`) | bar |
| 6 | volume | Å³ |
| 7–9 | 格子ベクトルの長さ `a`, `b`, `c` | Å |
| 10 | potential virial `W_pot` | atomic unit |

この出力は圧力の観測値です。現在の通常PIMD/CLMDでは、圧力に応じてセルを時間発展させるNPT動力学を意味しません。

## 7. restartファイル

`restart.dat` の先頭行はstepです。その後に、ビーズ・原子順でnormal-mode座標、質量スケールした運動量、力、必要に応じてthermostat変数が続きます。人手で編集するための形式ではありません。

更新は次の順で行われます。

1. 新しい状態を `restart.tmp` へ書く。
2. 既存の `restart.dat` を `restart1.dat` へ移す。
3. `restart.tmp` を `restart.dat` へ移す。

再開時は配列サイズとthermostat構成を作成時と一致させます。特に `Natom`, `Ndim`, `Nbead`, `Ncent`, `Nnhc` を変更しないでください。

実装上の出力元は [`print_ini.F90`](../print_ini.F90)、[`print_result.F90`](../print_result.F90)、[`print_ham.F90`](../print_ham.F90)、[`Restart.F90`](../Restart.F90) です。
