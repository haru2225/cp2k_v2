# cp2k_v2: 粘土エッジの原子電荷を、DZVPで最適化した構造の上でCP2Kで計算する

ピロフィライト(粘土)のエッジ付きリボンの局所環境ごとの原子電荷(Hirshfeld)を、CP2K(PBE、GTH、DZVP-MOLOPT-SR-GTH、400 Ry)で作る。
SOAP → 電荷、のちに SOAP → χ → QEq モデルの教師データが目的。

## v1(`haru2225/cp2k`)からの変更点: 最適化の基底を SZV → DZVP
v1は構造最適化を SZV-MOLOPT-SR-GTH(最小基底、300 Ry)で行い、電荷だけ DZVP で計算していた。
その最適化構造は結合が長すぎた(バルクの平均: Si–O 1.729 Å、O–H 1.02〜1.05 Å)。同じバルクを DZVP で最適化し直すと
Si–O 1.642 Å、Al–O 1.921 Å、O–H 0.969 Å と、PBEで期待される値に近くなった(典型値との比較は私の知識による。出典は未確認)。
結合長は電荷に強く効くので、v1の電荷ラベルは不正確な構造の上で計算されていた。**v2は最適化も DZVP/400 Ry で行う**。v1の結果は同梱しない。
`warmstart/` に、v1のSZV緩和構造(5構造)を入れてあり、DZVP最適化の出発点にする(収束が早いはず)。`rib_x_o00_al`、`rib_y_o37_si` は元の構造から。

## 使い方(スパコン)
```bash
git clone https://github.com/haru2225/cp2k_v2.git && cd cp2k_v2
qsub run_cp2k.pbs
```
これだけ。**構造を1つずつ順番に**処理する(ジョブ配列 `-J` は使わない): `bulk` → `rib_y_o00_si` → `rib_y_o00_al` → `rib_x_o00_si` → `rib_x_o00_al` → `rib_y_o37_al` → `rib_y_o37_si`
(最初の3つは、湿った系 `cp2k-water` の土台なので先に)。各構造で「DZVP構造最適化 → 電荷計算 → `runs/<名前>/charges.dat`」。
- 完了済み(`charges.dat`あり)はスキップ。**walltime(既定24時間)で途中で切れても、同じ `qsub run_cp2k.pbs` で続きから進む**(最適化は最後の構造から再開)。全部終わるまで何回かの再投入が要る見込み。
- 1構造だけ: `qsub -v NAME=rib_x_o00_al run_cp2k.pbs`。構造ごとに別ジョブにしたい場合は `bash submit_all.sh`。
- SCFが収束しない場合は自動でリトライ(最適化3回、1点計算2回): 1回目は速いOT、失敗したら直前の構造から、対角化 + Broyden混合 + 300 K Fermi–Dirac smearing の頑健な設定で続ける(smearingの電荷への影響は未検証)。
- `#PBS` の既定: `-q sc16`、`select=2:ncpus=16:mpiprocs=16`(2ノード)、`walltime=24:00:00`(仮置き)。コマンドラインで上書き可: `qsub -l walltime=48:00:00 run_cp2k.pbs`。
- 環境は自動: センターのサンプル(`/home/center/app/CP2K/cp2k_20251.sh`)の `module load` / `source` / `export PATH|LD_LIBRARY_PATH...` を再生、
  CP2K実行ファイルと `BASIS_MOLOPT` / `GTH_POTENTIALS` を探す(`-v CP2K_ENV_SCRIPT=...`、`CP2K_EXE=...`、`CP2K_DATA_DIR=...` で指定可)。計算ノードにpythonは不要。
- ログ: ジョブログの先頭に使った環境が出る。CP2Kの出力は `runs/<名前>/{opt_a*,sp_a*}.{out,stdout}`、失敗時はジョブログに末尾が出る。

計算時間の見込み(未実測、v1のSZVでは最適化に約2時間、電荷計算に約15分/16コア): DZVPは1ステップあたり数倍重いが、warm startで数十ステップのはず。1構造あたり数時間〜10時間程度か。

## 並列度(1構造あたりのノード数)
1つのジョブの中で、構造(系)は**1つずつ順番に**進める。ただし各構造では、CP2Kを**複数ノードでMPI並列**に走らせる。
MPIランク数は `$PBS_NODEFILE` の行数(全ノード分)から自動で決まり、複数ノードのときはMPIの種類(Open MPI / Intel MPI・MPICH)を見て `--hostfile` か `-f` を付ける
(`qsub -v MPI_OPTS="..."` で上書き、`MPI_OPTS=none` で無し)。
```bash
qsub run_cp2k.pbs            # 既定: 2ノード x 16 = 32ランク
bash go.sh 4                 # 4ノード x 16 = 64ランク(速いが、要求が大きいぶん待ち時間が伸びる)
bash go.sh 4 48:00:00        # walltimeも指定。以降の引数は qsub にそのまま渡す
```
`sc16` を16コア/ノードと仮定している。違うなら `qsub -l select=<N>:ncpus=<コア数>:mpiprocs=<コア数> run_cp2k.pbs`。
90〜150原子のCP2Kは64ランク前後までは伸びるはず(それ以上は効率が落ちる。未実測)。

## 次の段階: 水・イオン入りの系(リボンの最適化が終わってから)
湿った系は、DZVPで緩和したリボンを土台に作る必要がある。
```bash
# 1) 乾いた7構造が終わったら、結果を手元に持ってくる(runs/*/{relaxed.xyz,sp.out,opt_a*.out,charges.dat})
python3 collect_relaxed.py         # runs/<rib>/relaxed.xyz -> structures_relaxed/
python3 build_wet.py               # 要 numpy + LAMMPSのpythonモジュール。structures_wet/ に12系
python3 make_cp2k.py               # runs_wet/ の入力を作る
git add -A && git commit -m wet && git push
# 2) 湿った系の計算は cp2k-water で(独立しているので、この結果を待たなくてよい)
(cp2k-waterは独立しているので、このステップは不要。必要なら)
```
`build_wet.py` は、緩和したリボンのエッジ間の真空部にSPC水(約0.85 g/cm³、12〜23分子)を置き、リボンを固定した古典MD(ClayFF電荷、300 K、12 ps)でなじませ、6 ps・12 psの2スナップショットを出す。
`w` は水のみ(中性)、`wNa` は内部のAlを1つMgに置換してNa⁺を1つ置く(層電荷 −1、元のSi8Al3.5Mg0.5モデルと同密度)。Mgサイトは再緩和しない。

## 構造
- `structures/`(元の構造)と `runs/<名前>/`(入力): `bulk`(3D周期)、y法線リボン4(`rib_y_o00_{si,al}`、`rib_y_o37_{si,al}`)、x法線リボン2(`rib_x_o00_{si,al}`)。
  `o00`/`o37` は切る位置、`si`/`al` はプロトンの配分(Si–OH優先/Al–OH₂優先。カオリナイトのエッジの安定形とされる Si–OH と Al(OH)(OH₂) に対応)。形式電荷で中性。全方向周期セル + 真空。

## スクリプト
`run_cp2k.pbs`(ジョブ本体。構造を1つずつ)、`go.sh`(ノード数を指定して投入)、`submit_all.sh`(構造ごとに別ジョブ)、`make_cp2k.py`(入力の再生成)、`collect_relaxed.py`、`build_wet.py`、`build_edges.py`(`C2000.gro`、ClayCode、MIT、からリボンを作る)、
`hirshfeld_to_charges.py`、`analyze_charges.py`(SOAPモデルの検証: 要 numpy, ase, dscribe, scikit-learn)、`analyze_wet.py`(水によるリボン原子の電荷変化)。

## 確認の範囲
- ローカルのCP2K 2026.2で、バルクを DZVP で最適化(7ステップまで構造が正常に動くことを確認)。
- `run_cp2k.pbs` の流れ(最適化 → 電荷計算 → `charges.dat`、リトライ、再開)と `submit_all.sh` は偽のCP2K/qsubで確認。**センターの実機では未実行**。
- Hirshfeld電荷はClayFFの約1/3.8のスケール。MDに使うにはスケール換算が必要。
