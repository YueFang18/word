# PVDF 多结构/多元序列：LAMMPS 计算模板包

**当前状态：计算协议和后处理已整理；你的真实分子尚未接入、力场尚未确定。** 默认生产参数文件会退出。请先读 `model/README.md` 完成模型/力场接入，不能把此包当作已赋参的PVDF模型。

- `METHODS.zh.md`：五类目标的具体指标、公式、参数、收敛判据与序列对照设计。
- `lammps.in`：统一入口。
- `00_equilibrate/`：bulk NPT密度与参考构型。
- `01_adhesion/`：界面相互作用、法向/剪切脱粘、umbrella、交叉参数扫描。
- `02_mechanics/`：结构、单轴拉伸、剪切和完整6×6有限温度弹性矩阵的分量任务。
- `03_rheology/`：溶液结构、颗粒吸附、Green–Kubo和SLLOD流变。
- `04_transport/`：RDF/配位、Li MSD、完整轨迹和集体电荷矩。
- `05_oxidation/`：QM构型采样、外部DFT流程、带验证门槛的ReaxFF/NEB示例。
- `analysis/`：Numpy后处理：模量、分离功、WHAM、粘度、多原点扩散/电导、链统计与配位寿命、QM输入导出。
- `batch/`：序列描述符、每结构/工况/seed的批量命令与清单。
- `tests/`：数值测试和可在你安装LAMMPS后运行的玩具模型冒烟测试。玩具势不是PVDF力场。
- `VALIDATION.md`：本次实际验证范围。

依赖：LAMMPS（建议近期稳定版，所用styles/包需可用），Python>=3.10、Numpy>=2.0。ReaxFF和NEB另需REAXFF/REPLICA。ORCA是可选外部QM示例，不是LAMMPS依赖。

## 1. 第一轮接入

在包根目录运行所有命令。把真实文件放 `models/S01/`；用CLI覆盖不同模型，不必复制每个计算脚本：

```bash
lmp -in lammps.in -var task equilibrate \
  -var data models/S01/bulk.data \
  -var styles models/S01/styles.in -var coeffs models/S01/coeffs.in \
  -var groups models/S01/groups.in -var hooks model/hooks.in \
  -var T 300 -var seed 192845 -var out results/S01_eq
```

所有路径建议无空格/引号/`$`；LAMMPS变量替换与shell引号规则不同。`initvel yes`默认重新生成速度，`initvel no`要求data已经有适当Velocities；重新读取data不等于原轨迹无缝续跑。断点生产应另用read_restart并保持相关compute/fix状态，避免把重新初始化的MSD或相关函数拼接。

默认输出路径会被重用；**每次指定新的out目录**。批量生成器会在结果目录已存在时拒绝覆盖。固定键模板使用real/full；ReaxFF要求另一个charge-style无永久拓扑模型。

## 2. 五类目标运行示例

下面为了便于阅读省略 `-var styles ... -var coeffs ... -var groups ...`；实际必须提供，或配置model/下的生产接口。

```bash
# (1) 已建好的slab，不能直接拿bulk.data用
lmp -in lammps.in -var task adhesion -var mode equilibrium \
  -var data models/S01/interface.data -var out results/S01_adsorb
lmp -in lammps.in -var task adhesion -var mode pull_normal \
  -var data models/S01/interface.data -var pull_speed 1e-4 -var out results/S01_pull
# umbrella窗口示例：数据已经在该窗口附近预平衡；每窗不同初态和out
lmp -in lammps.in -var task adhesion -var mode umbrella \
  -var data models/S01/window_z10.data -var zcenter 10.0 -var Kspring 2.0 \
  -var out results/S01_u10

# (2) 用已完成bulk平衡且盒体积合适的模型
lmp -in lammps.in -var task mechanics -var mode tensile \
  -var data models/S01/bulk_reference.data -var strain_rate 1e-7 -var out results/S01_tension
lmp -in lammps.in -var task mechanics -var mode elastic \
  -var data models/S01/bulk_reference.data -var component 1 -var epsilon 0.005 \
  -var out results/S01_C1plus
# 对1..6分别运行+0.005/-0.005；同一参考构型与盒，不接力变形。

# (3) 溶液/真实浆料组成模型，先完成NPT平衡
lmp -in lammps.in -var task rheology -var mode gk \
  -var data models/S01/solution_reference.data -var out results/S01_gk
lmp -in lammps.in -var task rheology -var mode nemd \
  -var data models/S01/solution_reference.data -var shear_rate 1e-7 -var out results/S01_nemd

# (4) 完整电中性PVDF/盐/溶剂体系
lmp -in lammps.in -var task transport -var data models/S01/electrolyte_reference.data \
  -var out results/S01_transport

# (5) 仅采样，不会计算氧化电位
lmp -in lammps.in -var task oxidation -var data models/S01/electrolyte_reference.data \
  -var dump_every 2000 -var out results/S01_QM
```

## 3. 后处理示例

```bash
python analysis/postprocess.py traction results/S01_pull/traction.dat \
  --mode normal --output results/S01_pull/metrics.json
python analysis/postprocess.py modulus results/S01_tension/stress_strain.dat \
  --fit-min 0.002 --fit-max 0.01 --output results/S01_tension/metrics.json
python analysis/postprocess.py gk results/S01_gk/stress_raw.dat \
  --temperature 300 --maxlag-fs 20000 --blocks 5 \
  --curve results/S01_gk/eta_integral.dat --output results/S01_gk/metrics.json
python analysis/postprocess.py transport results/S01_transport/trajectory.dump \
  --li-ids results/S01_transport/Li_ids.dump --anion-ids results/S01_transport/anion_ids.dump \
  --dt-fs 0.5 --temperature 300 --maxlag-fs 250000 \
  --fit-min 100000 --fit-max 200000 \
  --curve results/S01_transport/msd_multi_origin.dat --output results/S01_transport/metrics.json
python analysis/association.py results/S01_transport/trajectory.dump \
  --li-ids results/S01_transport/Li_ids.dump --donor-types 18 --cutoff 3.0 --dt-fs 0.5 \
  --output results/S01_transport/association.csv --survival results/S01_transport/survival.dat
python analysis/structure.py results/S01_tension/trajectory.dump models/S01/chain_map.json \
  --output results/S01_tension/chains.csv --ppa-export results/S01_tension/backbones.jsonl
```

以上拟合区间、配位截断和相关窗口都只是命令示例，必须从你的轨迹判定，不能直接用其输出作结论。gk工具要求每统计块长于5个相关窗口，仍需你检查平台。transport工具拒绝非中性、可变盒和变电荷轨迹。`--discard`在表格分析中是数据行数，在transport中是轨迹帧数。

常见列：LAMMPS ave/time第一列都是step，第二列通常time_fs，后处理列索引从0开始。

| 文件 | 其余列 |
|---|---|
| bulk.dat | rho_g_cm3, V_A3, PE_kcal_mol, enthalpy_kcal_mol |
| interaction.dat | Ucross_kcal_mol, Wcross_J_m2, contact_atom_fraction, A_A2 |
| traction.dat | displacement_A, normal_MPa, shear_MPa, A_A2 |
| stress_strain.dat | ex, ey, ez, sigma_x_MPa, sigma_y_MPa, sigma_z_MPa |
| elastic_stress.dat | sigma_xx, yy, zz, yz, xz, xy (MPa) |
| stress_raw.dat | Pxy, Pxz, Pyz (atm), V_A3 |
| nemd.dat | Pxy_atm, eta_Pa_s, thermal_T_K |
| transport.dat | MSDx,y,z,total (Å²), V_A3, CN_poly, CN_solvent, CN_anion |
| charge_moment.dat | Mx, My, Mz (eÅ), V_A3 |
| window.dat | z_COM_minus_anchor_A, spring_energy_kcal_mol |

WHAM配置示例：

```json
[
  {"file":"results/S01_u10/window.dat","center_A":10.0,"K_kcal_mol_A2":2.0,"discard":100,"stride":10},
  {"file":"results/S01_u105/window.dat","center_A":10.5,"K_kcal_mol_A2":2.0,"discard":100,"stride":10}
]
```

`python analysis/postprocess.py wham windows.json --temperature 300 --bin-width 0.2`。实际需覆盖吸附盆地到脱附平台的完整窗口；stride按相关时间选择；需要独立窗口组/块bootstrap评估PMF误差，程序不伪造误差条。

完整弹性配置 `elastic_jobs.json`：键"1"到"6"，每项为 `{"plus":"results/C1plus/elastic_stress.dat","minus":"results/C1minus/elastic_stress.dat","epsilon":0.005}`。随后 `python analysis/postprocess.py elastic elastic_jobs.json`。输出原始C、对称C、反对称性和稳定性；输入的12个任务必须来自同一参考状态。

## 4. 多序列批量比较

```bash
# 编辑manifest中的序列、模型路径、物理状态、cases和seeds。
python batch/generate.py batch/study.example.json --output batch/generated
# 真正生产前开启模型路径检查：
python batch/generate.py my_study.json --output batch/production --require-models
bash batch/production/run_all.sh
```

生成：`sequence_features.csv`、每次任务的run_manifest.json（包括模型元数据/输入文件哈希）和run_all.sh。三种速度seed只是最初重复设计；独立packing、晶区形态和随机序列实现需要分别提供模型记录。输出目录唯一；不自动跨物理状态把上一步的data拿给下一步。

## 5. 验证与限制

```bash
python -m unittest discover -s tests -v
python tests/smoke_lammps.py --lmp /path/to/lmp --work-dir /tmp/pvdf_smoke
```

数值测试验证公式、单位、文件解析、WHAM和序列描述符；玩具模型冒烟测试仅验证LAMMPS脚本分支，不验证PVDF物理正确性。此交付环境没有LAMMPS，实际运行情况见VALIDATION.md。

**需要外部步骤的性质**：晶型/结晶度的参考结构标定、PPA/Z1+缠结、混溶自由能/χ、DFT电离能与氧化态电荷/反应能垒。本包提供相应采样/导出及方法说明，未把它们标成已自动算出的结果。
