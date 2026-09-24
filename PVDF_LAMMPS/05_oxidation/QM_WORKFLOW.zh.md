# 氧化反应：MD 采样之后怎样得到目标参数

## 1. 先定义电子转移对象

研究本征序列效应：使用相同端基、足够长的低聚物，比较相同组成但不同排列；做聚合度收敛。研究服役稳定性：加入代表性溶剂、阴离子、Li+，必要时建立实际荷电态/表面态的周期电极模型。上述两种结果不能混作同一标签。

`lammps.in -var task oxidation` 仅产生经典构型。`q` 是力场部分电荷，`pe` 是势能；二者都不是氧化态电荷或电离能。

## 2. 垂直电离能 VIP

针对每一个独立构型 R，在完全相同几何与环境下分别计算原电荷 Q 和 Q+1：

`VIP_eV = [E_DFT(Q+1,R) - E_DFT(Q,R)]_Hartree * 27.211386246`

多个构型报告均值、标准差和分布，按链端/支化点/不同共聚二元或三元片段分层。跨序列保持同样的低聚物长度策略、端基、溶剂化方法、基组和泛函；较低尾部可能揭示易氧化环境，但不能直接用分布最小值作为实验氧化电位。

导出示例（替换 molecule-ID、center-ID 和电荷/多重度）：

```bash
python analysis/export_qm.py results/QM/qm_sampling.dump \
  --type-map models/S01/type_to_element.json --frame 20 \
  --molecules 1 51 62 --center-id 12 \
  --charge 0 --mult-initial 1 --mult-oxidized 2 --output results/QM/frame20
```

生成 `snapshot.xyz` 与两个 ORCA 单点输入。只有知道具体体系后才能选择正确整数电荷和自旋多重度；0/1 和 +1/2 只适用于相应闭壳层初态/双重态氧化态。默认 `wB97X-D3 def2-TZVP TightSCF` 是气相示例，不构成适用于所有化学环境的推荐，更不是自动的溶液氧化电位。必须检查 SCF 收敛、自旋污染、空穴定位及方法/基组敏感性。

导出器只保留完整分子、不切断聚合物，不做自动加氢封端。很长链的 QM 计算可能不可承受，应另建经验证的低聚物片段并明确端基效应。完整分子的 image flags 必须正确。近表面或延伸网络用周期 DFT 或专门 QM/MM 建模；不要用本脚本把基底裁成任意簇。

## 3. 绝热电离能与氧化自由能

分别优化两个电荷态，`AIP = E(Q+1,R_ox_opt)-E(Q,R_init_opt)`；加入匹配的环境自由能/热修正或适当自由能取样才能得到 ΔG。VIP、AIP 和 ΔG 均不是同一个量。溶液/界面环境可用显式溶剂结合嵌入、经验证的连续介质或周期方法；注意垂直跃迁与绝热过程的溶剂响应不同。

欲报告相对 Li/Li+ 的电位，明确完整电子转移反应、电子化学势和参考电极；建议用同一计算水平的已知氧化还原对标定。不能直接把 `-HOMO` 或未标定的 VIP 当作 V vs Li/Li+。

## 4. 氧化态电荷定位

在相同几何输出初态/氧化态电子密度和自旋密度：

- 差分密度 `rho_initial - rho_oxidized`，其全空间积分应对应失去 1 个电子；
- Hirshfeld/Bader 等一致分区方法，计算每单体/片段的空穴权重；
- 比较端基、VDF–共聚单体连接处、官能团的空穴/自旋富集程度。

差分密度可能局部为负，原子电荷变化也不一定都为正；不能随意截断后称作概率。若计算参与率/局域化指数，必须说明使用的非负空穴密度或轨道投影定义。仅凭固定经典电荷不能做本分析。

## 5. 反应能与能垒

选择具体反应，例如氧化后脱氢、脱氟、C–C 断裂、与表面氧/溶剂反应。固定原子数、总电荷、自旋态及环境的可比性：

- `DeltaE_rxn = E_products - E_reactants`；自由能同理为 `DeltaG_rxn`；
- `DeltaE_dagger = E_TS - E_reactants`；有限温度关注 `DeltaG_dagger`；
- DFT NEB 或过渡态搜索后，检查唯一相关虚频并沿反应坐标验证连通两个端点；
- 若有不同电荷/自旋面交叉，需要显式电子转移或非绝热分析，单条同电荷 NEB 不足以解决。

`neb_reaxff.in` 是已验证反应势下的同电荷经典势能 NEB 示例。运行需要 LAMMPS REPLICA/REAXFF 包、多个 replicas、同一 atom-ID 的反应两端点。`finalcoords` 文件按 LAMMPS NEB 格式：第一行 N，随后 `id x y z`；初末坐标/周期像需构造合理，不能让插值穿过其他分子。

`reaxff.in` 可记录键级和产物，但必须已有覆盖全部元素、交叉反应和相关环境的参数集。不要拼接独立来源参数，也不要认为 C/H/O/F/Al 参数自动覆盖 Li/P/Ni/Co/Mn。它没有自动的电位控制和量子电子转移；高温事件率不直接等于室温高电压氧化速率。

## 6. 建议最终数据字段

`sequence_id, motif, environment, snapshot_id, oligomer_DP, terminal_groups, Q, multiplicity, QM_method, basis, solvation, VIP_eV, AIP_eV, DeltaG_ox_eV, reference_electrode, DeltaE_rxn_eV, DeltaG_dagger_eV, hole_localization_definition, uncertainty`

对没有真正算出的字段留空，不用经典 MD 描述符代填。

## 方法依据

- [Yu 等，2013，PVDF 等粘结剂的量子化学/LSV 比较](https://doi.org/10.5229/JKES.2013.16.3.177)
- [Fadel 等，2019，MD 构型的电离能分布与环境耦合](https://www.nature.com/articles/s41467-019-11317-3)
- [Gao 等，2022，C/H/O/F/Al ReaxFF](https://doi.org/10.1021/acs.jpcc.2c02043)
- [ORCA 坐标、电荷和多重度输入](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html)
- [LAMMPS NEB](https://docs.lammps.org/neb.html)
