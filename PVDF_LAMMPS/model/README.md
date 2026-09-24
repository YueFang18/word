# 模型与力场接入（当前没有给用户结构自动赋参）

当前只有分子结构、尚无确定力场，因此生产参数文件会主动退出，防止把未知参数算出的数当成 PVDF 性能。这里不提供虚构的 PVDF/共聚单体参数。

每种序列、每种物理状态需要：

1. `system.data`：`units real`、`atom_style full` 对应 data；完整 Masses、Atoms、Bonds、Angles、Dihedrals、必要的 Impropers。正交初始盒。Atoms 列为 `id mol type q x y z [ix iy iz]`。**image flags 保证整条链/每个离子连续展开**，不能对跨界长链随意清零。
2. `styles.in`：与实际参数完全匹配的 pair/bond/angle/dihedral/improper styles、special_bonds、mixing rule、kspace。两个 example 只给函数形式示例。
3. `coeffs.in`：所有 bonded/nonbonded/cross coefficients。代码 `read_data ... nocoeff` 后从此文件载入，避免来源混杂。若导出器把系数写在 data，可改用你自己的加载方案，但必须同步验证全部分支。CLASS2 的交叉项不可遗漏。
4. `groups.in`：按真实原子/分子映射选择组；看 groups.example.in。不同模型 type 编号可以不同，但化学含义与参数来源必须一致。
5. `hooks.in`：可选 SHAKE 等约束；**禁止**在此额外积分、加恒温/恒压或固定原子；各任务自己定义积分器。
6. `type_to_element.json`：例如 `{"1":"C","2":"H","3":"F"}`，供 QM 输出。实际包含所有选中原子类型。
7. `chain_map.json`：例如 `{"chains":[{"molecule_id":1,"backbone_ids":[1,4,7,10]}]}`。主链 ID 必须按连接顺序排列，不是数值排序；branched/network 架构须另定义每条主路径，本包链统计和 PPA 接口默认线性链。

## 按任务所需 groups

| 任务 | 必需组 |
|---|---|
| bulk equilibrate / oxidation | 无额外强制组 |
| mechanics / rheology | polymer、backbone；polymer 含完整链，backbone 仅主链骨架原子 |
| adhesion | polymer、substrate、anchor；pull 另需 grip；umbrella 另需完整分子 probe |
| cross_scan | probe、substrate（无共价跨组键，优先无溶剂小模型） |
| rheology adsorption | polymer、backbone、particles |
| transport | Li、anion（完整阴离子）、donor_poly、donor_solvent、donor_anion；及四个 RDF type/type-range 变量 |

不存在的非必需组可设 empty。输运 RDF 示例一次一个 type 或连续 `a*b` 类型范围；非连续类型用多对 RDF 或修改输出。所有配位组必须按化学归属定义，不能把共聚单体中的 F 与阴离子 F 混为同一种配位。

## 你现在应如何选力场

- 把每个共聚单体、端基、连接方式、头头/尾尾缺陷、支化/交联列出。检查二元和三元连接环境的参数覆盖，不能只看各孤立单体能否赋参。
- 可先比较兼容的 OPLS-AA 方案与 PCFF/COMPASS/专用 PVDF 参数。对主要未知连接处计算带环境考量的电荷与扭转能扫描；验证相对构象能、链尺寸、密度和溶剂化，再确定一套生产方案。
- 同一个力场函数形式不意味着所有参数可以直接混用；特别注意 1–4 缩放、LJ 9–6 与 12–6、混合规则和电荷标定。
- 用相同化学原子类型/电荷规则跨序列建模。不要为每个候选拟合完全不同精度的势后直接排名。
- 表面参数与交叉作用另做验证；`cross_scan.in` 可提供固定构型吸附曲线和法向力，与同构型 DFT 对比。
- 不同参数来源下出现的性能排名差异，本身就是应报告的模型不确定性。

只有物理模型、力场参数和充分采样都就绪，模板输出才可用于结构–性能比较。
