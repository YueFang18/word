# 实现与方法依据

检索/核对日期：2026-09-25。LAMMPS站点文档显示release 2Sep2026；请对照实际安装版本检查功能/包。

## LAMMPS官方文档

- [group/group：相互作用能与力、非键/多体限制](https://docs.lammps.org/compute_group_group.html)
- [spring：对完整group COM施加umbrella偏置](https://docs.lammps.org/fix_spring.html)
- [gyration/chunk：分链回转半径与张量](https://docs.lammps.org/compute_gyration_chunk.html)
- [dihedral/local：主链二面角及单位](https://docs.lammps.org/compute_dihedral_local.html)
- [弹性计算](https://docs.lammps.org/Howto_elastic.html)
- [Green–Kubo、非平衡粘度](https://docs.lammps.org/Howto_viscosity.html)
- [SLLOD](https://docs.lammps.org/fix_nvt_sllod.html)
- [deform与remap](https://docs.lammps.org/fix_deform.html)
- [MSD与COM扣除](https://docs.lammps.org/compute_msd.html)
- [RDF与配位积分](https://docs.lammps.org/compute_rdf.html)
- [property/atom与展开坐标](https://docs.lammps.org/compute_property_atom.html)
- [ReaxFF电荷平衡](https://docs.lammps.org/fix_qeq_reaxff.html)
- [ReaxFF物种分析](https://docs.lammps.org/fix_reaxff_species.html)
- [NEB](https://docs.lammps.org/neb.html)
- [GPU/KOKKOS加速包比较](https://docs.lammps.org/Speed_compare.html)

## 与本体系直接相关的原始研究

1. Lee et al. (2014), Molecular Dynamics Simulations of the Traction-Separation Response at the Interface between PVDF Binder and Graphite in the Electrode of Li-Ion Batteries. [DOI](https://doi.org/10.1149/2.0051409jes)
2. Lee (2016), Molecular Dynamics Study of the Separation Behavior at the Interface between PVDF Binder and Copper Current Collector. [DOI](https://doi.org/10.1155/2016/4253986)
3. Park et al. (2023), Crystallization behavior of polyvinylidene fluoride (PVDF) in NMP/DMF solvents: a molecular dynamics study. [DOI](https://doi.org/10.1039/D3RA00549F)
4. Investigation of molecular mechanisms of polyvinylidene fluoride under the effects of temperature, electric poling, and mechanical stretching using molecular dynamics simulations (2022). [DOI](https://doi.org/10.1016/j.polymer.2022.124691)
5. Vijayakumar et al. (2025), New Quantum Mechanics Based Force Field for Describing Dynamics and Piezoelectric Properties for Various Phases of Poly(vinylidene) Fluoride Polymer. [DOI](https://doi.org/10.1021/acs.jpcc.5c00409)
6. Gao et al. (2022), C/H/O/F/Al ReaxFF Force Field Development and Application to Study the Condensed-Phase Poly(vinylidene fluoride) and Reaction Mechanisms with Aluminum. [DOI](https://doi.org/10.1021/acs.jpcc.2c02043)
7. Yu et al. (2013), The Study on Prediction of Oxidative Decomposition Potential by Comparison between Simulation and Electrochemical Methods to Develop the Binder for High-voltage Lithium-ion Batteries. [DOI](https://doi.org/10.5229/JKES.2013.16.3.177)

## 可迁移方法（不自动等同于PVDF电极验证）

- Fadel et al. (2019), Role of solvent-anion charge transfer in oxidative degradation of battery electrolytes. [原文](https://www.nature.com/articles/s41467-019-11317-3)
- The Z1+ package: Shortest multiple disconnected path for the analysis of entanglements in macromolecular systems (2023). [DOI](https://doi.org/10.1016/j.cpc.2022.108567)
- [ORCA电荷/多重度/外部XYZ输入](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html)

脚本是依据这些方法新整理的模板，不是上述论文作者代码的逐字复现。不能引用论文的验证结果来代替你使用本模板的验证。
