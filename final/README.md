# Three-Phase — Final Deliverables 索引

## 数据文件 (`clean_data/`)

每个体系 4 个 CSV：

| 文件 | 内容 | Origin 用途 |
|---|---|---|
| `binodal_{体系}_concentrated.csv` | 浓相（给体富集） | 三元图的"左臂" |
| `binodal_{体系}_dilute.csv` | 稀相（受体富集） | 三元图的"右臂" |
| `spinodal_{体系}.csv` | 旋节线 | 三元图虚线 |
| `critical_{体系}.csv` | 临界点 | 三元图标记点 |

**六体系**：PM6_Tol / PM6_OXy / D18_Tol / D18_OXy / PM6_CF / D18_CF

**Origin 三元图列映射**：X = phi1_L8Bo（受体）、Y = phi3_Donor（给体）、Z = phi2_Solvent（溶剂）

---

## 图片文件

### MATLAB 生成（`output/figures/`）

| 图片 | 说明 |
|---|---|
| `ternary_four_systems.png` | 四体系综合三元相图（PM6/D18 × Tol/OXy） |
| `ternary_grid_2x2.png` | 2×2 分面三元相图 |
| `ratio_grid_2x2.png` | 2×2 分面 Ratio-Solvent 图 |
| `ratio_donor_compare.png` | 同溶剂下给体对比（PM6 vs D18 / Tol） |
| `ratio_solvent_compare.png` | 同给体下溶剂对比（Tol vs OXy / PM6） |
| `*_HR.png` | 300 DPI 高清版本 |

### Python 生成（`final/`）

| 图片 | 说明 |
|---|---|
| `fig6_drying_trajectory.png` | 干燥轨迹概念示意图 |
| `fig6_underlying_data.csv` | 概念图对应的临界点数据 |

---

## 物理解读

| 文档 | 说明 |
|---|---|
| `reports/相图指导溶剂选择_完整解读.md` | 为什么 CF 最好、OXy 绿溶剂最优 |
| `final/Fig6_TieLine_解释文档.md` | 临界点和结线斜率的通俗解释 |

---

## 重新生成

- **MATLAB 图**：运行 `scripts/main/FOUR_SYSTEM_DELIVERY.m`
- **分离数据**：运行 `final/export_all_clean.py`
- **fig6 概念图**：运行 `Project-1/generate_figures.py`
