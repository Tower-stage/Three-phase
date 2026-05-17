# 维基
https://deepwiki.com/KaihangShi/Ternary-Phase-Diagram

# Ternary-Phase-Diagram

这段代码是基于**Flory-Huggins高分子溶液理论**，用于计算和绘制**三元相图（聚合物 - 溶剂 - 非溶剂/小分子）**的MATLAB脚本。

当前仓库已针对 **PM6（聚合物给体）- L8-Bo（小分子受体）- 溶剂（Toluene / o-Xylene）** 体系进行了重构与参数化。

---

## 项目结构（重构后）

```
src/
├── binodal-PM6-L8Bo-Tol/      # 甲苯体系 Binodal + Tie line 计算
│   ├── main.m
│   └── fun.m
├── binodal-PM6-L8Bo-OXy/      # 邻二甲苯体系 Binodal + Tie line 计算
│   ├── main.m
│   └── fun.m
├── spinodal-gemini/           # 旋节线计算（网格算法，已解决触底问题）
│   ├── main.m
│   └── fun.m
├── Critical_point/            # 临界点计算
│   ├── find_critical_point_Tol.m
│   └── find_critical_point_OXy.m
├── binodal-Gemini/            # Gemini 重构版模板（DMF-THF-PIM1 参考体系）
├── legacy/                    # 旧版代码存档（已弃用）
│   ├── main_original.m / fun_original.m
│   ├── Spinodal_original/
│   ├── LLE_original/
│   └── Critical_point_original_*.m
└── readme.md                  # 参数速查表与调试指南
```

> **注意**：旧版 `src/Spinodal/`、`src/LLE/` 已移入 `legacy/`。前者已被 `spinodal-gemini` 的网格算法完全替代；后者的单次调试功能已包含在 Binodal 循环中。

---

## 快速开始：计算完整四图元相图

你的目标是得到 **Critical Point + Spinodal + Binodal + Tie line**。

### 1. Spinodal（旋节线）

```matlab
cd src/spinodal-gemini
main
```

采用**网格算法**直接提取 `det(G'')=0` 的等高线，无需担心"触底"问题。运行后数据存储在变量 `xsol` 中。

### 2. Binodal + Tie line（双节点线 + 结线）

```matlab
cd src/binodal-PM6-L8Bo-Tol   % 或 binodal-PM6-L8Bo-OXy
main
```

运行后生成 `binodal_PM6_L8Bo_*.csv`。**CSV 中每相邻两行即构成一条 Tie line**（第 `2i-1` 行为浓相，第 `2i` 行为稀相），无需单独计算 Tie line。

### 3. Critical Point（临界点）

```matlab
cd src/Critical_point
find_critical_point_Tol   % 或 find_critical_point_OXy
```

策略：**先计算 Spinodal 曲线，再在 Spinodal 上寻找三阶条件变号点**，通过线性插值得到精确临界点。结果保存为 `critical_point_*.csv`。

---

## 关键参数速查

| 参数 | Tol 体系 | O-Xy 体系 | 物理意义 |
|---|---|---|---|
| `v1` | 1132.1 | 1132.1 | L8-Bo 摩尔体积 [cm³/mol] |
| `v2` | 106.3 | 120.6 | 溶剂摩尔体积 |
| `v3` | 1743900 | 1743900 | PM6 摩尔体积 (DP≈1500) |
| `X13` | **0.62** | **0.65** | L8-Bo/PM6 相互作用参数（已调优） |
| `g23` (bb) | 0.3852 | 0.4120 | 溶剂/PM6 相互作用参数 |
| `p1-p5` | [见下 Tol] | [见下 O-Xy] | L8-Bo/溶剂 g₁₂ 多项式系数 |

**Tol 体系 g₁₂ 系数：** `0, 0, 0, -0.2000, 0.6000`（HSP 估算 + 线性组成依赖）

**O-Xy 体系 g₁₂ 系数：** `0, 0, 0, -0.2000, 0.7500`（HSP 估算 + 线性组成依赖）

---

## 调试指南

- **Binodal 大量迭代失败**：微调 `x0_default`（初始猜测值），或检查 `iX13` 是否合理。
- **`iX13` 取值**：Binodal 计算推荐 **Tol 0.62 / O-Xy 0.65**。`spinodal-gemini` 中实际调试用 **0.447**。若相图形状异常可尝试在 [0.55, 0.70] 范围内微调。
- **临界点未找到**：若 `crit_val` 未变号，说明当前参数下临界点可能落在网格边缘，可尝试增大网格分辨率（`linspace(..., 1200)`）或微调 `X13`。
- **可视化**：所有 CSV 均可直接导入 **Origin** → Plot → Ternary 绘制三元相图。

---

## 历史修改记录

| Parameter | File | Description |
|---|---|---|
| `v1`, `v2`, `v3` | `fun.m` | 摩尔体积 |
| `g23` (bb) | `fun.m` | 溶剂–聚合物作用参数 |
| `p1`–`p5` | `fun.m` | 非溶剂–溶剂 `g12` 多项式系数 |
| `iX13` | `main.m` | 非溶剂–聚合物作用参数 |
| `x0` | `main.m` | 初始猜测（不收敛时调整） |

---

## 基团数量收集

groups_Tol = {9: 5, 10: 1}  # Toluene: 5*AC, 1*ACCH3
groups_OXY = {9: 4, 10: 2}  # o-Xylene: 4*AC, 2*ACCH3
groups_L8-Bo = {1: 8, 2: 28, 1: 8, 3: 4, 9: 12, 11: 4, 18: 2, 25: 4, 39: 5, 49: 4}

SMILES_L8-Bo = CCCCCCC(CCCC)Cc1c(C=C2C(=O)c3cc(F)c(F)cc3C2=C(C#N)C#N)sc2c1sc1c3c4nsnc4c4c5sc6c(CC(CCCC)CCCCCC)c(C=C7C(=O)c8cc(F)c(F)cc8C7=C(C#N)C#N)sc6c5n(CC(CC)CCCC)c4c3n(CC(CC)CCCC)c21

## P值收集
L8 : Toluene
21.0142 -53.8880 50.9171 -24.4012 -1.8878
L8 : o-Xylene
16.6738 -40.8682 39.6207 -19.3621 -1.5993
