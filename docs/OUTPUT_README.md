# Three-Phase 三元相图计算 — 最终交付说明

## 运行方式

```matlab
cd d:\Deskkk\First CC\Three-phase
FINAL_DELIVERY
```

运行后自动完成 Spinodal → Critical Point → Binodal → Tie Lines → 出图。

## 输出文件

### CSV 数据 (可直接导入 Origin)

| 文件 | 内容 | 格式 |
|------|------|------|
| `output/data/binodal_Tol.csv` | Toluene 体系双节点线 | 每相邻两行 = 一条 tie line（浓相 / 稀相） |
| `output/data/binodal_OXy.csv` | o-Xylene 体系双节点线 | 同上 |
| `output/data/spinodal_Tol.csv` | Toluene 体系旋节线 | φ1, φ3, φ2 |
| `output/data/spinodal_OXy.csv` | o-Xylene 体系旋节线 | φ1, φ3, φ2 |
| `output/data/critical_Tol.csv` | Toluene 体系临界点 | φ1, φ2, φ3 |
| `output/data/critical_OXy.csv` | o-Xylene 体系临界点 | φ1, φ2, φ3 |

### 图片

| 文件 | 说明 |
|------|------|
| `output/figures/ternary_combined.png` | 双体系综合三元相图 |
| `output/figures/ternary_Tol.png` | Tol 体系单独三元相图 |
| `output/figures/ternary_OXy.png` | O-Xy 体系单独三元相图 |
| `output/figures/ratio_Tol.png` | Tol 体系 Ratio-Solvent 相图（新可视化） |
| `output/figures/ratio_OXy.png` | O-Xy 体系 Ratio-Solvent 相图 |
| `output/figures/*_HR.png` | 以上所有图的 300 DPI 高清版本 |

## 优化后的参数

| 参数 | Toluene | o-Xylene | 物理意义 |
|------|---------|----------|---------|
| v3 | **100000** | **120000** | PM6 摩尔体积（DP≈86 / 103） |
| X13 | **0.95** | **0.80** | L8-Bo / PM6 相互作用 |
| g23 | **0.385** | **0.412** | 溶剂 / PM6 相互作用 |
| v1 | 1132.1 | 1132.1 | L8-Bo 摩尔体积 |
| v2 | 106.3 | 120.6 | 溶剂摩尔体积 |
| g12 系数 | [0,0,0,-0.20,0.60] | [0,0,0,-0.20,0.75] | L8-Bo/溶剂多项式 |

## 计算质量

| 指标 | Toluene | o-Xylene |
|------|:---:|:---:|
| Binodal 成功率 | 95% (114/120) | 99% (119/120) |
| 浓相单调性 | ✓ | ✓ |
| 稀相单调性 | ✓ | ✓ |
| 残差范围 | 0 ~ 2.4e-24 | 0 ~ 2.1e-23 |
| 浓相 PM6 跨 | 0.384 → 0.049 | 0.505 → 0.050 |
| 临界点 φ₃ | 0.0506 | 0.0522 |

## 参数来源

基于 192 组三维网格扫描 (v3 × X13 × g23)，以 Tol 体系为基准：

- **原始参数** (v3=1,743,900, X13=0.62): 相图过度挤压在 L8-Bo 侧
- **优化后** (v3=100,000, X13=0.95): 相图向中心展开，两臂分布更平衡
- **O-Xy** 通过 60 组针对性扫描优化至 v3=120,000, X13=0.80

### 调节规律

| 操作 | 效果 |
|------|------|
| 降低 v3 | 相图向中心拉宽，两臂远离 L8-Bo 边 |
| 提高 X13 | 相分离区扩大，tie line 变长 |
| 提高 g23 | 补偿 v3 降低时损失的分相驱动力 |

## 文件结构

```
Three-phase/
├── scripts/
│   ├── main/                 # 主运行脚本
│   │   └── FOUR_SYSTEM_DELIVERY.m
│   ├── scan/                 # 参数扫描脚本
│   │   ├── param_sweep.m
│   │   ├── param_sweep_v2.m
│   │   ├── find_best_OXy.m
│   │   ├── scan_D18_Tol.m
│   │   └── scan_D18_OXy.m
│   ├── legacy/               # 旧版交付脚本
│   │   ├── FINAL_DELIVERY.m
│   │   └── FINAL_deliver.m
│   └── debug/                # 调试与测试脚本
│       └── test_D18_OXy.m
├── docs/
│   ├── OUTPUT_README.md      # 本文件
│   ├── README.md
│   ├── LEARNING_GUIDE.md
│   └── 调试说明.md
├── output/
│   ├── data/                 # CSV 数据文件
│   ├── figures/              # PNG 图片输出
│   └── temp/                 # MAT 临时结果
├── src/                      # 原始源码与工具
│   ├── binodal-PM6-L8Bo-Tol/
│   ├── binodal-PM6-L8Bo-OXy/
│   ├── spinodal-gemini/
│   ├── Critical_point/
│   └── debug-tools/
└── reports/                  # 物理解读报告
```
