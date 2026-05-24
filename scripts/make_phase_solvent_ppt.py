#!/usr/bin/env python3
"""
Generate "Ternary Phase Diagrams Guide Green Solvent Selection" PPT
Uses all figures from Three-phase/output/figures/
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.enum.shapes import MSO_SHAPE
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.dirname(SCRIPT_DIR)
FIG_DIR = os.path.join(PROJ_DIR, 'output', 'figures')
OUT_PATH = os.path.join(PROJ_DIR, 'output', 'Phase_Solvent_Selection_Guide.pptx')

DARK  = RGBColor(0x1a, 0x1a, 0x2e)
PRIMARY = RGBColor(0x4f, 0x46, 0xe5)
GREEN  = RGBColor(0x10, 0xb9, 0x81)
AMBER  = RGBColor(0xf5, 0x9e, 0x0b)
RED    = RGBColor(0xef, 0x44, 0x44)
WHITE  = RGBColor(0xff, 0xff, 0xff)
GRAY   = RGBColor(0x66, 0x66, 0x66)
DARK_TEXT = RGBColor(0x1f, 0x29, 0x37)
LIGHT_BG = RGBColor(0xf8, 0xf9, 0xfa)
CYAN   = RGBColor(0x08, 0x91, 0xb2)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

def title_slide(title, subtitle=""):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    bg = slide.background; fill = bg.fill; fill.solid(); fill.fore_color.rgb = DARK
    tb = slide.shapes.add_textbox(Inches(1.5), Inches(2.0), Inches(10.3), Inches(1.5))
    p = tb.text_frame.paragraphs[0]; p.text = title
    p.font.size = Pt(38); p.font.bold = True; p.font.color.rgb = WHITE
    if subtitle:
        tb2 = slide.shapes.add_textbox(Inches(1.5), Inches(3.8), Inches(10.3), Inches(1.0))
        p2 = tb2.text_frame.paragraphs[0]; p2.text = subtitle
        p2.font.size = Pt(18); p2.font.color.rgb = RGBColor(0xaa,0xaa,0xbb)
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(1.5), Inches(3.5), Inches(2), Pt(4))
    bar.fill.solid(); bar.fill.fore_color.rgb = PRIMARY; bar.line.fill.background()
    return slide

def content_slide(title):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), prs.slide_width, Inches(0.9))
    bar.fill.solid(); bar.fill.fore_color.rgb = DARK; bar.line.fill.background()
    tb = slide.shapes.add_textbox(Inches(0.8), Inches(0.12), Inches(12), Inches(0.7))
    p = tb.text_frame.paragraphs[0]; p.text = title
    p.font.size = Pt(26); p.font.bold = True; p.font.color.rgb = WHITE
    line = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0), Inches(0.9), prs.slide_width, Pt(3))
    line.fill.solid(); line.fill.fore_color.rgb = PRIMARY; line.line.fill.background()
    return slide

def add_img(slide, filename, left, top, width, height=None):
    path = os.path.join(FIG_DIR, filename)
    if os.path.exists(path):
        if height is None:
            return slide.shapes.add_picture(path, Inches(left), Inches(top), Inches(width))
        else:
            return slide.shapes.add_picture(path, Inches(left), Inches(top), Inches(width), Inches(height))
    else:
        tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(0.5))
        tb.text_frame.paragraphs[0].text = f"[{filename}]"
        return tb

def bullet_box(slide, left, top, width, height, items, fs=15):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame; tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i==0 else tf.add_paragraph()
        p.text = item; p.font.size = Pt(fs); p.font.color.rgb = DARK_TEXT; p.space_after = Pt(5)
    return tf

def table(slide, left, top, col_w, headers, rows, fs=11):
    nr = len(rows)+1; nc = len(headers); tw = sum(col_w)
    ts = slide.shapes.add_table(nr, nc, Inches(left), Inches(top), Inches(tw), Inches(0.35*nr))
    tbl = ts.table
    for j, h in enumerate(headers):
        c = tbl.cell(0, j); c.text = h
        for p in c.text_frame.paragraphs: p.font.size=Pt(fs); p.font.bold=True; p.font.color.rgb=WHITE; p.alignment=PP_ALIGN.CENTER
        c.fill.solid(); c.fill.fore_color.rgb = DARK
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            c = tbl.cell(i+1, j); c.text = str(v)
            for p in c.text_frame.paragraphs: p.font.size=Pt(fs); p.font.color.rgb=DARK_TEXT; p.alignment=PP_ALIGN.CENTER
            if i%2==0: c.fill.solid(); c.fill.fore_color.rgb = LIGHT_BG
    for i, w in enumerate(col_w): tbl.columns[i].width = Inches(w)

# ============================================================
# SLIDE 1: TITLE
# ============================================================
title_slide("三元相图指导绿色溶剂选择",
            "从 Flory-Huggins 相图理解氯仿与绿色溶剂的效率差异  |  2026-05-22")

# ============================================================
# SLIDE 2: THE QUESTION
# ============================================================
s = content_slide("核心问题：绿色溶剂 vs 氯仿，差异从何而来？")
bullet_box(s, 0.8, 1.3, 5.8, 3.0, [
    "氯仿 (CF) 是 OPV 领域的\"黄金标准\"加工溶剂",
    "绿色溶剂 (Toluene, o-Xylene) 替代后 PCE 下降",
    "关键追问:",
    "  1. 效率衰减的本质是什么？热力学还是动力学？",
    "  2. 为什么 o-Xylene 在绿色溶剂中表现最好？",
    "  3. 三元相图能告诉我们什么？",
], fs=16)

bullet_box(s, 7.5, 1.3, 5.2, 5.5, [
    "器件效率数据:",
    "",
    "PM6:L8-Bo / Tol:   18.75%",
    "PM6:L8-Bo / OXy:  19.50%",
    "PM6:D18:L8-Bo / Tol:  19.58%",
    "PM6:D18:L8-Bo / OXy: 20.25%",
    "",
    "o-Xylene 在两个体系中均优于 Toluene",
    "D18 三元策略额外增益 +0.8%",
    "最佳体系: D18+OXy = 20.25%",
], fs=14)

# ============================================================
# SLIDE 3: ALL 6 SYSTEMS DATA
# ============================================================
s = content_slide("六体系临界点数据总览")
table(s, 0.8, 1.5, [1.8, 1.6, 1.4, 1.4, 1.4, 1.4, 2.0],
    ['体系', '给体', '溶剂', 'CP phi1', 'CP phi2', 'CP phi3', 'PCE %'],
    [
        ['PM6-CF',  'PM6', 'CF',     '0.455', '0.485', '0.061', '— (ref)'],
        ['PM6-Tol', 'PM6', 'Toluene','0.231', '0.718', '0.051', '18.75'],
        ['PM6-OXy','PM6', 'o-Xylene','0.272', '0.676', '0.052', '19.50'],
        ['D18-CF',  'D18', 'CF',     '0.390', '0.537', '0.073', '— (ref)'],
        ['D18-Tol', 'D18', 'Toluene','0.127', '0.819', '0.055', '19.58'],
        ['D18-OXy','D18', 'o-Xylene','0.503', '0.428', '0.069', '20.25'],
    ], fs=13)

bullet_box(s, 0.8, 4.5, 12, 2.5, [
    "CP φ₂ = 相分离起始门槛。φ₂ 越高 → 成膜过程中越早进入两相区 → 畴区粗化时间越长",
    "CF 的 CP φ₂ (0.485/0.537) 显著低于 Tol (0.718/0.819) → CF 延迟相分离 → 更可控形貌",
    "o-Xylene 的 CP φ₂ (0.676/0.428) 介于 CF 和 Tol 之间 → 比 Tol 更接近 CF 的相分离行为",
    "D18-OXy 的 CP φ₂=0.428 是所有体系中最低的 → 极端延迟相分离 + D18预聚集 = 20.25%",
], fs=14)

# ============================================================
# SLIDE 4: CP COMPARISON VISUALLY
# ============================================================
s = content_slide("临界点 φ₂ 对比：越低越接近 CF，越高越偏离")
add_img(s, 'fig5_phase_metrics.png', 0.8, 1.2, 5.8, 3.0)
add_img(s, 'fig4_red_pce.png', 7.0, 1.2, 5.8, 3.0)
bullet_box(s, 0.8, 4.5, 12, 2.5, [
    "左图: 四体系相图量化指标。D18-Tol 相分离最早 (CP=0.819), D18-OXy 最晚 (CP=0.428)",
    "右图: RED 溶剂质量 + PCE。o-Xylene 对 PM6 是更好的溶剂 (RED=0.63 vs 0.98)",
    "核心规律: CP φ₂ 越低 → 器件 PCE 越高 (D18-OXy 例外需要预聚集补偿)",
], fs=14)

# ============================================================
# SLIDE 5: TERNARY PHASE DIAGRAMS — Tol/OXy
# ============================================================
s = content_slide("三元相图：Toluene vs o-Xylene (四体系)")
add_img(s, 'ternary_real_2x2.png', 0.3, 1.1, 12.7, 6.0)

# ============================================================
# SLIDE 6: TERNARY PHASE DIAGRAMS — CF
# ============================================================
s = content_slide("三元相图：Chloroform 体系 (新增)")
add_img(s, 'ternary_real_CF_2x1.png', 1.0, 1.1, 11.3, 5.5)
bullet_box(s, 0.8, 6.5, 12, 0.7, [
    "PM6-CF: CP φ₂=0.485 (38 tie-lines) | D18-CF: CP φ₂=0.537 (23 tie-lines) | 参数为 HSP 估算值，未经完整扫描优化",
], fs=11)

# ============================================================
# SLIDE 7: THE PHYSICAL MECHANISM
# ============================================================
s = content_slide("物理机制：临界点决定相分离窗口宽度")
add_img(s, 'fig6_drying_trajectory.png', 0.5, 1.1, 7.0, 4.2)

bullet_box(s, 8.0, 1.3, 5.0, 5.5, [
    "成膜过程 = 溶剂挥发 → 浓度↑",
    "→ φ₂ 下降 → 触及 CP → LLPS",
    "",
    "CF (CP≈0.49):",
    "  触发晚 + 挥发快 (bp 61°C)",
    "  → 极速穿越两相区",
    "  → 形貌被\"冻结\"在最优状态",
    "",
    "Tol (CP≈0.77):",
    "  触发极早 + 挥发中等",
    "  → LLPS 区停留极长",
    "  → 畴区过度粗化",
    "",
    "OXy + D18 (CP≈0.43):",
    "  触发最晚 + D18 预聚集",
    "  → 即使挥发慢，畴区仍可控",
], fs=13)

# ============================================================
# SLIDE 8: WHY o-XYLENE WINS
# ============================================================
s = content_slide("为什么 D18 + o-Xylene = 20.25%？")
add_img(s, 'fig7_four_model.png', 0.5, 1.1, 8.5, 4.5)

bullet_box(s, 9.5, 1.3, 3.5, 5.5, [
    "四层协同优势:",
    "",
    "1. HSP dh=3.10",
    "   更匹配 L8-Bo",
    "",
    "2. CP φ₂=0.428",
    "   最晚相分离",
    "",
    "3. D18 预聚集",
    "   跳过 LLPS 成核",
    "",
    "4. 给体相纯度",
    "   φ₃=0.530",
    "",
    "→ 20.25% PCE",
], fs=13)

# ============================================================
# SLIDE 9: THREE-LAYER FRAMEWORK
# ============================================================
s = content_slide("三元统一框架：HSP → 相图 → 动力学")
table(s, 0.8, 1.5, [1.5, 2.5, 3.5, 3.5],
    ['层次', '工具', 'CF 表现', '绿色溶剂设计启示'],
    [
        ['1. 热力学', 'HSP / Flory-Huggins',
         'δh=5.70 匹配 L8-Bo', '选择 δh=3-6 的溶剂或共混物'],
        ['2. 相图', '三元相图 CP φ₂',
         'CP≈0.49 延迟相分离', '目标 CP φ₂=0.40-0.55'],
        ['3. 动力学', '挥发速率 + 预聚集',
         '快挥发天然短LLPS', '慢挥发需预聚集补偿'],
    ], fs=13)

add_img(s, 'fig8_design_compass.png', 1.5, 4.0, 10.0, 3.0)

# ============================================================
# SLIDE 10: GREEN SOLVENT DESIGN RULES
# ============================================================
s = content_slide("绿色溶剂理性设计判据")
table(s, 0.8, 1.5, [2.2, 2.0, 3.0, 4.5],
    ['判据', '推荐范围', 'Tol 评估', 'OXy 评估'],
    [
        ['CP φ₂', '0.40–0.55', '0.72–0.82 太高 ✗', '0.43–0.68 接近 ✓'],
        ['浓相 φ₃ max', '> 0.45', '0.38 不足 ✗', '0.50–0.53 达标 ✓'],
        ['g23 (溶剂-给体)', '0.20–0.35', '0.39–0.45 偏高', '0.30–0.41 适中'],
        ['溶剂 δh', '3.0–6.0', '2.0 太低 ✗', '3.1 达标 ✓'],
        ['PM6/溶剂 Ra', '< 8.0', '7.8 边缘 ✓', '5.1 达标 ✓'],
    ], fs=13)

bullet_box(s, 0.8, 5.0, 12, 2.0, [
    "理想绿色溶剂画像: CP φ₂ ≈ 0.45–0.50 + bp 80–120°C + δh ≈ 3–6",
    "o-Xylene + D18 预聚集是目前最接近这一画像的绿色溶剂方案",
], fs=14)

# ============================================================
# SLIDE 11: CONCLUSIONS
# ============================================================
s = content_slide("结论")
bullet_box(s, 0.8, 1.3, 6.0, 5.5, [
    "1. 三元相图的 CP φ₂ 是预测溶剂适用性的核心参数",
    "   → 直接决定相分离触发时机 → 控制畴区尺寸",
    "",
    "2. CF 优势 = 低 CP + 快挥发的完美配合",
    "   → 相分离晚触发 + 快速穿越 = 最优形貌",
    "",
    "3. Tol 劣势 = 相分离过早触发 (CP≈0.77)",
    "   → LLPS 区停留过长 → 畴区过度粗化",
    "",
    "4. o-Xy + D18 通过\"晚触发+空间补偿\"逼近 CF",
    "   → 最低 CP (0.428) + 预聚集补偿慢挥发",
    "   → 20.25% = 绿色溶剂最高效率之一",
], fs=15)

bullet_box(s, 7.5, 1.3, 5.0, 5.5, [
    "5. 理性溶剂设计三原则:",
    "",
    "  (a) HSP 初筛:",
    "      溶剂混合物落入 CF 邻域",
    "",
    "  (b) 相图精选:",
    "      CP φ₂ 在 0.40–0.55 范围",
    "",
    "  (c) 动力学补偿:",
    "      刚性给体第三组分 + 预聚集",
    "",
    "  → 可推广至任意新型绿色溶剂",
], fs=15)

# ============================================================
# SLIDE 12: HSP CONTEXT
# ============================================================
s = content_slide("HSP 分析：热力学基础")
add_img(s, 'fig1_hsp_bubble.png', 0.3, 1.1, 6.5, 3.0)
add_img(s, 'fig3_ra_heatmap.png', 7.0, 1.1, 5.5, 3.0)
bullet_box(s, 0.8, 4.5, 12, 2.5, [
    "左: HSP 二维气泡图。CF (δh=5.70) 与 L8-Bo (δh=5.93) 完美匹配。o-Xy (δh=3.10) > Tol (δh=2.00)",
    "右: Ra 距离矩阵。PM6/L8-Bo Ra=6.31 → 适中相容性差异利于 BHJ 形貌",
], fs=14)

# ============================================================
# SLIDE 13: ACKNOWLEDGMENT
# ============================================================
title_slide("谢谢！",
            "计算方法: Flory-Huggins 三元相平衡 (MATLAB + Python)\n"
            "HSP: Hoftyzer-van Krevelen 基团贡献法\n"
            "数据: 六体系 × 临界点 × Binodal × Spinodal × Tie-lines")

# ============================================================
# SAVE
# ============================================================
prs.save(OUT_PATH)
print(f"PPT saved: {OUT_PATH}")
print(f"Slides: {len(prs.slides)}")
