#!/usr/bin/env python3
"""Consolidate all ternary phase diagram data into one Excel file."""
import os, csv
import openpyxl
from openpyxl.utils import get_column_letter

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'output', 'data')
OUT_PATH = os.path.join(os.path.dirname(__file__), 'Ternary_All_Data_v2.xlsx')

SYSTEMS = [
    ('PM6_Tol', 'PM6 / L8-Bo / Toluene'),
    ('PM6_OXy', 'PM6 / L8-Bo / o-Xylene'),
    ('D18_Tol', 'D18 / L8-Bo / Toluene'),
    ('D18_OXy', 'D18 / L8-Bo / o-Xylene'),
    ('PM6_CF',  'PM6 / L8-Bo / Chloroform'),
    ('D18_CF',  'D18 / L8-Bo / Chloroform'),
]

wb = openpyxl.Workbook()
# Remove default sheet
wb.remove(wb.active)

for tag, title in SYSTEMS:
    ws = wb.create_sheet(title=tag)

    # ===== BINODAL =====
    bin_file = os.path.join(DATA_DIR, f'binodal_{tag}.csv')
    bino_data = []
    if os.path.exists(bin_file):
        with open(bin_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
            # Only concentrated phase: odd rows (0, 2, 4, ...)
            for i in range(0, len(rows)-1, 2):
                row = rows[i]
                if len(row) >= 4:
                    bino_data.append([float(row[1]), float(row[3]), float(row[2])])  # phi1, phi2, phi3

    # ===== SPINODAL =====
    spin_file = os.path.join(DATA_DIR, f'spinodal_{tag}.csv')
    spin_data = []
    if os.path.exists(spin_file):
        with open(spin_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            for row in reader:
                if len(row) >= 4:
                    spin_data.append([float(row[1]), float(row[3]), float(row[2])])  # phi1, phi2, phi3

    # ===== CRITICAL =====
    crit_file = os.path.join(DATA_DIR, f'critical_{tag}.csv')
    crit_vals = None
    if os.path.exists(crit_file):
        with open(crit_file, 'r') as f:
            reader = csv.reader(f)
            next(reader)
            row = next(reader)
            crit_vals = [float(row[0]), float(row[1]), float(row[2])]

    # ===== WRITE HEADERS =====
    # Row 1: Title
    ws.merge_cells('A1:I1')
    ws['A1'] = title
    ws['A1'].font = openpyxl.styles.Font(bold=True, size=14)

    # Row 2: System info
    ws.merge_cells('A2:I2')
    n_tls = len(bino_data)  # concentrated phase only
    n_spin = len(spin_data)
    info_parts = [f'{n_tls} tie-line pairs', f'{n_spin} spinodal points']
    if crit_vals:
        info_parts.append(f'CP: phi1={crit_vals[0]:.4f}, phi2={crit_vals[1]:.4f}, phi3={crit_vals[2]:.4f}')
    ws['A2'] = ' | '.join(info_parts)
    ws['A2'].font = openpyxl.styles.Font(color='666666', size=10)

    # Row 4: Column headers
    headers = ['phi1_L8Bo', 'phi2_Solvent', 'phi3_Donor',
               '',  # spacer
               'phi1_L8Bo', 'phi2_Solvent', 'phi3_Donor']
    header_notes = ['BINODAL', '', '', '', 'SPINODAL', '', '']
    for j, h in enumerate(headers):
        if h:
            cell = ws.cell(row=4, column=j+1, value=h)
            cell.font = openpyxl.styles.Font(bold=True, size=11)
            cell.fill = openpyxl.styles.PatternFill(start_color='4472C4' if j < 3 else 'ED7D31',
                                                     end_color='4472C4' if j < 3 else 'ED7D31',
                                                     fill_type='solid')
            cell.font = openpyxl.styles.Font(bold=True, color='FFFFFF', size=10)
    # Section labels
    ws.merge_cells('A3:C3')
    ws['A3'] = 'BINODAL (concentrated phase only)'
    ws['A3'].font = openpyxl.styles.Font(bold=True, color='4472C4', size=11)
    ws.merge_cells('F3:H3')
    ws['F3'] = 'SPINODAL'
    ws['F3'].font = openpyxl.styles.Font(bold=True, color='ED7D31', size=11)

    # ===== WRITE DATA =====
    # Binodal (cols A-C, starting row 5)
    for i, row in enumerate(bino_data):
        for j, val in enumerate(row):
            ws.cell(row=i+5, column=j+1, value=round(val, 6))

    # Spinodal (cols F-H, starting row 5)
    for i, row in enumerate(spin_data):
        for j, val in enumerate(row):
            ws.cell(row=i+5, column=j+6, value=round(val, 6))

    # Critical point note at bottom
    if crit_vals:
        last_row = max(len(bino_data), len(spin_data)) + 6
        ws.merge_cells(f'A{last_row}:C{last_row}')
        ws.cell(row=last_row, column=1,
                value=f'Critical Point: phi1={crit_vals[0]:.4f}  phi2={crit_vals[1]:.4f}  phi3={crit_vals[2]:.4f}')
        ws.cell(row=last_row, column=1).font = openpyxl.styles.Font(bold=True, color='C00000')

    # Column widths
    ws.column_dimensions['A'].width = 14
    ws.column_dimensions['B'].width = 14
    ws.column_dimensions['C'].width = 16
    ws.column_dimensions['D'].width = 3
    ws.column_dimensions['F'].width = 14
    ws.column_dimensions['G'].width = 14
    ws.column_dimensions['H'].width = 16

    print(f"  Sheet '{tag}': {n_tls} TL pairs, {n_spin} spin pts, CP phi2={crit_vals[1]:.4f}" if crit_vals else f"  Sheet '{tag}': {n_tls} TL pairs, {n_spin} spin pts")

wb.save(OUT_PATH)
print(f"\nSaved: {OUT_PATH}")
