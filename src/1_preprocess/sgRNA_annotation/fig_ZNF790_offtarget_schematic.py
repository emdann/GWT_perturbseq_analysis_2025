"""
Schematic: ZNF790-2 guide off-target binding at CD3D promoter region.

Key data:
  ZNF790-2 spacer:  5'-GCGGAACCTTCCACCATCAA-3'  (20 nt, on-target chr19:ZNF790)
  Seed (3' 10 nt):     CCACCATCAA
  Off-target site:  chr11:118,342,752, minus strand
  CD3D TSS:         chr11:118,342,705, minus strand  (47 bp away)
  CD3G TSS:         chr11:118,344,344, plus strand   (~1.6 kb away)
  Seed match length: 10 bp,  Hamming distance: 7 (full 20-mer)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# ── Coordinates ──
CD3D_TSS = 118_342_705
CD3G_TSS = 118_344_344
OT_POS = 118_342_752  # off-target match position (start of seed on genome)

# Display window
WIN_START = 118_341_500
WIN_END = 118_345_000

# ── Figure setup ──
fig, ax = plt.subplots(figsize=(10, 5.5))
ax.set_xlim(WIN_START, WIN_END)
ax.set_ylim(-1.5, 4.5)
ax.axis('off')

# ── Color palette ──
COL_CHROM = '#b0b0b0'
COL_CD3D = '#2166ac'
COL_CD3G = '#b2182b'
COL_GUIDE = '#e08214'
COL_SEED = '#d62728'
COL_MISMATCH = '#999999'
COL_PROMOTER = '#d1e5f0'

# ── 1. Chromosome backbone ──
chrom_y = 1.5
ax.plot([WIN_START, WIN_END], [chrom_y, chrom_y], color=COL_CHROM,
        linewidth=3, solid_capstyle='round', zorder=1)
ax.text((WIN_START + WIN_END) / 2, chrom_y - 0.35, 'chr11', fontsize=9,
        ha='center', va='top', color='#666666', style='italic')

# ── 2. Gene bodies ──
gene_h = 0.22

# CD3D (minus strand) — draw from TSS leftward
cd3d_body_start = CD3D_TSS - 1200
cd3d_body_end = CD3D_TSS
ax.add_patch(FancyBboxPatch((cd3d_body_start, chrom_y - gene_h / 2),
             cd3d_body_end - cd3d_body_start, gene_h,
             boxstyle="round,pad=30", facecolor=COL_CD3D, edgecolor='none',
             alpha=0.7, zorder=2))
# Arrow for strand direction (minus → leftward)
ax.annotate('', xy=(cd3d_body_start - 80, chrom_y),
            xytext=(cd3d_body_start + 150, chrom_y),
            arrowprops=dict(arrowstyle='->', color=COL_CD3D, lw=2), zorder=3)
ax.text((cd3d_body_start + cd3d_body_end) / 2, chrom_y + 0.35,
        'CD3D', fontsize=12, fontweight='bold', ha='center', va='bottom',
        color=COL_CD3D, zorder=4)
ax.text((cd3d_body_start + cd3d_body_end) / 2, chrom_y + 0.15,
        '(− strand)', fontsize=7, ha='center', va='bottom',
        color=COL_CD3D, alpha=0.7, zorder=4)

# CD3G (plus strand) — draw from TSS rightward
cd3g_body_start = CD3G_TSS
cd3g_body_end = CD3G_TSS + 1000
ax.add_patch(FancyBboxPatch((cd3g_body_start, chrom_y - gene_h / 2),
             cd3g_body_end - cd3g_body_start, gene_h,
             boxstyle="round,pad=30", facecolor=COL_CD3G, edgecolor='none',
             alpha=0.7, zorder=2))
ax.annotate('', xy=(cd3g_body_end + 80, chrom_y),
            xytext=(cd3g_body_end - 150, chrom_y),
            arrowprops=dict(arrowstyle='->', color=COL_CD3G, lw=2), zorder=3)
ax.text((cd3g_body_start + cd3g_body_end) / 2, chrom_y + 0.35,
        'CD3G', fontsize=12, fontweight='bold', ha='center', va='bottom',
        color=COL_CD3G, zorder=4)
ax.text((cd3g_body_start + cd3g_body_end) / 2, chrom_y + 0.15,
        '(+ strand)', fontsize=7, ha='center', va='bottom',
        color=COL_CD3G, alpha=0.7, zorder=4)

# ── 3. TSS markers ──
tss_tick = 0.3
for tss, col, label in [(CD3D_TSS, COL_CD3D, 'TSS'),
                         (CD3G_TSS, COL_CD3G, 'TSS')]:
    ax.plot([tss, tss], [chrom_y - tss_tick, chrom_y + tss_tick],
            color=col, linewidth=2, zorder=3)

# ── 4. Promoter region (±2 kb around CD3D TSS) ──
prom_start = CD3D_TSS - 2000
prom_end = CD3D_TSS + 2000
ax.add_patch(FancyBboxPatch((prom_start, chrom_y - 0.45), prom_end - prom_start, 0.9,
             boxstyle="round,pad=40", facecolor=COL_PROMOTER, edgecolor=COL_CD3D,
             alpha=0.25, linestyle='--', linewidth=1.2, zorder=0))
ax.text(prom_start + 100, chrom_y - 0.55, 'CD3D promoter region (±2 kb)',
        fontsize=7, color=COL_CD3D, alpha=0.7, va='top')

# ── 5. Off-target binding site ──
ot_width = 200  # visual width for the binding site
ot_y = chrom_y + 0.6
ax.add_patch(FancyBboxPatch((OT_POS - ot_width / 2, ot_y), ot_width, 0.25,
             boxstyle="round,pad=20", facecolor=COL_SEED, edgecolor=COL_SEED,
             alpha=0.3, zorder=2))

# ── 6. Guide RNA ──
guide_y = 3.2
guide_width = 600  # visual width
guide_x = OT_POS - guide_width / 2

# Guide body
ax.add_patch(FancyBboxPatch((guide_x, guide_y - 0.15), guide_width, 0.3,
             boxstyle="round,pad=20", facecolor=COL_GUIDE, edgecolor='#c06000',
             alpha=0.85, linewidth=1.5, zorder=3))

# Seed region (3' end, right portion)
seed_frac = 10 / 20  # 10 of 20 nt
seed_x = guide_x + guide_width * (1 - seed_frac)
ax.add_patch(FancyBboxPatch((seed_x, guide_y - 0.15), guide_width * seed_frac, 0.3,
             boxstyle="round,pad=20", facecolor=COL_SEED, edgecolor='#a01010',
             alpha=0.6, linewidth=1.5, zorder=4))

# Labels on guide
ax.text(guide_x + guide_width * 0.25, guide_y, "5'", fontsize=8,
        ha='center', va='center', color='white', fontweight='bold', zorder=5)
ax.text(guide_x + guide_width * 0.75, guide_y, "3'", fontsize=8,
        ha='center', va='center', color='white', fontweight='bold', zorder=5)
ax.text(OT_POS, guide_y + 0.35, 'ZNF790-2 sgRNA', fontsize=11,
        fontweight='bold', ha='center', va='bottom', color='#915b00', zorder=5)

# ── 7. Sequence annotation ──
seq_y = 2.55
spacer = "5'-GCGGAACCTTCCACCATCAA-3'"
ax.text(OT_POS, seq_y, spacer, fontsize=8, ha='center', va='center',
        fontfamily='monospace', color='#333333', zorder=5,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff8e8',
                  edgecolor='#e0c080', alpha=0.9))

# Seed label
ax.annotate('10-nt seed match', xy=(seed_x + guide_width * seed_frac / 2, guide_y - 0.15),
            xytext=(OT_POS + 500, guide_y - 0.5),
            fontsize=8, color=COL_SEED, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=COL_SEED, lw=1.2),
            ha='center', va='top', zorder=5)

# Mismatch region label
ax.annotate('7 mismatches\n(PAM-distal)',
            xy=(guide_x + guide_width * 0.25, guide_y - 0.15),
            xytext=(OT_POS - 600, guide_y - 0.5),
            fontsize=8, color=COL_MISMATCH, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=COL_MISMATCH, lw=1.2),
            ha='center', va='top', zorder=5)

# ── 8. Binding arrow (guide → genome) ──
ax.annotate('', xy=(OT_POS, ot_y + 0.25),
            xytext=(OT_POS, guide_y - 0.15),
            arrowprops=dict(arrowstyle='->', color=COL_SEED, lw=2,
                            linestyle='--', alpha=0.7),
            zorder=3)

# ── 9. Distance annotation ──
dist_y = chrom_y + 1.05
ax.annotate('', xy=(CD3D_TSS, dist_y), xytext=(OT_POS, dist_y),
            arrowprops=dict(arrowstyle='<->', color='#444444', lw=1))
ax.text((CD3D_TSS + OT_POS) / 2, dist_y + 0.08, '47 bp from TSS',
        fontsize=7, ha='center', va='bottom', color='#444444')

# ── 10. Scale bar ──
scale_y = chrom_y - 0.85
scale_start = WIN_END - 1500
scale_end = WIN_END - 500
ax.plot([scale_start, scale_end], [scale_y, scale_y], color='black', linewidth=1.5)
ax.plot([scale_start, scale_start], [scale_y - 0.05, scale_y + 0.05], color='black', linewidth=1.5)
ax.plot([scale_end, scale_end], [scale_y - 0.05, scale_y + 0.05], color='black', linewidth=1.5)
ax.text((scale_start + scale_end) / 2, scale_y - 0.1, '1 kb', fontsize=8,
        ha='center', va='top')

# ── 11. Legend ──
legend_elements = [
    mpatches.Patch(facecolor=COL_GUIDE, edgecolor='#c06000', alpha=0.85,
                   label='sgRNA (PAM-distal)'),
    mpatches.Patch(facecolor=COL_SEED, edgecolor='#a01010', alpha=0.6,
                   label='Seed match region (PAM-proximal)'),
    mpatches.Patch(facecolor=COL_CD3D, alpha=0.7, label='CD3D gene body'),
    mpatches.Patch(facecolor=COL_CD3G, alpha=0.7, label='CD3G gene body'),
    mpatches.Patch(facecolor=COL_PROMOTER, edgecolor=COL_CD3D, alpha=0.3,
                   linestyle='--', label='CD3D promoter (±2 kb)'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=7,
          framealpha=0.9, edgecolor='#cccccc')

# ── Title ──
ax.set_title('ZNF790-2 sgRNA off-target binding at CD3D promoter\n'
             '(chr11: CD3D–CD3G locus)',
             fontsize=13, fontweight='bold', pad=10)

plt.tight_layout()
plt.savefig('results/fig_ZNF790_offtarget_schematic.pdf', dpi=300, bbox_inches='tight')
plt.savefig('results/fig_ZNF790_offtarget_schematic.png', dpi=300, bbox_inches='tight')
plt.show()
print('Saved to results/fig_ZNF790_offtarget_schematic.{pdf,png}')
