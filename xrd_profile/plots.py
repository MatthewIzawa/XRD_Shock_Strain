"""
Shared plotting primitives for shock-XRD manuscript figures.

Extracted from Paper1_JAC/regenerate_all_figures_v3.py. Every figure script
in Paper1_JAC/ and HED_XRD_Shock/ should use these constants and helpers —
if a colour, marker, or default changes, it changes here and both papers
pick it up.

This module is intentionally NOT exported via xrd_profile.__init__.__all__;
import explicitly as `from xrd_profile.plots import ...`.
"""
import csv
import numpy as np
from matplotlib.lines import Line2D
from adjustText import adjust_text


# ─── Canonical group colour / marker assignments ────────────────────────
GROUP_COLOUR = {
    'Lunar':      '#2166ac',
    'HED':        '#b2182b',
    'OC':         '#4dac26',
    'CC':         '#7570b3',
    'Achondrite': '#e08214',
    'Martian':    '#1b7837',
    # Ungrouped achondrites (e.g. NWA 6693, CC-associated impact melt).
    # Colour + marker match the dark-grey 'X' used in fig7's NWA 6693
    # PDF panel so the encoding is consistent across the figure set.
    'Ungrouped':  '#4d4d4d',
}
GROUP_MARKER = {
    'Lunar':      'o',
    'HED':        's',
    'OC':         '^',
    'CC':         'D',
    'Achondrite': 'P',
    'Martian':    '*',
    'Ungrouped':  'X',
}
# Source → marker edge colour (legacy; encoded instrument provenance via
# black-edge-vs-no-edge — too subtle at publication scale, per Jenkins's
# Fig 1 review comment 2026-05-16). Retained as alias for any external
# code that still imports it; new code should use SOURCE_FILL below.
SOURCE_EDGE = {
    'Lab Cu': 'none',
    'Lab Co': 'none',
    'I11':    'k',
}
# Source → marker fill style (current; encodes instrument provenance via
# hollow-vs-filled, far more visually distinct at publication scale).
# Laboratory = hollow (group colour in edge, no fill);
# Synchrotron = filled (group colour interior, black edge).
SOURCE_FILL = {
    'Lab Cu': 'none',      # hollow
    'Lab Co': 'none',      # hollow
    'I11':    'filled',    # filled
}
GROUP_ORDER = ['Lunar', 'HED', 'OC', 'CC', 'Achondrite', 'Martian', 'Ungrouped']

# Short aliases matching the v3 script for drop-in compatibility.
GC = GROUP_COLOUR
GM = GROUP_MARKER
SE = SOURCE_EDGE
SF = SOURCE_FILL


# ─── Legend handles ─────────────────────────────────────────────────────
def legend_handles(include_source=True, groups=None):
    """Legend handle list. Group handles first, optional Lab/Synchrotron
    handles after. `groups` overrides the default ordering."""
    h = []
    for g in (groups or GROUP_ORDER):
        h.append(Line2D([0], [0], marker=GROUP_MARKER[g], color='w',
                        markerfacecolor=GROUP_COLOUR[g], markersize=7,
                        markeredgecolor='grey', markeredgewidth=0.4, label=g))
    if include_source:
        # Lab: hollow marker (group colour would go in the edge; legend
        # shows grey edge for the generic Lab handle).
        h.append(Line2D([0], [0], marker='o', color='w',
                        markerfacecolor='none', markeredgecolor='grey',
                        markeredgewidth=1.2, markersize=7, label='Lab'))
        # Synchrotron: filled marker.
        h.append(Line2D([0], [0], marker='o', color='w',
                        markerfacecolor='grey', markeredgecolor='k',
                        markeredgewidth=0.4, markersize=7, label='Synchrotron'))
    return h


# ─── Scatter a single survey row ────────────────────────────────────────
def scatter_sample(ax, row, x, y, label=False, label_fontsize=5, **kw):
    """Plot one survey row as a scatter point. Returns the label Text
    object when label=True so the caller can collect them for adjustText.

    Lab samples render as hollow markers (group colour in the edge);
    Synchrotron samples render as filled markers (group colour interior,
    black edge). This replaces the legacy edge-only Lab/Synch distinction
    that was too subtle at publication scale (per Jenkins review comment,
    2026-05-16)."""
    group = row.get('Group', '')
    source = row.get('Source', '')
    group_colour = GROUP_COLOUR.get(group, 'grey')
    marker = GROUP_MARKER.get(group, 'o')
    fill = SOURCE_FILL.get(source, 'filled')  # default = filled (safe fallback)
    if fill == 'none':
        scatter_kw = dict(
            facecolors='none',
            edgecolors=group_colour,
            marker=marker,
            s=50, linewidths=1.2, zorder=5,
        )
    else:
        scatter_kw = dict(
            c=group_colour,
            marker=marker,
            edgecolors='k',
            s=50, linewidths=0.4, zorder=5,
        )
    scatter_kw.update(kw)
    ax.scatter(x, y, **scatter_kw)
    if label:
        return ax.text(x, y, row['Sample'],
                       fontsize=label_fontsize, va='center', zorder=10)
    return None


# ─── adjustText wrapper ─────────────────────────────────────────────────
# Defaults tuned on the JAC 29-sample figures. Override per-figure via kw.
ADJUST_DEFAULTS = dict(
    arrowprops=dict(arrowstyle='-', color='grey', lw=0.3),
    fontsize=5,
    expand=(1.8, 2.0),
    force_text=(0.6, 0.9),
    force_static=(0.3, 0.5),
    max_move=40,
    ensure_inside_axes=True,
)


def do_adjust(texts, ax, **kw):
    """adjustText with JAC-tuned defaults. Silently no-ops on empty input
    or adjustText exception (which can happen on pathological layouts)."""
    if not texts:
        return
    params = dict(ADJUST_DEFAULTS)
    params['ax'] = ax
    params.update(kw)
    try:
        adjust_text(texts, **params)
    except Exception:
        pass


# ─── Survey CSV loader ──────────────────────────────────────────────────
NUMERIC_COLS = ('Q_max', 'WA_D_median_A', 'WA_rms_strain',
                'PDF_pk1_r_A', 'PDF_pk1_FWHM_A')


def load_survey(csv_path):
    """Load survey CSV → list[dict] with numeric casting.

    Float columns (NaN on failure): Q_max, WA_D_median_A, WA_rms_strain,
                                    PDF_pk1_r_A, PDF_pk1_FWHM_A.
    Integer column (0 on failure):  WA_families.
    All other columns remain strings.
    """
    rows = []
    with open(csv_path, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            for k in NUMERIC_COLS:
                try:
                    row[k] = float(row[k])
                except (ValueError, TypeError, KeyError):
                    row[k] = np.nan
            try:
                row['WA_families'] = int(float(row['WA_families']))
            except (ValueError, TypeError, KeyError):
                row['WA_families'] = 0
            rows.append(row)
    return rows


def filter_survey(rows, *, group=None, source=None, groups=None, sources=None):
    """Subset helper. Any of group/source (single) or groups/sources (list)."""
    out = rows
    if group is not None:
        out = [r for r in out if r.get('Group') == group]
    if groups is not None:
        gs = set(groups)
        out = [r for r in out if r.get('Group') in gs]
    if source is not None:
        out = [r for r in out if r.get('Source') == source]
    if sources is not None:
        ss = set(sources)
        out = [r for r in out if r.get('Source') in ss]
    return out


__all__ = [
    'GROUP_COLOUR', 'GROUP_MARKER', 'SOURCE_EDGE', 'GROUP_ORDER',
    'GC', 'GM', 'SE',
    'legend_handles', 'scatter_sample', 'do_adjust',
    'ADJUST_DEFAULTS', 'NUMERIC_COLS',
    'load_survey', 'filter_survey',
]
