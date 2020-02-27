import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch


BASES = ["A", "T", "C", "G"]
MISMATCHES = ["C>T", "C>A", "C>G", "A>G", "A>C", "A>T"]
COLORS = ["#FFDF00", "#039AAC", "#DF1F36", "#038A3E", "#183D81", "#714883"]
CONTEXTS = [
    ["T_G", "T_A", "T_C", "T_T"],
    ["C_G", "C_A", "C_C", "C_T"],
    ["A_G", "A_A", "A_C", "A_T"],
    ["G_G", "G_A", "G_C", "G_T"],
]


def create_context_df():
    """ Create dataframe of all context and mismatch options
    with x and y positions in lego plot """
    positions = [[0, 0], [4, 0], [8, 0], [0, 4], [4, 4], [8, 4], [0, 8]]
    mismatches = MISMATCHES + ["NA"]  # NA is for context bars
    data = []
    for position, mismatch in zip(positions, mismatches):
        xo, yo = position
        x = np.tile(range(xo, xo + 4), 4)
        y = np.repeat(range(yo, yo + 4), 4)
        i = 0
        for row in CONTEXTS:
            for context in row:
                data.append([context, mismatch, x[i], y[i]])
                i += 1
    df = pd.DataFrame(data, columns=["context", "mismatch", "x", "y"])
    return df


def complement(sequence):
    """ Complement an ACTG sequence. If non-ACTG base,
    character is kept in output """
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complements.get(base, base) for base in sequence)


def create_mismatch_map(bases):
    """ Mismatches on opposites strands should be unified for
    creation of lego plot. This function creates a mapping between
    mismatch and unified mismatch """
    mismatch_map = {}
    for base1 in bases:
        for base2 in bases:
            if base1 == base2:
                continue
            else:
                mismatch = f"{base1}>{base2}"
                res = min([mismatch, complement(mismatch)])
                mismatch_map[mismatch] = res
    return mismatch_map


MISMATCH_MAP = create_mismatch_map(BASES)


def unify_mismatches(mismatch):
    """ For a given mismatch, return unified mismatch """
    return MISMATCH_MAP.get(mismatch)


def create_plot_df(counts_df):
    """ Expecting counts_df to have columns: context, mismatch, counts """
    counts_df = counts_df.copy()
    context_df = create_context_df()
    counts_df["mismatch"] = counts_df["mismatch"].map(unify_mismatches)
    counts_df = (
        counts_df.groupby(["context", "mismatch"])["counts"].sum().reset_index()
    )
    res = context_df.merge(counts_df, how="left").replace(np.nan, 0)
    return res


def create_lego_plot(df, title="", colors=COLORS, width=0.8):
    """
    df: Dataframe from create_plot_df
    """
    # Set scaling and viewing angle
    plt.figure(figsize=(14, 10))
    ax = plt.axes(projection="3d")
    ax.view_init(elev=60, azim=50)

    # Add colors to plotting dataframe
    df = df.copy()
    df.loc[df.mismatch == "NA", "counts"] = 0  # 0 counts on context bars
    colormap = dict(zip(MISMATCHES, colors))
    df["color"] = df["mismatch"].map(colormap).replace(np.nan, "white")

    # Create plot
    padding = (1 - width) / 2
    for gid, gdf in df.groupby("y"):
        z = np.zeros(len(gdf))
        dx = dy = [width] * len(gdf)
        shade = True if gid < 8 else False
        ax.bar3d(
            gdf.x + padding,
            gdf.y + padding,
            z,
            dx,
            dy,
            gdf["counts"],
            edgecolor="k",
            color=gdf.color,
            zsort="max",
            shade=shade,
        )
        ax.collections[gid].set_sort_zpos(gid)

    # Add legend
    legend_items = [
        Patch(facecolor=color, label=label)
        for color, label in zip(colors, MISMATCHES)
    ]
    ax.legend(handles=legend_items, loc=(0.8, 0.6))
    ax.text2D(0.1, 0.8, title, transform=ax.transAxes)

    # Add Context Info
    y = 8.5
    for row in CONTEXTS:
        x = 0.48
        for context in row:
            ax.text(
                x,
                y,
                0,
                context,
                ha="center",
                va="center",
                zdir=(0.3, -0.05, 0),
                fontsize=9,
                zorder=100,
            )
            x += 1.0
        y += 1.0

    # Configure axes
    ax.set(xticks=[], yticks=[], xticklabels=[], yticklabels=[])
    ax.set_ylim3d(bottom=0, top=12)
    ax.set_xlim3d(left=0, right=12)
    transparent = (1, 1, 1, 0)
    ax.w_xaxis.line.set_color(transparent)
    ax.w_yaxis.line.set_color(transparent)
    ax.xaxis.set_pane_color(transparent)
    ax.yaxis.set_pane_color(transparent)
    ax.zaxis.set_pane_color(transparent)
    return ax
