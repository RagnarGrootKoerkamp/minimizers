#!/usr/bin/env python3

# Plot for the minimizer density lower bound paper and thesis for small parameters.

import math
import gzip
import json
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import itertools
from functools import lru_cache
from itertools import chain, zip_longest
import pickle, random
import copy
from glob import glob
from pprint import pprint
from copy import deepcopy
import pickle as pck
from functools import cache
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from fractions import Fraction


import matplotlib.style
import matplotlib as mpl

from header import *

mpl.style.use("default")

title_font = 22
subplot_font = 14
plot_font = 18
legend_font = 16
col_font = subplot_font + 4

matplotlib.rc("font", size=subplot_font)

ylim_pad_factor = 0.05

# Same colours as mod-mini paper
scheme_to_color = {
    "Random minimizer": "orange",
    "Double decycling": "black",
    "Miniception": "#00dd00",
    "Mod-minimizer": "blue",
}
scheme_names = [
    "Random minimizer",
    "Double decycling",
    "Miniception",
    "Mod-minimizer",
]


# Uniform plotting style
def style(
    ax,
    w=None,
    sigma=None,
    k=1,
    title=None,
    ymax_pad_factor=ylim_pad_factor,
    miniception=None,
    ilp=False,
):
    ax.grid(False)
    ax.grid(True, axis="x", color="#ccc", linewidth=0.5)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    if w:
        ax.spines["bottom"].set_visible(False)
        ax.set_ylim(
            ymin=1 / w * (1 - ylim_pad_factor), ymax=(1 + ymax_pad_factor) * 2 / (w + 1)
        )

        rnge = 2 / (w + 1) - 1 / w
        dst = rnge / 6

        ticks = [1 / w]
        if sigma and ilp:
            extra = float(gp(sigma, w, k))
            ticks.append(extra)
        else:
            ticks.append(2 / (w + 1))

        tick_range = ticks[-1] - ticks[0]
        ticks.extend([ticks[0] + tick_range / 2])
        # extra = ragnar_WABI_LB(w, k)
        # if extra > 1/w + dst:
        #     ticks.append(extra)
        # ticks.append((ticks[0] + ticks[1]) / 2)
        if miniception:
            extra = miniception
            if all(abs(extra - t) > dst for t in ticks):
                ticks.append(extra)
        ax.set_yticks(ticks)
    else:
        ax.spines["bottom"].set_visible(True)
        ax.set_ylim(ymin=0, ymax=0.8)
        w = 2
        if sigma:
            g_sigma = float(g(w=2, k=1, sigma=sigma))
            ax.set_yticks([0, g_sigma / 2, g_sigma])
        else:
            ax.set_yticks([0, 1 / 4, 1 / 2, 2 / (w + 1)])

    if title:
        ax.set_title(title, loc="right", y=1.0, pad=-14, fontsize=subplot_font + 1)


def row_label(ax, label):
    ax.annotate(
        label,
        xy=(0, 0.5),
        xytext=(0, 0),
        xycoords=ax.yaxis.label,
        textcoords="offset points",
        size="large",
        ha="right",
        va="center",
        rotation=90,
    )


def col_label(ax, label):
    ax.annotate(
        label,
        xy=(0.5, 1),
        xytext=(0, 5),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="baseline",
    )


def hide_x_ticks(ax):
    ax.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False,
    )  # labels along the bottom edge are off


def hide_y_ticks(ax):
    ax.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        right=False,  # ticks along the top edge are off
        labelleft=False,
    )  # labels along the bottom edge are off


# lower bound style
# TODO: separate g and g'
lbs = {
    "linestyle": "solid",
    "linewidth": 1.5,
}

lower_bounds = [
    (lambda w, k, sigma: 1 / w, "$\\frac{1}{w}$ trivial bound", "black"),
    (
        lambda w, k, sigma: marcais(w, k),
        "$\\frac{1.5+1/2w+..}{w+k}$ Marçais et al.",
        "slategray",
    ),
    (
        lambda w, k, sigma: marcais_improved(w, k),
        "$\\frac{1.5}{w+k-0.5}$ Improved",
        "green",
    ),
    (
        lambda w, k, sigma: lb(w, k),
        "$\\frac{\\lceil\\frac{w+k}{w}\\rceil}{w+k}$ Simple",
        "royalblue",
    ),
    (lambda w, k, sigma: g(sigma, w, k), "$g_\\sigma$ Improved", "purple"),
    (lambda w, k, sigma: gp(sigma, w, k), "$g^\\prime_\\sigma$ Monotone", "red"),
]


def plot_lower_bounds(
    ax, *, w_range=None, k_range=None, w=None, k=None, sigma, ilp=False, subset=False
):
    # clone style
    my_lbs = lbs.copy()
    if subset:
        my_lower_bounds = [lower_bounds[-1], lower_bounds[0]]
        my_lbs["linewidth"] = 1.5
    else:
        my_lower_bounds = lower_bounds

    if k:
        for f, label, color in my_lower_bounds:
            if f is None:
                ax.plot([], [], label=" ", alpha=0)
            else:
                ax.plot(
                    w_range,
                    [f(w, k, sigma) for w in w_range],
                    label=label,
                    color=color,
                    **my_lbs,
                )
    else:
        for f, label, color in my_lower_bounds:
            if f is None:
                ax.plot([], [], label=" ", alpha=0)
            else:
                ax.plot(
                    k_range,
                    [f(w, k, sigma) for k in k_range],
                    label=label,
                    color=color,
                    **my_lbs,
                )

    if ilp:
        match_k = []
        match_dens = []
        nonmatch_k = []
        nonmatch_dens = []
        if k:
            for w in w_range:
                if (w, k, sigma) in fwd_wksigma_to_dens:
                    hasmatch = gp(sigma, w, k) == fwd_wksigma_to_dens[(w, k, sigma)]
                    if hasmatch:
                        match_k.append(w)
                        match_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])
                    else:
                        nonmatch_k.append(w)
                        nonmatch_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])
        else:
            for k in k_range:
                if (w, k, sigma) in fwd_wksigma_to_dens:
                    hasmatch = gp(sigma, w, k) == fwd_wksigma_to_dens[(w, k, sigma)]
                    if hasmatch:
                        match_k.append(k)
                        match_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])
                    else:
                        nonmatch_k.append(k)
                        nonmatch_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])
        # Plot transparent circle
        ax.scatter(
            nonmatch_k,
            nonmatch_dens,
            label="Optimal solution",
            color="black",
            s=30,
            facecolors="none",
            lw=1.3,
            zorder=100,
        )
        ax.scatter(
            match_k,
            match_dens,
            label="Optimum matches bound",
            color="black",
            s=30,
            zorder=100,
        )
        ax.set_ylim(ymax=(1 + ylim_pad_factor) * match_dens[0])


# MAIN PLOT

ax_width = 3
ax_height = 3
# ax_width, ax_height = ax_height, ax_width
fig, axs = plt.subplots(nrows=ax_height, ncols=ax_width)
fig.set_figheight(6)
fig.set_figwidth(12)
fig.subplots_adjust(wspace=0.25)

w_range = list(range(3, 6))
sigma = 2
for idx, w in enumerate(w_range):
    k_range = list(range(1, 15))
    ax = axs[idx, 0]
    style(ax, w, sigma, title=f"$w={w}$  ", ilp=True)
    plot_lower_bounds(ax, k_range=k_range, w=w, sigma=sigma, ilp=True)

    if w != w_range[-1]:
        hide_x_ticks(ax)
    else:
        ax.set_xlabel("$k$")

sigma_range = list(range(2, 5))
w = 2
for idx, sigma in enumerate(sigma_range):
    k_range = list(range(1, 15))
    ax = axs[idx, 1]
    style(ax, w, sigma, title=f"$\\sigma={sigma}$  ", ilp=True)
    plot_lower_bounds(ax, k_range=k_range, w=w, sigma=sigma, ilp=True)

    if sigma != sigma_range[-1]:
        hide_x_ticks(ax)
    else:
        ax.set_xlabel("$k$")

k = 1
for idx, sigma in enumerate(sigma_range):
    w_range = list(range(2, 13))
    ax = axs[idx, 2]
    style(ax, None, sigma, title=f"$\\sigma={sigma}$  ", ilp=True)
    plot_lower_bounds(ax, w_range=w_range, k=k, sigma=sigma, ilp=True)

    if sigma != sigma_range[-1]:
        hide_x_ticks(ax)
    else:
        ax.set_xlabel("$w$")


col_label(axs[0][0], f"$\\sigma=2$")
col_label(axs[0][1], f"$w=2$")
col_label(axs[0][2], f"$k=1$")

handles, labels = ax.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    bbox_to_anchor=(0.5, -0.24),
    loc="lower center",
    prop={"size": legend_font},
    ncols=3,
)
fig.supylabel("Density", x=0.04, fontsize=plot_font)
# fig.suptitle(
#     "Density lower bounds and optimal solutions for forward schemes",
#     fontsize=title_font - 1,
# )
fig.savefig(f"thesis/lower-bound.svg", bbox_inches="tight")
# fig.show()
# plt.close(fig)
