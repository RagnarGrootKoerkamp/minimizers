#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import log
import matplotlib as mpl
import header as h
from header import plt
from sympy import Symbol, Lambda, Max, Min
from matplotlib.ticker import MaxNLocator


## PLOTTING


def plot(
    name,
    tps,
    add=0,
    bold_last=False,
    plot_w=False,
    plot_t=False,
    ymin=0.0405,
    ymax=0.082,
    rol=True,
    df=False,
    ncols=4,
    k=None,
    height=4.8,
    **kwargs,
):
    global w
    data = []
    # if plot_w:
    #     title = f"Densities and lower bounds (σ={sigma})"
    # else:
    #     title = f"Densities and lower bounds (σ={sigma}, w={w})"
    ax = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(6.4, height)
    h.style(ax, w, ks, plot_w=plot_w, plot_t=plot_t, df=df)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.grid(True, axis="y", color="#ccc", linewidth=0.5)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    if plot_t:
        wks = [(w, k) for _ in ts]
        xs = ts
    elif plot_w:
        wks = [(w, 1) for w in ws]
        xs = ws
    else:
        wks = [(w, k) for k in ks]
        xs = ks

    for i, tp_args in enumerate(tps):
        if isinstance(tp_args, tuple):
            (tp, args) = tp_args
        else:
            tp = tp_args
            args = {}

        last = 1.5 if bold_last and i == len(tps) - 1 else 0

        print(f"{tp}: {args}", flush=True)
        ds = []
        if plot_t:
            for t in ts:
                my_args = {k: v for (k, v) in args.items()}
                if "lc" in my_args:
                    my_args.pop("lc")
                if "ls" in my_args:
                    my_args.pop("ls")
                if "o" in args:
                    my_args["offset"] = args["o"](k)
                if "k0" in args and not isinstance(args["k0"], int):
                    my_args["k0"] = my_args["k0"](k)
                my_args["sampling"] = t
                d = h.density(tp, w, k, sigma, **my_args)
                if df:
                    d *= w + 1
                my_args["k"] = k
                my_args["w"] = w
                my_args["d"] = d
                my_args["tp"] = tp
                my_args["sigma"] = sigma

                data.append(my_args)
                ds.append(d)
        else:
            for w, k in wks:
                my_args = {k: v for (k, v) in args.items()}
                if "lc" in my_args:
                    my_args.pop("lc")
                if "ls" in my_args:
                    my_args.pop("ls")
                if "o" in args:
                    my_args["offset"] = args["o"](k)
                if "k0" in args and not isinstance(args["k0"], int):
                    my_args["k0"] = my_args["k0"](k)
                d = h.density(tp, w, k, sigma, **my_args)
                if df:
                    d *= w + 1
                my_args["k"] = k
                my_args["w"] = w
                my_args["d"] = d
                my_args["tp"] = tp
                my_args["sigma"] = sigma

                data.append(my_args)
                ds.append(d)
        print(ds[0])
        lc = args.get("lc", None)
        ls = args.get("ls", None)
        label = str(args.get("label", tp))
        if label is None:
            args.pop("label")
            if "lc" in args:
                args.pop("lc")
            if "ls" in args:
                args.pop("ls")
            label = f"{tp}: {args}"
        if not args.get("ao", False) and not args.get("aot", False):
            ls = "solid"
        if args.get("ao", False) and not args.get("aot", False):
            label += " ak"
            if not ls:
                ls = "dotted"
        if not args.get("ao", False) and args.get("aot", False):
            label += " at"
            if not ls:
                ls = "dashed"
        if args.get("ao", False) and args.get("aot", False):
            label += " at,ak"
            if not ls:
                ls = "dashdot"
        lw = 1 + last
        if args.get("modulo", False):
            lw -= 0.5
        alpha = 1
        if args.get("thin", False):
            alpha = 0.4
            lw = 0.5
            last = -1
            label = None

        (line,) = plt.plot(
            xs,
            ds,
            label=None,
            linestyle=ls,
            color=lc,
            marker="o",
            markersize=2 + last,
            lw=0,
            alpha=alpha,
        )

        # Rolling minimum of ds
        dm = 1
        dms = []
        for d in ds:
            if rol and not plot_t:
                dm = min(dm, d)
            else:
                dm = d
            dms.append(dm)
        plt.plot(
            xs,
            dms,
            label=None,
            color=line.get_color(),
            markersize=0,
            lw=lw,
            alpha=alpha,
        )
        # Add unified legend item.
        plt.plot(
            [],
            [],
            label=label,
            color=line.get_color(),
            marker="o",
            markersize=2 + last,
            lw=lw,
        )

    for _ in range(add):
        plt.plot([], [], label=" ", alpha=0)

    h.plot_lower_bounds(sigma, xs, wks, df=df, **kwargs)

    if plot_w:
        loc = "upper center"
    else:
        loc = "lower center"
    plt.legend(
        loc=loc, bbox_to_anchor=(0, 0.03, 1, 1), ncols=ncols, mode="expand", fontsize=9
    )
    plt.savefig(f"{name}.png", bbox_inches="tight", dpi=400)
    plt.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close()
    return data


# THESIS PLOTS

sigma = 4
w = 24
k = Symbol("k")
ks = None
t = 4

n = 10000000
h.gen(n, sigma)

# Plot 0: t value
# k =
w = 24
k = 60
ts = range(1, k + 1)
fs = [
    ("Random", {"sampling": 0, "lc": "orange", "label": "Random"}),
    ("Random", {"mod": 1, "sampling": 0, "lc": "blue", "label": "Random mod-mini"}),
]
plot("thesis/0-mod-t", fs, plot_t=True, k=k, ymin=0.0380, ymax=0.115, height=2.8)

# Plot 1: Random and lexicographic minimizers
ks = range(1, 32)
fs = [
    ("Lex", {"label": "Lex"}),
    ("Random", {"lc": "orange", "label": "Random"}),
    ("Alternating", {"label": "Alternating"}),
    ("AntiLex", {"label": "AntiLex"}),
    ("ABB", {"label": "ABB"}),
    ("ABB2", {"lc": "teal", "label": "ABB+"}),
    # (
    #     "Threshold",
    #     {
    #         "r": t,
    #         "t": 3,
    #         "open": 0,
    #         "h": 1,
    #         "loose": 1,
    #         "mod": 1,
    #         "lc": "teal",
    #         "label": "Threshold",
    #     },
    # ),
]
plot("thesis/1-lex", fs, ymax=0.09)

# Plot 2: Pure schemes
ks = range(1, 60)
fs = [
    ("Random", {"lc": "orange", "label": "Random"}),
    ("ABB2", {"lc": "teal", "label": "ABB+"}),
    ("Decycling", {"lc": "grey", "label": "Decycling"}),
    ("DoubleDecycling", {"lc": "black", "label": "Double Decycling"}),
    ("OpenClosed", {"r": t, "open": 1, "lc": "#bb0066", "label": "O"}),
    (
        "OpenClosed",
        {"r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
    ),
    (
        "OpenClosed",
        {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "label": "OC"},
    ),
]
plot("thesis/2-ext", fs, ymin=0.0485, add=0, ncols=4, trivial=False)

# Plot 3: Mod schemes
ks = range(1, 83)


def fs_mod(t, sus=False):
    l = [
        ("Random", {"lc": "orange", "label": "Random", "thin": True}),
        ("Random", {"mod": 1, "r": t, "lc": "orange", "label": "Random-mod"}),
        ("ABB2", {"lc": "teal", "label": "ABB+", "thin": True}),
        ("ABB2", {"mod": 1, "r": t, "lc": "teal", "label": "ABB+-mod"}),
        (
            "DoubleDecycling",
            {"lc": "black", "label": "Double Decycling", "thin": True},
        ),
        (
            "DoubleDecycling",
            {"mod": 1, "r": t, "lc": "black", "label": "Double Decycling-mod"},
        ),
        (
            "OpenClosed",
            {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "thin": True},
        ),
        (
            "OpenClosed",
            {
                "mod": 1,
                "r": t,
                "open": 1,
                "closed": 1,
                "lc": "#6600bb",
                "label": "OC-mod",
            },
        ),
    ]
    if sus:
        l += [
            (
                "SusAntiLex",
                {"lc": "#0000ff", "label": "Sus-anchor AntiLex", "thin": True},
            ),
            (
                "SusAntiLex",
                {
                    "mod": 1,
                    "r": 1,
                    "lc": "#0000ff",
                    "label": "Sus-anchor AntiLex mod",
                },
            ),
            # (
            #     "Threshold",
            #     {
            #         "r": t,
            #         "t": 192,
            #         "open": 0,
            #         "h": 1,
            #         "loose": 1,
            #         "mod": 1,
            #         "lc": "teal",
            #         "label": "Threshold",
            #     },
            # ),
        ]
    return l


plot("thesis/3-mod", fs_mod(t), add=0, ncols=3, trivial=True)

# Plot 3: k=1 anchors
k = 1
ws = range(2, 21)
fs = [
    ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
    ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
    ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
    ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
    ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
    ("SusABB", {"lc": "#8010b0", "label": "Sus-anchor ABB"}),
    ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
]
plot(
    "thesis/4-selection",
    fs,
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
)

# Plot 4: selection for alphabet 2
sigma = 2
h.gen(n, sigma)
plot(
    "thesis/5-selection-s2",
    fs,
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
)

# Plot 5: sampling with selection for alphabet 4
sigma = 4
w = 24
ks = range(1, 40)
h.gen(n, sigma)
plot("thesis/6-small-k", fs_mod(t, sus=True), add=0, ncols=3, trivial=True)

# Plot 6: sampling with selection for alphabet 256
# sigma = 256
# t = 1
# h.gen(n, sigma)
# plot("thesis/7-large-alphabet", fs_mod(t, sus=True), add=0, ncols=3, trivial=True)
