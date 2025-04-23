#!/usr/bin/env python3

# Plots for DSB 2025 slides on minimizers: https://curiouscoding.nl/slides/minimizers-dsb25/

import header as h
from header import plot
from sympy import Symbol, Lambda, Max, Min

sigma = 4
w = 24
k = Symbol("k")
ks = None
t = 4

n = 10000000
h.gen(n, sigma)

# Plot 2: Pure schemes
ks = range(1, 83)
fs = [
    ("Random", {"lc": "orange", "label": "Random"}),
    # ("ABB2", {"lc": "teal", "label": "ABB+"}),
    ("Decycling", {"lc": "grey", "label": "Decycling"}),
    ("DoubleDecycling", {"lc": "black", "label": "Double Decycling"}),
    # ("OpenClosed", {"r": t, "open": 1, "lc": "#bb0066", "label": "O"}),
    (
        "OpenClosed",
        {"r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
    ),
    # (
    #     "OpenClosed",
    #     {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "label": "OC"},
    # ),
]
plot(
    "dsb-slides/1-before",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=0,
    ncols=4,
    trivial=True,
    tight=False,
    marcais=True,
    mini=True,
    ks=ks,
)

# Add mod
fs = [
    # ("Random", {"lc": "orange", "label": "Random"}),
    ("Random", {"lc": "orange", "label": "Random", "thin": True}),
    ("Random", {"mod": 1, "r": t, "lc": "orange", "label": "Random-mod"}),
    # ("ABB2", {"lc": "teal", "label": "ABB+"}),
    ("Decycling", {"lc": "grey", "label": "Decycling"}),
    ("DoubleDecycling", {"lc": "black", "label": "Double Decycling"}),
    # ("OpenClosed", {"r": t, "open": 1, "lc": "#bb0066", "label": "O"}),
    (
        "OpenClosed",
        {"r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
    ),
    # (
    #     "OpenClosed",
    #     {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "label": "OC"},
    # ),
]
plot(
    "dsb-slides/2-mod",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=0,
    ncols=4,
    trivial=True,
    tight=False,
    marcais=True,
    mini=True,
    ks=ks,
)
plot(
    "dsb-slides/3-lb",
    sigma,
    w,
    fs,
    ymin=0.040,
    add=0,
    ncols=4,
    trivial=True,
    tight=True,
    marcais=True,
    mini=True,
    ks=ks,
)


# Plot 3: Mod schemes


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


plot(
    "dsb-slides/4-full",
    sigma,
    w,
    fs_mod(t),
    add=1,
    ncols=3,
    trivial=True,
    mini=True,
    ymin=0.04,
    marcais=True,
    ks=ks,
)

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
    "dsb-slides/10-bd-anchor",
    sigma,
    w,
    fs[0:4],
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=False,
    ws=ws,
)
plot(
    "dsb-slides/11-sus",
    sigma,
    w,
    fs[0:5],
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
    ws=ws,
)
plot(
    "dsb-slides/12-abb",
    sigma,
    w,
    fs[0:6],
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
    ws=ws,
)
plot(
    "dsb-slides/13-asus",
    sigma,
    w,
    fs[0:7],
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
    ws=ws,
)


# Plot 4: selection for alphabet 2
sigma = 2
h.gen(n, sigma)
plot(
    "dsb-slides/14-asus-s2",
    sigma,
    w,
    fs,
    plot_w=True,
    tight=True,
    ymin=0,
    ymax=0.85,
    ncols=3,
    add=1,
    bold_last=True,
    ws=ws,
)
