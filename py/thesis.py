#!/usr/bin/env python3

# Plots for my thesis chapters on minimizers: https://curiouscoding.nl/posts/minimizers

import header as h
from header import plot
from sympy import Symbol, Lambda, Max, Min

sigma = 4
w = 15
k = Symbol("k")
ks = None
t = 4

n = 2000000
h.gen(n, sigma)

# Plot 0: t value
# k =
# w = 24
# k = 60
# ts = range(1, k + 1)
# fs = [
#     ("Random", {"sampling": 0, "lc": "orange", "label": "Random"}),
#     ("Random", {"mod": 1, "sampling": 0, "lc": "blue", "label": "Random mod-mini"}),
# ]
# plot(
#     "thesis/0-mod-t",
#     sigma,
#     w,
#     fs,
#     plot_t=True,
#     k=k,
#     ymin=0.0380,
#     ymax=0.115,
#     height=2.8,
#     ts=ts,
# )

# Plot 1: Random and lexicographic minimizers
ks = range(1, 9)
fs = [
    ("Lex", {"label": "Lex"}),
    ("Random", {"lc": "orange", "label": "Random"}),
    ("Alternating", {"label": "Alternating"}),
    ("AntiLex", {"label": "AntiLex"}),
    ("ABB", {"label": "ABB"}),
    ("ABB2", {"lc": "teal", "label": "ABB+"}),
    ("RcAntiLex", {"label": "RcAntiLex"}),
    ("RcAntiLexMax", {"label": "RcAntiLexMax"}),
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
plot("test", sigma, w, fs, ymin=0.09, ymax=0.14, ks=ks)

# # Plot 2: Pure schemes
# ks = range(1, 60)
# fs = [
#     ("Random", {"lc": "orange", "label": "Random"}),
#     ("ABB2", {"lc": "teal", "label": "ABB+"}),
#     ("Decycling", {"lc": "grey", "label": "Decycling"}),
#     ("DoubleDecycling", {"lc": "black", "label": "Double Decycling"}),
#     ("OpenClosed", {"r": t, "open": 1, "lc": "#bb0066", "label": "O"}),
#     (
#         "OpenClosed",
#         {"r": t, "closed": 1, "lc": "#00dd00", "label": "C=Miniception"},
#     ),
#     (
#         "OpenClosed",
#         {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "label": "OC"},
#     ),
# ]
# plot("thesis/2-ext", sigma, w, fs, ymin=0.0485, add=0, ncols=4, trivial=False, ks=ks)

# # Plot 3: Mod schemes
# ks = range(1, 83)


# def fs_mod(t, sus=False):
#     l = [
#         ("Random", {"lc": "orange", "label": "Random", "thin": True}),
#         ("Random", {"mod": 1, "r": t, "lc": "orange", "label": "Random-mod"}),
#         ("ABB2", {"lc": "teal", "label": "ABB+", "thin": True}),
#         ("ABB2", {"mod": 1, "r": t, "lc": "teal", "label": "ABB+-mod"}),
#         (
#             "DoubleDecycling",
#             {"lc": "black", "label": "Double Decycling", "thin": True},
#         ),
#         (
#             "DoubleDecycling",
#             {"mod": 1, "r": t, "lc": "black", "label": "Double Decycling-mod"},
#         ),
#         (
#             "OpenClosed",
#             {"r": t, "open": 1, "closed": 1, "lc": "#6600bb", "thin": True},
#         ),
#         (
#             "OpenClosed",
#             {
#                 "mod": 1,
#                 "r": t,
#                 "open": 1,
#                 "closed": 1,
#                 "lc": "#6600bb",
#                 "label": "OC-mod",
#             },
#         ),
#     ]
#     if sus:
#         l += [
#             (
#                 "SusAntiLex",
#                 {"lc": "#0000ff", "label": "Sus-anchor AntiLex", "thin": True},
#             ),
#             (
#                 "SusAntiLex",
#                 {
#                     "mod": 1,
#                     "r": 1,
#                     "lc": "#0000ff",
#                     "label": "Sus-anchor AntiLex mod",
#                 },
#             ),
#             # (
#             #     "Threshold",
#             #     {
#             #         "r": t,
#             #         "t": 192,
#             #         "open": 0,
#             #         "h": 1,
#             #         "loose": 1,
#             #         "mod": 1,
#             #         "lc": "teal",
#             #         "label": "Threshold",
#             #     },
#             # ),
#         ]
#     return l


# plot("thesis/3-mod", sigma, w, fs_mod(t), add=0, ncols=3, trivial=True, ks=ks)

# # Plot 3: k=1 anchors
# k = 1
# ws = range(2, 21)
# fs = [
#     ("BdAnchor", {"r": 0, "lc": "#00ff00", "label": "Bd-anchor r=0"}),
#     ("BdAnchor", {"r": 1, "lc": "#00dd00", "label": "Bd-anchor r=1"}),
#     ("BdAnchor", {"r": 2, "lc": "#00bb00", "label": "Bd-anchor r=2"}),
#     ("BdAnchor", {"r": 3, "lc": "#009900", "label": "Bd-anchor r=3"}),
#     ("SusLex", {"lc": "#a00080", "label": "Sus-anchor Lex"}),
#     ("SusABB", {"lc": "#8010b0", "label": "Sus-anchor ABB"}),
#     ("SusAntiLex", {"lc": "#0000ff", "label": "Sus-anchor AntiLex"}),
# ]
# plot(
#     "thesis/4-selection",
#     sigma,
#     w,
#     fs,
#     plot_w=True,
#     tight=True,
#     ymin=0,
#     ymax=0.85,
#     ncols=3,
#     add=1,
#     bold_last=True,
#     ws=ws,
# )

# # Plot 4: selection for alphabet 2
# sigma = 2
# h.gen(n, sigma)
# plot(
#     "thesis/5-selection-s2",
#     sigma,
#     w,
#     fs,
#     plot_w=True,
#     tight=True,
#     ymin=0,
#     ymax=0.85,
#     ncols=3,
#     add=1,
#     bold_last=True,
#     ws=ws,
# )

# # Plot 5: sampling with selection for alphabet 4
# sigma = 4
# w = 24
# ks = range(1, 40)
# h.gen(n, sigma)
# plot(
#     "thesis/6-small-k",
#     sigma,
#     w,
#     fs_mod(t, sus=True),
#     add=0,
#     ncols=3,
#     trivial=True,
#     ks=ks,
# )

# # Plot 6: sampling with selection for alphabet 256
# # sigma = 256
# # t = 1
# # h.gen(n, sigma)
# # plot(
# #     "thesis/7-large-alphabet",
# #     sigma,
# #     w,
# #     fs_mod(t, sus=True),
# #     add=0,
# #     ncols=3,
# #     trivial=True,
# #     ks=ks,
# # )
