import minimizers.minimizers
import matplotlib.pyplot as plt
import sympy
from sympy import totient, isprime, divisors
from sympy.functions.combinatorial.numbers import mobius
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from functools import cache
import fractions
import json

## PLOT STYLING


def style(ax, w, ks, title=None, plot_w=False, plot_t=False, df=False):

    ax.grid(False)
    ax.grid(True, axis="x", color="#ccc", linewidth=0.5)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    if not plot_w:
        ax.spines["bottom"].set_visible(False)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    df = w + 1 if df else 1

    if not plot_w and not plot_t:
        ylim_pad_factor = 0.02

        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
        ymin = 1 / w * (1 - ylim_pad_factor) * df
        ymax = (1 + ylim_pad_factor) * 2 / (w + 1) * df
        ax.set_ylim(ymin=ymin, ymax=ymax)
        ax.set_yticks([1 / w * df, 1.5 / (w + 0.5) * df, 2 / (w + 1) * df])
        ax.set_xlim(xmin=ks[0] / 2, xmax=ks[-1] + ks[0] / 2)
    if plot_w:
        ax.set_xlabel("w")
    elif plot_t:
        ax.set_xlabel("t")
    else:
        ax.set_xlabel("k")
    ax.set_ylabel("Density")
    if title:
        ax.set_title(title)


def plot_lower_bounds(
    sigma,
    xs,
    wks,
    df=False,
    loose=False,
    tight=True,
    ctd=False,
    marcais=False,
    trivial=True,
    mini=False,
    ilp=False,
):
    # Fix for density factor
    f = lambda w: w + 1 if df else 1
    if mini:
        plt.plot(
            xs,
            [1 / sigma**k * f(w) for (w, k) in wks],
            color="red",
            linewidth=1,
            ls="--",
            label="1/σ^k",
        )
    if ctd:
        plt.plot(
            xs,
            [lb_continuation(w, k) * f(w) for (w, k) in wks],
            color="black",
            linewidth=1,
            label="continuation",
        )
    if tight:
        plt.plot(
            xs,
            [gp(sigma, w, k) * f(w) for (w, k) in wks],
            color="red",
            linewidth=1.5,
            label="New l.b. (g')",
        )
    if loose:
        plt.plot(
            xs,
            [lbp(w, k) * f(w) for (w, k) in wks],
            color="royalblue",
            linewidth=1.5,
            label="⌈(w+k)/w⌉/(w+k)",
        )
    if marcais:
        plt.plot(
            xs,
            [marcais_improved(w, k) * f(w) for (w, k) in wks],
            color="grey",
            linewidth=1.5,
            label="1.5/(w+k-0.5)",
        )
    if trivial:
        plt.plot(
            xs,
            [1 / w * f(w) for (w, k) in wks],
            color="black",
            linewidth=1.5,
            label="1/w",
        )
    if ilp:
        match_k = []
        match_dens = []
        nonmatch_k = []
        nonmatch_dens = []
        for (w, k), x in zip(wks, xs):
            if (w, k, sigma) in fwd_wksigma_to_dens:
                hasmatch = gp(sigma, w, k) == fwd_wksigma_to_dens[(w, k, sigma)]
                if hasmatch:
                    match_k.append(x)
                    match_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])
                else:
                    nonmatch_k.append(x)
                    nonmatch_dens.append(fwd_wksigma_to_dens[(w, k, sigma)])

        # Plot transparent circle
        ax = plt.gca()
        ax.scatter(
            nonmatch_k,
            nonmatch_dens,
            # label="Optimal solution",
            label=None,
            color="black",
            s=30,
            facecolors="none",
            lw=1.3,
            zorder=100,
        )
        ax.scatter(
            match_k,
            match_dens,
            # label="Optimum matches bound",
            label=None,
            color="black",
            s=30,
            zorder=100,
        )
        # ax.set_ylim(ymax=(1 + ylim_pad_factor) * match_dens[0])


## LOWER BOUNDS


def marcais(w, k):
    return (1.5 + max(0, (k - w) // w) + 1 / (2 * w)) / (w + k)


def marcais_improved(w, k):
    return 1.5 / (w + k - 0.5)


# Simplified lb
def lb(w, k):
    return 1 / (w + k) * ((w + k + (w - 1)) // w)


# Suffix-max version (lb')
def lbp(w, k):
    k2 = 1 + (k - 1 + w - 1) // w * w
    return max(lb(w, k), lb(w, k2))


# Continutation
def lb_continuation(w, k):
    return 1 / (w + k) * ((w + k + (w - 1)) / w)


# Helper
def aperiodic_necklaces(n, sigma=2):
    return int((sum(mobius(d) * sigma ** (n // d) for d in divisors(n))) // n)


# Precise lower bound
def g(sigma, w, k):
    s = 0
    for p in divisors(w + k):
        s += aperiodic_necklaces(p, sigma) * ((p + w - 1) // w)
    return fractions.Fraction(s, (sigma ** (w + k)))


# g'
def gp(sigma, w, k):
    k2 = 1 + (k - 1 + w - 1) // w * w
    return max(g(sigma, w, k), g(sigma, w, k2))


## Cache rust calls


from diskcache import Cache

diskcache = Cache("cache")


@cache
def gen_inner(text_len, sigma):
    return minimizers.generate_random_string(text_len, sigma)


_text = []


def gen(text_len, sigma):
    global _text
    _text = gen_inner(text_len, sigma)


@diskcache.memoize(tag="density")
def density_inner(tp, text_len, w, k, sigma, **args):
    return minimizers.density(tp, _text, w, k, sigma, **args)


def density(tp, w, k, sigma, **args):
    return density_inner(tp, len(_text), w, k, sigma, **args)


fwd_wksigma_to_dens = {}
fwd = json.load(open("fwd.json"))
for w, k, sigma, d1, d2 in fwd:
    fwd_wksigma_to_dens[(w, k, sigma)] = fractions.Fraction(d1, d2)

# The main plotting function.


def plot(
    name,
    sigma,
    w,
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
    ks=None,
    ws=None,
    ts=None,
    height=4.8,
    title=None,
    **kwargs,
):
    data = []
    ax = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(6.4, height)
    style(ax, w, ks, plot_w=plot_w, plot_t=plot_t, df=df, title=title)
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
                if "label" in my_args:
                    my_args.pop("label")
                if "thin" in my_args:
                    my_args.pop("thin")
                my_args["sampling"] = t
                print("sigma", sigma)
                d = density(tp, w, k, sigma, **my_args)
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
                if "label" in my_args:
                    my_args.pop("label")
                if "thin" in my_args:
                    my_args.pop("thin")
                d = density(tp, w, k, sigma, **my_args)
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

    plot_lower_bounds(sigma, xs, wks, df=df, **kwargs)

    if plot_w:
        loc = "upper center"
    else:
        loc = "lower center"
    plt.legend(
        loc=loc, bbox_to_anchor=(0, 0.03, 1, 1), ncols=ncols, mode="expand", fontsize=9
    )
    # plt.savefig(f"{name}.png", bbox_inches="tight", dpi=400)
    plt.savefig(f"{name}.svg", bbox_inches="tight")
    plt.close()
    return data
