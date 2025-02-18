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


def style(ax, w, ks, title=None, plot_w=False, df=False):

    ax.grid(False)
    ax.grid(True, axis="x", color="#ccc", linewidth=0.5)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    ax.spines["bottom"].set_visible(False)

    df = w + 1 if df else 1

    if not plot_w:
        ylim_pad_factor = 0.02

        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
        ymin = 1 / w * (1 - ylim_pad_factor) * df
        ymax = (1 + ylim_pad_factor) * 2 / (w + 1) * df
        ax.set_ylim(ymin=ymin, ymax=ymax)
        ax.set_yticks([1 / w * df, 1.5 / (w + 0.5) * df, 2 / (w + 1) * df])
        ax.set_xlim(xmin=ks[0] / 2, xmax=ks[-1] + ks[0] / 2)
    if plot_w:
        ax.set_xlabel("w")
    else:
        ax.set_xlabel("k")
    ax.set_ylabel("Density")
    if title:
        ax.set_title(title)


def plot_lower_bounds(
    sigma, xs, wks, loose=True, tight=False, ctd=False, marcais=False, df=False
):
    f = lambda w: w + 1 if df else 1
    if tight:
        plt.plot(
            xs,
            [gp(sigma, w, k) * f(w) for (w, k) in wks],
            color="red",
            linewidth=1.5,
            label="g'",
        )
    if loose:
        plt.plot(
            xs,
            [lbp(w, k) * f(w) for (w, k) in wks],
            color="red",
            linewidth=1.5,
            label="⌈(w+k)/w⌉/(w+k)",
        )
    if ctd:
        plt.plot(
            xs,
            [lb_continuation(w, k) * f(w) for (w, k) in wks],
            color="black",
            linewidth=1,
            label="continuation",
        )
    plt.plot(
        xs, [1 / w * f(w) for (w, k) in wks], color="black", linewidth=1.5, label="1/w"
    )


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
