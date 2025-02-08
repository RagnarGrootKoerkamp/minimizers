#!/usr/bin/env python3

import pickle as pck
import itertools
import sys
from math import log

w = 3
k = 1
s = 3
if len(sys.argv) == 4:
    k, s, w = map(int, sys.argv[1:])

l = w + k - 1
ll = l + 1
with open(f"fwd/sols/w{w}-k{k}-s{s}.pck", "rb") as pck_in:
    scheme = pck.load(pck_in)[(w, k, s)][0]


def find(g, kmer):
    for i, gk in enumerate(g):
        eq = True
        for x, y in zip(kmer, gk):
            if x == y or x == "?" or y == "?":
                continue
            eq = False
        if eq:
            return i
    return None


def greedy(scheme, g=[]):
    dg = {}
    for i, kmer in enumerate(g):
        dg[kmer] = i

    def f(lmer):
        best = (10000, 0)
        for idx in range(w):
            kmer = lmer[idx : idx + k]
            kmer = "".join(str(x) for x in kmer)
            x = find(g, kmer)
            if x is not None:
                best = min(best, (x, idx))
        if best[0] < 10000:
            return best[1]
        return scheme[lmer]

    return f


SUF = [0, -0.5]


def asus(lmer):
    global SUF
    suffixes = sorted(
        [
            (
                [ord(lmer[i]) - ord("0")]
                + [-(ord(x) - ord("0")) for x in lmer[i + 1 :]]
                + SUF,
                i,
            )
            for i in range(len(lmer))
        ]
    )
    return suffixes[0][1]


def merge(pre, suf):
    P = 0.5
    S = 0.5
    l = []
    for i in range(max(len(pre), len(suf))):
        l.append(+(ord(pre[len(pre) - 1 - i]) - ord("0") if i < len(pre) else P))
        l.append(-(ord(suf[i]) - ord("0") if i < len(suf) else S))
    # print("merge", pre, suf, l)
    return l


def lr(lmer):
    global SUF
    suffixes = sorted(
        [(merge(lmer[: i + 1], lmer[i + 1 :]), i) for i in range(len(lmer))]
    )
    # for x in suffixes: print(x)
    return suffixes[0][1]


# Optimal for k=1,s=2,w<=9
def lr2(lmer, all=False):
    l = len(lmer)
    # Find position of last 0 in the longest 000..0011..1111 run.
    best = (0, 0, -100, -100)
    bests = [l - 1]
    for i in range(l - 1):
        if lmer[i] == "0" and lmer[i + 1] == "1":
            c0 = 0
            d0 = 0
            for j in range(i, -1, -1):
                if lmer[j] == "0":
                    c0 += 1
                else:
                    break
            for j in range(j, -1, -1):
                if lmer[j] == "1":
                    d0 += 1
                else:
                    break
            c1 = 0
            d1 = 0
            for j in range(i + 1, l):
                if lmer[j] == "1":
                    c1 += 1
                else:
                    break
            for j in range(j, l):
                if lmer[j] == "0":
                    d1 += 1
                else:
                    break
            # best = max(best, (c0 + c1, c1, d0 + d1, +d1, -i, i))
            val = (c0 + c1, c1, 0, 0)
            if val > best:
                best = val
                bests = [i]
            elif val == best:
                bests.append(i)
            # print(lmer, i, best)

    if len(bests) == 1 or False:
        if all:
            return (bests[0], bests)
        return bests[0]

    x = bests[-1]
    for i in range(len(bests) - 1):
        if (bests[i + 1] - bests[i]) % 2 == 0:
            x = bests[i]
            break

    if all:
        return (x, bests)
    return x

    # Left: find first 1; pos before that
    # f1 = lmer.index("1")
    # if f1 < bests[0]:
    #     fl = f1 - 1
    # else:
    #     fl = -2
    # x = [fl] + bests + [ll - 1 - int(lmer[-1] == "0")]
    # print(x)
    # b = (-1, -1, -1)
    # for i in range(1, len(x) - 1):
    #     b = max(b, (x[i + 1] - x[i], x[i] - x[i - 1], x[i]))
    # print(b)
    # if all:
    #     return (b[2], x)
    return b[2]


lr = lr2

# lr("1010100")
# exit(0)


def density(scheme, ll=None):
    if isinstance(scheme, dict):
        f = lambda key: scheme[key]
    else:
        f = scheme
    if ll is None:
        ll = w + k
    # Iterate over all ll'mers
    count = 0
    for llmer in itertools.product(range(s), repeat=ll):
        llmer2 = llmer + llmer
        poss = set()
        for i in range(ll):
            lmer = llmer2[i : i + l]
            poss.add((i + f("".join(str(x) for x in lmer))) % ll)
        # print("".join(str(x) for x in llmer), poss)
        count += len(poss)
    total_lmers = s**ll * ll
    density = count / total_lmers
    return density


def cycles(s1, s2):
    print()
    print("BROKEN CYCLES")
    for llmer in itertools.product(range(s), repeat=ll):
        llmer = "".join(str(x) for x in llmer)
        ismin = True
        llmer2 = llmer * 2
        p1 = set()
        p2 = set()
        for i in range(ll):
            llmer3 = llmer2[i : i + ll]
            if llmer3 < llmer:
                ismin = False
                break
        if not ismin:
            continue

        for i in range(ll):
            lmer = llmer2[i : i + l]
            p1.add((i + s1[lmer]) % ll)
            p2.add((i + s2[lmer]) % ll)
        if len(p1) == len(p2):
            continue

        print(llmer, p1, p2)
        for i in range(ll):
            lmer = llmer2[i : i + l]
            idx1 = s1[lmer]
            idx2 = s2[lmer]
            p1.add((i + idx1) % ll)
            p2.add((i + idx2) % ll)
            print(
                f'{" " * i}{lmer}{" " * (ll - i)}',
                f"{(i + idx2) % ll:>2}",
                f"{lmer[idx2:]:<14}",
                list((i + x) for x in lr2(lmer, all=True)[1]),
            )


# for lmer in scheme:
#     idx = scheme[lmer]
#     print(lmer, idx, lmer[idx:])

# k=1 s=2 w=4
alt = [("1110", 3), ("0100", 0), ("1000", 3), ("0000", 3), ("0101", 2)]
# k=1 s=2 w=5
alt = [
    ("00000", 0),
    ("00001", 3),
    ("00010", 2),
    ("00011", 2),
    ("00100", 1),
    ("00101", 1),
    ("01000", 0),
    ("01001", 0),
    ("01010", 0),
    ("01011", 2),
    ("01100", 0),
    ("01101", 0),
    ("01110", 0),
    ("10000", 1),
    ("10001", 3),
    ("10010", 2),
    ("10011", 2),
    ("10100", 1),
    ("10101", 1),
    ("11000", 2),
    ("11001", 3),
    ("11100", 3),
]
# k=1 s=2 w=6
# Avoid 01001
alt = {"010010": 3, "101001": 4}
# k=1 s=2 w=7
# Avoid 01001 and 010001
alt = {
    # "1101001": 5,
    "1010010": 4,
    "0100100": 3,
    # "0100010": 4,
    # "1010001": 5,
}
# k=1 s=2 w=8
# Avoid 01001 and 011001
# Prefer 010110 over 0110
alt8 = {
    "01100100": 4,
    "01100101": 4,
    "10110010": 5,
    "01011001": 6,
    "11011001": 6,
    "00101100": 1,
    "10010110": 2,
    "01001011": 3,
    "01001000": 3,
    "10100100": 4,
    "01010010": 5,
    "10101001": 6,
    "11010010": 5,
    "01000100": 4,
}

alt = {}

d1 = density(scheme)
d2 = density(lr)
print(f"{" "*ll}      {d1:<5.4} | {d2:<5.4}    {(d2-d1)*2**ll}")


# TRY DIFFERENT SUFFIXES
if True:
    for sl in range(0):
        for suf1 in itertools.product(range(s), repeat=sl):
            for suf2 in [[], [-0.5], [0.5], [1.5]]:
                SUF = list(suf1) + suf2
                d1 = density(scheme)
                d2 = density(lr)
                print(f"            {d1:<5.4} | {d2:<5.4}    {(d2-d1)*2**ll}   {SUF}")
    print()

lr = {x: lr(x) for x in scheme}

print("BEFORE UPDATES")
for lmer in scheme:
    idx = scheme[lmer]
    idx2 = lr[lmer]
    if idx != idx2:
        print(
            f"{lmer} {idx} {lmer[:idx]:>10}|{lmer[idx:]:<10} | {idx2} {lmer[:idx2]:>10}|{lmer[idx2:]:<10}"
        )
print()

# if d1 != d2:
#     print("diff: ", d2 - d1, (d2 - d1) * 2 ** (w + k))
for x, y in alt.items():
    assert scheme[x] == y
    # scheme[x] = y
    lr[x] = y
    d1 = density(scheme)
    d2 = density(lr)
    print(f"{x} {y}     {d1:<5.4} | {d2:<5.4}    {(d2-d1)*2**ll}")
print()


if alt:
    print("AFTER UPDATES")
    for lmer in scheme:
        idx = scheme[lmer]
        idx2 = lr[lmer]
        if idx != idx2:
            print(
                f"{lmer} {idx} {lmer[:idx]:>10}|{lmer[idx:]:<10} | {idx2} {lmer[:idx2]:>10}|{lmer[idx2:]:<10}"
            )
    print()
print(f"{' '*l}       {d1:<5.4} | {d2:<5.4}    {(d2-d1)*2**ll}")

cycles(scheme, lr)
exit(0)

g2 = [
    "1000",
    "1001",
    "1011",
    "1111",
    # "0001",
    "1101",
    "0101",
    "0011",
    "0010",
    "0111",
    # "0000",
    # asd
]

g3 = [
    "0020",
    "1020",
    "2020",
    "0120",
    "1120",
    "2120",
    "0220",
    "1220",
    "2220",
    #
    "0010",
    "1010",
    "2010",
    "0110",
    "1110",
    "0210",
    "1210",
    "2210",
    # "2110",
    "0211",
    "2110",  #
    #
    "0000",
    "1002",
    # "1221",
    # "0212",
    # "0211",
    # "0000",
    # "0111",
    # "0012",
]

g = ["011", "010", "110", "111", "000", "001"]

g = ["0001", "0101", "1101", "0100", "0111", "1100", "0011"]

g = ["110??", "010??", "0111?"]

g = ["?10", "001"]

g = ["1000", "1001", "1011"]

g = ["100011", "101111", "100001", "100110"]

g = ["?20", "02?", "01?"]
g = ["0001?", "?100?", "?101?", "01?1?", "?1?10", "00110"]
g = ["1000", "1001", "1011", "1010"]
g = ["0101", "0100", "1110", "0110", "0001", "0000", "0111", "1001", "1111"]
g = [
    "11001",
    "01001",
    "11011",
    "10001",
    "00001",
    "11110",
    "01011",
    "11111",
    "11100",
    "11000",
    "11010",
    "01010",
    "10000",
]
g = ["1110", "0110", "0100", "0001", "0000", "0101", "0111"]

g = ["1000", "1001", "1011", "1010"]
# g = []

d = None
for i in range(len(g) + 1):
    ss = greedy(scheme, g[:i])
    d2 = density(ss)
    if d is not None and d2 > d:
        print("              Increase at", g[i - 1])
    d = d2

all_kmers = set(lmer[idx : idx + k] for lmer, idx in scheme.items())
for x in g:
    all_kmers.add(x)
kmers = sorted(list(all_kmers))

# Reverse-engineer a greedy scheme
ans = []
for kmer in kmers:
    # Test if kmer is always preferred
    num_found = 0
    num_greedy = 0
    for lmer2 in scheme:
        if kmer not in lmer2:
            continue
        num_found += 1
        idx2 = ss(lmer2)
        kmer2 = lmer2[idx2 : idx2 + k]
        if kmer2 == kmer:
            num_greedy += 1
    if num_found == 0:
        continue
    ans.append((num_greedy / num_found, num_greedy, num_found, kmer))

for a, b, c, d in sorted(ans):
    gg = "g" if find(g, d) is not None else ""
    print(f"{d}: {b:>2}/{c:>2} {gg}")
