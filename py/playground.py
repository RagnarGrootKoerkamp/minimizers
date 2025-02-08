#!/usr/bin/env python3
# Trying out optimal-density schemes for k=1.

import itertools
import sys


def fmt(window):
    return "".join(str(x) for x in window)


def min_rot(window):
    return window == min(window[i:] + window[:i] for i in range(len(window)))


def density(f, w):
    ll = w + 1

    total_sampled = 0

    failing = 0

    # Iterate all length 2^ll 01 strings
    for s in itertools.product([0, 1], repeat=ll):
        s = list(s)
        ss = s + s
        sampled = set()
        for i in range(ll):
            window = ss[i : i + w]
            sampled.add((i + f(window)) % ll)
        cnt = len(sampled)
        if cnt > 2 and min_rot(s):
            if ll % cnt != 0 or s[: ll // cnt] * cnt != s:
                failing += 1
                print(f"\nFAILING CYCLE {failing:>3}")
                print(fmt(s), sampled)
                # print()
                for i in range(ll):
                    window = ss[i : i + w]
                    sampled.add((i + f(window, p=True)) % ll)

        total_sampled += len(sampled)
    print("failing: ", failing)
    d = total_sampled / (ll * (2**ll))
    print("density=", d, "\n")
    return d


def runs(vals, p=False):
    m = max(vals)
    # Find max values, and count how many non-max values follow.
    runs = []
    for i, x in enumerate(vals):
        if x == m:
            runs.append([1, i])
        elif len(runs) > 0:
            runs[-1][0] += 1
    return runs


def scheme1(window, p=False):
    p = False
    w = len(window)
    rs = runs(window)
    x = max(rs)[1]
    if p:
        print(fmt(window), rs, x)
    return x


# 2 levels
def scheme2(window, p=False):
    p = False
    w = len(window)
    rs = runs(window)
    if len(rs) == 1:
        x = rs[0][1]
        if p:
            print(fmt(window), rs, x)
    else:
        lens = [r[0] for r in rs]
        if p:
            print("lens", lens)
        rs2 = runs(lens, p)
        if p:
            print("runs", rs2)
        x2 = max(rs2, key=lambda x: (x[0], x[1]))[1]
        x = rs[x2][1]
        if p:
            print(fmt(window), rs, rs2, x2, x)
            print()
    return x


# 3 levels
def scheme3(window, p=False):
    # p = False
    w = len(window)
    rs = runs(window)
    if len(rs) == 1:
        x = rs[0][1]
        if p:
            print(fmt(window), rs, x)
    else:
        lens = [r[0] for r in rs]
        if lens[-2] < max(lens):
            lens[-1] = max(lens)
        if p:
            print("lens", lens)
        rs2 = runs(lens, p)
        if p:
            print("runs", rs2)
        if len(rs2) == 1:
            x2 = max(rs2, key=lambda x: (x[0], x[1]))[1]
        else:
            lens2 = [r[0] for r in rs2]
            if lens2[-2] < max(lens2):
                lens2[-1] = max(lens2)
            if p:
                print("lens", lens2)
            rs3 = runs(lens2, p)
            if p:
                print("runs", rs3)
            x3 = max(rs3, key=lambda x: (x[0], x[1]))[1]
            x2 = rs2[x3][1]

        x = rs[x2][1]
        if p:
            print(fmt(window), rs, rs2, x2, x)
            print()
    return x


def step_runs(vals, p=False):
    m = max(vals)
    # [ones+zeros, zeros, pos]
    runs = [[0, 0, -1]]
    for i, x in enumerate(vals):
        if x == m and runs[-1][1] == 0:
            # one, extend
            runs[-1][0] += 1
            runs[-1][2] = i
        elif x == m:
            # drop leading zeros
            if runs[-1][0] == runs[-1][1]:
                runs.pop()

            # one, new
            runs.append([1, 0, i])
        else:
            # zero, extend
            runs[-1][0] += 1
            runs[-1][1] += 1
            if runs[-1][0] == runs[-1][1]:
                runs[-1][2] = i
    if len(runs) > 1 and runs[-1][1] == 0:
        runs.pop()

    for i in range(len(runs)):
        runs[i][1] = runs[i][0] - runs[i][1]
        runs[i][0] -= runs[i][1]

    return runs


def step1(window, p=False):
    # p = False
    w = len(window)
    rs = step_runs(window)
    x = max(rs)[2]
    if p:
        print(fmt(window), x, rs)
    return x


def step2(window, p=False):
    # p = False
    w = len(window)
    rs = step_runs(window)
    x = max(rs)
    cnt = sum(int(x[0:2] == r[0:2]) for r in rs)
    if cnt == 1:
        x = max(rs)[2]
        if p:
            print(fmt(window), x, rs)
        return x
    lens = [r[0:2] for r in rs]
    rs2 = step_runs(lens)
    x2 = max(rs2, key=lambda r: r[0:2])[2]
    x = rs[x2][2]
    if p:
        print(fmt(window), x, x2, "cnt=", cnt, rs, "  ", rs2)
    return x


def step3(window, p=False):
    p = False
    w = len(window)

    # 1
    rs = step_runs(window)
    mx = max(rs)
    cnt = sum(int(mx[0:2] == r[0:2]) for r in rs)
    if cnt == 1:
        x = mx[2]
        if p:
            print(fmt(window), x, rs)
        return x

    # 2
    lens2 = [r[0:2] for r in rs]
    rs2 = step_runs(lens2)
    mx2 = max(rs2, key=lambda r: r[0:2])
    cnt = sum(int(mx2[0:2] == r[0:2]) for r in rs2)
    if cnt == 1:
        x2 = mx2[2]
        x = rs[x2][2]
        if p:
            print(fmt(window), x, x2, "cnt=", cnt, rs, "  ", rs2)
        return x

    # 3
    lens3 = [r[0:2] for r in rs2]
    rs3 = step_runs(lens3)
    mx3 = max(rs3, key=lambda r: r[0:2])

    x3 = mx3[2]
    x2 = rs2[x3][2]
    x = rs[x2][2]
    if p:
        print(fmt(window), x, x2, x3, "cnt=", cnt, rs, "  ", rs2, "  ", rs3)
    return x


def rindex(l, x):
    return len(l) - l[::-1].index(x) - 1


def findl(text, pat):
    for i in range(len(text) - len(pat) + 1):
        if text[i : i + len(pat)] == pat:
            return i
    return -1


def sus(window, p=False):
    w = len(window)
    best = (1, True, True, True, True, [9], -1)
    first1 = window.index(1) if 1 in window else 0
    last1 = rindex(window, 1) if 1 in window else 0

    if findl(window, [1, 0]) == -1:
        return first1

    for i in range(first1, len(window)):
        subs = list(window[i:])
        idx = first1 + findl(window[first1:], subs)
        uniq = idx < first1 or idx == i

        prefix_repeats = False
        for r in range(1, len(subs) + 1):
            if i < first1:
                continue
            rep = subs[:r]
            reps = rep * (w // r + 1)

            # if 0 in rep and 1 in rep and rindex(rep, 0) < rep.index(1):
            #     break

            if window[first1:i] == reps[len(reps) - i + first1 :]:
                prefix_repeats = True
            # if the prefix of len r is unique, we stop.
            if (
                # findl(window[0 : i + r - 1], rep)
                # == -1
                # and
                findl(window[i + 1 :], rep)
                == -1
            ):
                break

        # Everything excluding trailing 0s
        safe_prefix = window[i : last1 + 1 if i < last1 else None]
        safe_prefix_uniq = (
            i > last1 or first1 + findl(window[first1:], safe_prefix) == i
        )

        # Avoid 0..1..0$
        ends_in_0 = 1 in subs and findl(window, subs[:-1]) < i and subs[-1] == 0

        best = min(
            best,
            (
                subs[0],
                True or not safe_prefix_uniq,
                prefix_repeats,
                True or not uniq,
                True or ends_in_0,
                subs + [9],
                i,
            ),
        )

    if p:
        print(
            fmt(window),
            f"{best[-1]:>2}",
            f"{fmt(best[-2][:-1]):{w}}",
            "\t",
            best[0:3],
        )
    assert best[-1] != -1
    return best[-1]


def ladder(window, p=False):
    w = len(window)
    i = 0
    leading_ones = 0
    while i < w and window[i] == 1:
        leading_ones += 1
        i += 1
    runs = []
    poss = []
    for i in range(i, w):
        if window[i] == 0:
            runs.append(0)
            poss.append(i)
        else:
            runs[-1] += 1

    rs = len(runs)

    # 000000
    if rs == 0:
        return 0

    if rs == 1:
        x = poss[0]
        if p:
            print(fmt(window), f"{x:>2}", "UNIQ", rs, f"({leading_ones})", runs, poss)
        return x

    # (0) 1 2 3 4 4 4 3 2 1
    #                 ^ prefix
    #           ^ suffix

    # Find the smallest substring of runs
    # Do not choose from an increasing prefix or decreasing suffix
    prefix = 0
    last = leading_ones
    while prefix < rs:
        if last <= runs[prefix]:
            last = runs[prefix]
            prefix += 1
        else:
            break

    # Only increasing: return last 0
    if prefix == rs:
        x = poss[-1]
        if p:
            print(fmt(window), f"{x:>2}", "INC", rs, f"({leading_ones})", runs, poss)
        return x

    suffix = rs
    last = 0
    while suffix - 1 >= 0:
        if runs[suffix - 1] >= last:
            last = runs[suffix - 1]
            suffix -= 1
        else:
            break

    # Only decreasing: return first 0
    # if suffix == 0:
    #     x = poss[prefix]
    #     if p:
    #         print(fmt(window), f"{x:>2}", "DEC", rs, f"({leading_ones})", runs, poss)
    #     return x

    # assert suffix != prefix, "HUH??"

    # Increasing, then decreasing, without local minimum: return leftmost 0
    if suffix < prefix and False:
        x = poss[prefix]
        # x = poss[0]
        if p:
            print(
                fmt(window),
                f"{x:>2}",
                "TOP",
                rs,
                f"({leading_ones})",
                runs,
                poss,
            )
        return x

    # We have a local minimum.
    # (0) 1 2 3 4 2 4 3 2 1
    #             ^ prefix
    #               ^ suffix
    # assert suffix > prefix

    # Pick the smallest unique substring starting in [prefix..suffix)
    best = ([999], -1)
    for i in range(prefix, w):
        # Avoid unstable minima
        if findl(runs[prefix:], runs[i:-1]) < i - prefix:
            continue
        best = min(best, (runs[i:] + [999], i))
    unstable = False
    if best[1] == -1:
        unstable = True
        for i in range(prefix, w):
            best = min(best, (runs[i:] + [999], i))

    assert best[1] != -1
    x = poss[best[1]]
    if p:
        print(
            fmt(window),
            f"{x:>2}",
            "MIN",
            rs,
            f"({leading_ones})",
            runs,
            poss,
            "UNSTABLE" if unstable else "",
        )
    return x


w = 6
if len(sys.argv) > 1:
    w = int(sys.argv[1])

# density(scheme1, w)
# density(scheme2, w)
# density(scheme3, w)
# density(step1, w)
# density(step3, w)
# density(step2, w)

# density(step2, w)
# density(sus, w) # w=12 -> 45
density(ladder, w)  # w=13 -> 40
