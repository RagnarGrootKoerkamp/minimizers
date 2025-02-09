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


def sus2(window, p=False):
    w = len(window)
    best = (1, True, True, True, True, [9], -1)
    first1 = window.index(1) if 1 in window else 0
    last1 = rindex(window, 1) if 1 in window else 0

    if findl(window, [1, 0]) == -1:
        return first1

    for i in range(first1, len(window)):
        subs = list(window[i:])
        # for j in range(len(subs) - 1, 0, -1):
        #     subs[j] ^= subs[j - 1]
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

    mode = "MIN"
    # Only decreasing: return first 0
    if suffix == 0:
        mode = "DEC"
    if suffix == 0 and False:
        x = poss[prefix - 1]
        if p:
            print(fmt(window), f"{x:>2}", mode, rs, f"({leading_ones})", runs, poss)
        return x

    # assert suffix != prefix, "HUH??"

    # Increasing, then decreasing, without local minimum: return leftmost 0
    if suffix < prefix:
        mode = "TOP"
    if suffix < prefix and False:
        x = poss[prefix - 1]
        # x = poss[0]
        if p:
            print(
                fmt(window),
                f"{x:>2}",
                mode,
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
            mode,
            rs,
            f"({leading_ones})",
            runs,
            poss,
            "UNSTABLE" if unstable else "",
        )
    # TODO:
    # Given some suffix S=00..1..,
    # say the length p prefix makes it unique.
    # Then no suffix (extended with zeros) of it may be less or equal than S[:p].

    return x


def lcp(a, b):
    i = 0
    while i < len(a) and i < len(b) and a[i] == b[i]:
        i += 1
    return i


# Position of smallest unique substring
def sus(window, p=False):
    w = len(window)
    best = ([999], -1)
    uniq_len = 0
    for i in range(w):
        s = window[i:] + [999]
        u = lcp(s, best[0]) + 1
        if s < best[0]:
            best = (s, i)
            uniq_len = u
        else:
            uniq_len = max(uniq_len, u)
    return (best[1], uniq_len)


# Find SUS of window+0+window+0.
# Then return that if it's inside window.
# Otherwise find SUS of the complement.
def rot_sus(window, p=False):
    if 1 not in window:
        return 0

    # Find SUS in window+0+window
    w = len(window)
    ww = window + [0] + window + [0]

    best = ([999], -1)
    uniq_len = 0
    for i in range(w + 1):
        s = ww[i : i + w + 1]
        l = lcp(s, best[0])
        if s < best[0]:
            best = (s, i)
            uniq_len = l + 1
        else:
            uniq_len = max(uniq_len, l + 1)
    start = best[1]
    end = start + uniq_len
    if end <= w:
        if p:
            print(
                fmt(window),
                f"{start:>2}",
                f"{uniq_len:>2}",
                "A",
                fmt(best[0][:uniq_len]),
            )
        return start

    # Plain sus is the same? Return it also.
    (s2, u2) = sus(window, p)
    if start == s2:
        if p:
            print(
                fmt(window),
                f"{start:>2}",
                f"{uniq_len:>2}",
                "B",
                fmt(best[0][:uniq_len]),
            )
        return start

    if p:
        print(
            fmt(window),
            f"{start:>2}",
            f"{uniq_len:>2}",
            "_",
            fmt(best[0][:uniq_len]),
            "i range",
            end,
            start + w + 1,
        )

    # SUS overlaps the appended 0.
    # Find the next-best SUS in the complement.
    best = ([999], -1)
    for i in range(0, start):
        s = ww[i : i + w + 1]
        l = lcp(s, best[0])
        if s < best[0]:
            best = (s, i)
            uniq_len = l + 1
        else:
            uniq_len = max(uniq_len, l + 1)
    assert best[1] != -1
    start = best[1] % (w + 1)
    end = (start + uniq_len - 1) % (w + 1) + 1
    if p:
        print(
            fmt(window), f"{start:>2}", f"{uniq_len:>2}", "C", fmt(best[0][:uniq_len])
        )
    return start


# Find smallest substring of window+1+window+1
# that starts inside the window, and skips the leading zeros.
def new(window, p=False):
    if 1 not in window:
        return 0
    suf = 1
    w = len(window)
    ww = window + [suf] + window + [suf]
    best = ([999], -1)
    uniq_len = 0
    first1 = window.index(1)
    for i in range(first1, w):
        s = ww[i : i + w + 1]
        u = lcp(s, best[0]) + 1
        if s <= best[0]:
            best = (s, i)
            uniq_len = u
        else:
            uniq_len = max(uniq_len, u)

    for r in range(uniq_len, 1, -1):
        step = False
        if best[1] + r < w and best[0][0:r] == best[0][r : 2 * r]:
            step = True
            if p:
                print(
                    fmt(window),
                    f"{best[1]:>2}",
                    f"{uniq_len:>2}",
                    fmt(best[0][:uniq_len]),
                    "before",
                )
            best = (best[0][r:], best[1] + r)
            uniq_len -= r
        if step:
            break

    if p:
        print(fmt(window), f"{best[1]:>2}", f"{uniq_len:>2}", fmt(best[0][:uniq_len]))
    return best[1]


# Find the smallest suffix > 000.
# Then pick the position of the first 1.
def new2(window, p=False):
    if 1 not in window:
        return 0
    w = len(window)
    best = ([999], -1)
    uniq_len = 0
    for i in range(w):
        s = window[i:] + [999]
        if 1 not in s:
            continue

        u = lcp(s, best[0]) + 1
        if s <= best[0]:
            best = (s, i)
            uniq_len = u
        else:
            uniq_len = max(uniq_len, u)

    x = best[1] + best[0].index(1)

    if best[0][:4] == [0, 1, 0, 1] or best[0] == [0, 1, 0, 999]:
        x += 2

    if p:
        print(
            fmt(window),
            f"{x:>2}",
            f"{best[1]:>2}",
            f"{uniq_len:>2}",
            fmt(best[0][:uniq_len]),
        )
    return x


# Optimal for k=1,s=2,w<=9
def lr2(lmer, p=False):
    l = len(lmer)
    # Find position of last 0 in the longest 000..0011..1111 run.
    best = (0, 0, -100, -100)
    bests = [l - 1]
    for i in range(l - 1):
        if lmer[i] == 0 and lmer[i + 1] == 1:
            c0 = 0
            d0 = 0
            for j in range(i, -1, -1):
                if lmer[j] == 0:
                    c0 += 1
                else:
                    break
            for j in range(j, -1, -1):
                if lmer[j] == 1:
                    d0 += 1
                else:
                    break
            c1 = 0
            d1 = 0
            for j in range(i + 1, l):
                if lmer[j] == 1:
                    c1 += 1
                else:
                    break
            for j in range(j, l):
                if lmer[j] == 0:
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

    if p:
        print(fmt(lmer), bests[0], fmt(lmer[bests[0] :]), bests)

    if len(bests) == 1 or False:
        return bests[0]

    x = bests[-1]
    for i in range(len(bests) - 1):
        if (bests[i + 1] - bests[i]) % 2 == 0:
            x = bests[i]
            break

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


# Smallest one, but prevent taking a substring that can be left-extended into something smaller.
def sus3(window, p=False):
    w = len(window)
    best = (1, True, [9], -1)
    first1 = window.index(1) if 1 in window else 0

    if findl(window, [1, 0]) == -1:
        return first1

    for i in range(len(window)):
        s = window[i:]
        idx = first1 + findl(window[first1:], s)
        uniq = idx < first1 or idx == i

        prefix_repeats = False or i == w - 1
        for r in range(1, len(s) + 1):
            rep = s[:r]
            reps = rep * (i // r + 1)
            # Check if window[:i] is just copies of rep.
            # If `rep` extends to the left, also check that it matches the end of the seq.
            if window[:i] != reps[len(reps) - i :]:
                continue
            tail = len(reps) - i - 1
            if tail > 0 and window[-tail:] != reps[:tail]:
                continue
            prefix_repeats = True

        best = min(
            best,
            (
                s[0],
                prefix_repeats,
                s + [9],
                i,
            ),
        )

    x = best[-1]
    s = best[-2][:-1]

    if p:
        print(
            fmt(window),
            f"{best[-1]:>2}",
            f"{fmt(best[-2][:-1]):{w}}",
        )
    assert best[-1] != -1
    return best[-1]


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
# density(sus2, w)  # w=12 -> 45
# density(ladder, w)  # w=12 -> 40
# density(rot_sus, w)
# density(new, w)
# density(new2, w)  # w=12 -> 100
# density(lr2, w)  # w=12 -> 100

density(sus3, w)  # w=12 -> 45
