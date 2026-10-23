"""Exact finite computations for fragile power domination. Standard library only."""
from collections import Counter
from math import comb, factorial, gcd


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def adjacency(n, edges):
    a = [0] * n
    for u, v in edges:
        if u == v or not (0 <= u < n and 0 <= v < n):
            raise ValueError("Expected a simple graph on range(n)")
        a[u] |= 1 << v
        a[v] |= 1 << u
    return a


def observed_mask(adj, active):
    """Initial closed neighbourhood followed by parallel forcing rounds."""
    seen = active
    for v in range(len(adj)):
        if active >> v & 1:
            seen |= adj[v]
    while True:
        extra = 0
        for v in range(len(adj)):
            if seen >> v & 1:
                missing = adj[v] & ~seen
                if missing and missing & (missing - 1) == 0:
                    extra |= missing
        if not extra:
            return seen
        seen |= extra


def coverage_table(n, edges, sensors):
    adj = adjacency(n, edges)
    return [observed_mask(adj, sum(1 << v for i, v in enumerate(sensors)
                                  if mask >> i & 1)).bit_count()
            for mask in range(1 << len(sensors))]


def mobius(values):
    n = len(values)
    if n == 0 or n & (n - 1):
        raise ValueError("Table length must be a power of two")
    result = list(values)
    for i in range(n.bit_length() - 1):
        for mask in range(n):
            if mask >> i & 1:
                result[mask] -= result[mask ^ (1 << i)]
    return result


def layer_sums(values):
    a = [0] * len(values).bit_length()
    for mask, value in enumerate(values):
        a[mask.bit_count()] += value
    return a


def diagonal_coefficients(values):
    a = layer_sums(values)
    s = len(a) - 1
    return [sum((-1)**(j-k) * choose(s-k, j-k) * a[k]
                for k in range(j+1)) for j in range(s+1)]


def alpha(s, ell, k, i):
    if not 1 <= i <= ell <= s or not 1 <= k <= s:
        raise ValueError("Expected 1 <= i <= ell <= s and 1 <= k <= s")
    if k <= ell:
        return int(k == i)
    return (-1)**(ell-i) * choose(s-i, k-i) * choose(k-i-1, ell-i)


def gadget(kind, d):
    """P: d sensor leaves plus an unsensored path of length two. N: star."""
    if d < 2 or kind not in ("P", "N"):
        raise ValueError("Expected kind P/N and d >= 2")
    edges = [(0, i) for i in range(1, d+1)]
    if kind == "P":
        edges += [(0, d+1), (d+1, d+2)]
    return d + (3 if kind == "P" else 1), edges, list(range(1, d+1))


def gadget_polynomial(kind, d):
    p = [0] * (d+1)
    p[1] = 2*d
    for j in range(2, d+1):
        p[j] = (-1)**(j+1) * comb(d, j)
    if kind == "P":
        p[d] += 2
    else:
        p[d-1] += d
        p[d] -= d
    return p


def forest_recipe(d, reduce_multipliers=True):
    """Compressed integer construction, with optional gcd cancellation.

    Set reduce_multipliers=False to reproduce the original construction.
    Neither version is claimed to minimize the number of vertices.
    """
    counts = Counter({("P", d): 1})
    polynomial = gadget_polynomial("P", d)
    for j in range(d, 1, -1):
        c = polynomial[j]
        if c == 0:
            continue
        kind = "N" if c > 0 else "P"
        g = gadget_polynomial(kind, j)
        scale, copies = abs(g[j]), abs(c)
        if reduce_multipliers:
            common = gcd(scale, copies)
            scale //= common
            copies //= common
        counts = Counter({key: value * scale for key, value in counts.items()})
        counts[(kind, j)] += copies
        polynomial = [x * scale for x in polynomial]
        for k, x in enumerate(g):
            polynomial[k] += copies * x
    return counts, polynomial


def forest_size_bound(d):
    """The proved upper bound (5/2) d! (d+1)!, valid for both recipes."""
    if d < 2:
        raise ValueError("Expected d >= 2")
    return 5 * factorial(d) * factorial(d+1) // 2


def recipe_statistics(counts):
    return {
        "vertices": sum(c*(d+(3 if kind == "P" else 1))
                        for (kind, d), c in counts.items()),
        "sensors": sum(c*d for (kind, d), c in counts.items()),
        "components": sum(counts.values()),
    }


def materialize_recipe(counts, max_vertices=10000):
    stats = recipe_statistics(counts)
    if stats["vertices"] > max_vertices:
        raise ValueError("Recipe exceeds materialization limit; use compressed counts")
    n, edges, sensors = 0, [], []
    for (kind, d), copies in sorted(counts.items(), reverse=True):
        nv, es, ss = gadget(kind, d)
        for _ in range(copies):
            edges.extend((u+n, v+n) for u, v in es)
            sensors.extend(v+n for v in ss)
            n += nv
    return n, edges, sensors
