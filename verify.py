"""Run exact checks and write verification.json. No third-party dependencies."""
import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from fragile import (choose, alpha, coverage_table, mobius, layer_sums,
                     diagonal_coefficients, gadget, gadget_polynomial,
                     forest_recipe, recipe_statistics, materialize_recipe)


def observed_reference(n, edges, sensors):
    """Independent implementation: adjacency sets and one force per iteration."""
    neighbours = [set() for _ in range(n)]
    for u, v in edges:
        neighbours[u].add(v)
        neighbours[v].add(u)
    seen = set(sensors)
    for v in sensors:
        seen.update(neighbours[v])
    while True:
        forced = None
        for v in sorted(seen):
            remaining = neighbours[v] - seen
            if len(remaining) == 1:
                forced = next(iter(remaining))
                break
        if forced is None:
            return len(seen)
        seen.add(forced)


def evaluate_direct(values, probabilities):
    result = Fraction(0)
    for mask, value in enumerate(values):
        weight = Fraction(1)
        for i, p in enumerate(probabilities):
            weight *= p if mask >> i & 1 else 1-p
        result += value * weight
    return result


def evaluate_mobius(mu, probabilities):
    result = Fraction(0)
    for mask, value in enumerate(mu):
        weight = Fraction(1)
        for i, p in enumerate(probabilities):
            if mask >> i & 1:
                weight *= p
        result += value * weight
    return result


def run():
    report = {"arithmetic": "Python exact integers and fractions.Fraction"}
    count = 0
    for s in range(1,31):
        for ell in range(1,s+1):
            for k in range(1,s+1):
                recurrent = {}
                for i in range(ell,0,-1):
                    recurrent[i] = choose(s-i,k-i) - sum(
                        recurrent[j]*choose(s-i,j-i) for j in range(i+1,ell+1))
                    assert recurrent[i] == alpha(s,ell,k,i)
                    count += 1
    report["alpha_checks"] = count
    report["alpha_max_s"] = 30

    graphs = placements = closures = hetero = 0
    for n in range(1,6):
        possible = list(combinations(range(n),2))
        for edge_mask in range(1 << len(possible)):
            edges = [e for j,e in enumerate(possible) if edge_mask >> j & 1]
            whole = coverage_table(n, edges, list(range(n)))
            graphs += 1
            for mask,value in enumerate(whole):
                assert value == observed_reference(n,edges,
                    [v for v in range(n) if mask >> v & 1])
                closures += 1
            for sensor_mask in range(1 << n):
                ss = [v for v in range(n) if sensor_mask >> v & 1]
                values = [whole[sum(1 << v for i,v in enumerate(ss) if sub >> i & 1)]
                          for sub in range(1 << len(ss))]
                mu = mobius(values)
                b = diagonal_coefficients(values)
                assert b == layer_sums(mu)
                a = layer_sums(values)
                s = len(ss)
                for ell in range(1,s+1):
                    lhs = all(x == 0 for x in b[ell+1:])
                    rhs = all(a[k] == sum(alpha(s,ell,k,i)*a[i]
                               for i in range(1,ell+1)) for k in range(1,s+1))
                    assert lhs == rhs
                ps = [Fraction(v+1,n+2) for v in ss]
                assert evaluate_direct(values,ps) == evaluate_mobius(mu,ps)
                hetero += 1
                placements += 1
    report.update(graphs=graphs,placements=placements,
                  independent_closure_checks=closures,heterogeneous_checks=hetero,
                  graph_scope="All labelled simple graphs on 1 through 5 vertices; all sensor subsets, including empty")

    gadget_cases = 0
    for d in range(2,9):
        for kind in ("P","N"):
            n,edges,sensors = gadget(kind,d)
            values = coverage_table(n,edges,sensors)
            for mask,value in enumerate(values):
                k = mask.bit_count()
                expected = (k + int(k>0) + 2*int(k==d)) if kind == "P" else (
                    0 if k==0 else d+1 if k>=d-1 else k+1)
                assert value == expected
                assert value == observed_reference(n,edges,
                    [v for i,v in enumerate(sensors) if mask >> i & 1])
                gadget_cases += 1
            assert diagonal_coefficients(values) == gadget_polynomial(kind,d)
    report["gadget_subset_checks"] = gadget_cases
    report["gadget_max_d"] = 8

    counts,polynomial = forest_recipe(2)
    n,edges,sensors = materialize_recipe(counts)
    vals = coverage_table(n,edges,sensors)
    assert n == 18 and len(sensors) == 8
    assert diagonal_coefficients(vals) == [0,18]+[0]*7
    mu = mobius(vals)
    assert max(m.bit_count() for m,x in enumerate(mu) if x) == 2
    assert sorted(x for m,x in enumerate(mu) if m.bit_count()==2 and x) == [-3,1,1,1]
    for mask,value in enumerate(vals):
        assert value == observed_reference(n,edges,
            [v for i,v in enumerate(sensors) if mask >> i & 1])
    report["example_survival_patterns"] = len(vals)
    report["example_vertices"] = n
    report["example_sensors"] = len(sensors)
    report["example_quadratic_interactions"] = [-3,1,1,1]

    recipes=[]
    for d in range(2,9):
        counts,p = forest_recipe(d)
        stats=recipe_statistics(counts)
        assert p == [0,stats["vertices"]]+[0]*(d-1)
        rebuilt=[0]*(d+1)
        for (kind,k),copies in counts.items():
            for j,x in enumerate(gadget_polynomial(kind,k)):
                rebuilt[j] += copies*x
        assert rebuilt == p
        assert counts[("P",d)] > 0
        recipes.append(dict(degree=d,**stats,coefficients=p,
            blocks=[dict(kind=kind,sensors=k,copies=c) for (kind,k),c in sorted(counts.items())]))
    report["compressed_forest_recipes"] = recipes
    report["status"] = "PASS"
    return report


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="verification.json")
    args=parser.parse_args()
    result=run()
    Path(args.output).write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "compressed_forest_recipes"},indent=2))
