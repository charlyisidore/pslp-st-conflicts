#!/usr/bin/env python3

import argparse
import json
import os
from functools import cmp_to_key


# Describe conflicts using a list of edges instead of a binary matrix
use_conflict_graph = False

# Fix conflicts between initial items of the same stack
fix_initial_conflicts = False


parser = argparse.ArgumentParser(description="Instance converter")
parser.add_argument(
    "-i", "--input", dest="input", type=str, required=True, help="original instance"
)
parser.add_argument("-o", "--output", dest="output", type=str, help="save to a file")
args = parser.parse_args()


# Problem
def read_slp(params, data):
    n_items = int(params[0])
    n_stacks = int(params[1])
    capacity = int(params[2])
    # print(f"slp n_items={n_items} n_stacks={n_stacks} capacity={capacity}")
    data["n_items"] = n_items
    data["n_stacks"] = n_stacks
    data["capacity"] = capacity
    data["arrivals"] = [0] * n_items
    data["departures"] = [0] * n_items
    data["initial_positions"] = [0] * n_items
    data["initial_levels"] = [0] * n_items
    data["sc"] = []


# Item
def read_item(params, data):
    i = int(params[0]) - 1
    arrival = int(params[1])
    departure = int(params[3]) + 1
    # print(f"item i={i} arrival={arrival} departure={departure}")
    assert i >= 0 and i < data["n_items"]
    assert arrival == int(params[2])
    data["arrivals"][i] = arrival
    data["departures"][i] = departure


# Initial items
def read_fix(params, data):
    i = int(params[0]) - 1
    k = int(params[1])
    l = int(params[2])
    # print(f"initial i={i} stack={k} level={l}")
    assert i >= 0 and i < data["n_items"]
    assert k > 0 and k <= data["n_stacks"]
    assert l > 0 and l <= data["capacity"]
    assert data["arrivals"][i] == 0
    data["initial_positions"][i] = k
    data["initial_levels"][i] = l


# Conflicts
def read_sc(params, data):
    i = int(params[0]) - 1
    j = int(params[1]) - 1
    # print(f"conflict i={i} j={j}")
    assert i >= 0 and i < data["n_items"]
    assert j >= 0 and j < data["n_items"]
    data["sc"].append((i, j))


data = {}

with open(args.input, mode="r") as fp:
    for line in fp.readlines():
        params = [x.strip() for x in line.split()]
        if len(params) == 0:
            continue
        switcher = {"slp": read_slp, "item": read_item, "fix": read_fix, "sc": read_sc}
        reader = switcher.get(params[0], lambda: "invalid line type")
        reader(params[1:], data)


def reorder(data):
    def compare(i, j):
        a_i = data["arrivals"][i]
        a_j = data["arrivals"][j]
        if a_i < a_j:
            return -1
        elif a_i > a_j:
            return 1
        elif a_i == 0:
            k_i = data["initial_positions"][i]
            k_j = data["initial_positions"][j]
            if k_i < k_j:
                return -1
            elif k_i > k_j:
                return 1
            else:
                l_i = data["initial_levels"][i]
                l_j = data["initial_levels"][j]
                if l_i < l_j:
                    return -1
                elif l_i > l_j:
                    return 1

        return -1 if i < j else 1 if i > j else 0

    order = sorted(range(data["n_items"]), key=cmp_to_key(compare))

    reverse = [0] * len(order)
    for i in range(len(order)):
        reverse[order[i]] = i

    n_items = data["n_items"]
    n_initial_items = len([1 for t in data["arrivals"] if t == 0])

    result = {
        "n_items": n_items,
        "n_stacks": data["n_stacks"],
        "capacity": data["capacity"],
        "arrivals": [0] * n_items,
        "departures": [0] * n_items,
        "initial_positions": [0] * n_initial_items,
        "initial_levels": [0] * n_initial_items,
        "indices": [i + 1 for i in order],
    }

    for i in range(n_items):
        k = order[i]
        result["arrivals"][i] = data["arrivals"][k]
        result["departures"][i] = data["departures"][k]
        if i < n_initial_items:
            result["initial_positions"][i] = data["initial_positions"][k]
            result["initial_levels"][i] = data["initial_levels"][k]

    # data["sc"] contains allowed pairs
    # conflict graph/matrix contains forbidden pairs
    matrix = [[1] * n_items for _ in range(n_items)]
    for (i, j) in data["sc"]:
        matrix[reverse[i]][reverse[j]] = 0

    # Fix forbidden stackings in initial items
    if fix_initial_conflicts:
        for i in range(n_initial_items):
            for j in range(i):
                assert result["arrivals"][i] == 0
                assert result["arrivals"][j] == 0
                assert result["initial_positions"][i] >= result["initial_positions"][j]
                if result["initial_positions"][i] == result["initial_positions"][j]:
                    assert result["initial_levels"][i] > result["initial_levels"][j]
                    matrix[i][j] = 0

    conflicts = [
        (i, j)
        for i in range(n_items)
        for j in range(n_items)
        if i != j and matrix[i][j] == 1
    ]

    if use_conflict_graph:
        result["conflict_graph"] = sorted((i + 1, j + 1) for (i, j) in conflicts)
    else:
        result["conflict_matrix"] = [[0] * n_items for _ in range(n_items)]
        for (i, j) in conflicts:
            result["conflict_matrix"][i][j] = 1

    return result


for i in range(data["n_items"]):
    assert (data["arrivals"] == 0) == (
        data["initial_positions"] == 0 and data["initial_levels"] == 0
    )

data = reorder(data)

if args.output:
    with open(args.output, mode="w") as fp:
        json.dump(data, fp, separators=(",", ":"))
else:
    print(json.dumps(data, separators=(",", ":")))
