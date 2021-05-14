#!/usr/bin/env python3

import argparse
import json
import os
import random
import sys

parser = argparse.ArgumentParser(description="Instance generator")
parser.add_argument(
    "-n", "--n-items", dest="n_items", type=int, required=True, help="number of items"
)
parser.add_argument(
    "-m",
    "--n-stacks",
    dest="n_stacks",
    type=int,
    required=True,
    help="number of stacks",
)
parser.add_argument(
    "-b", "--capacity", dest="capacity", type=int, help="capacity of stacks"
)
parser.add_argument("-o", "--output", dest="output", type=str, help="save to a file")
parser.add_argument(
    "--blocking-matrix",
    dest="blocking_density",
    type=float,
    help="generate a blocking matrix (with given density) instead of departures",
)
parser.add_argument(
    "--conflict-matrix",
    dest="conflict_density",
    type=float,
    help="generate a conflict matrix (with given density) instead of weights",
)
parser.add_argument("--random-seed", dest="random_seed", type=int, help="random seed")
args = parser.parse_args()

if args.random_seed is not None:
    random_seed = args.random_seed
else:
    random_seed = int.from_bytes(os.urandom(8), sys.byteorder)

random.seed(random_seed)

data = {
    "n_items": args.n_items,
    "n_stacks": args.n_stacks,
    "capacity": args.capacity if args.capacity else args.n_items,
    "random_seed": random_seed,
}

if args.blocking_density is not None:
    data["blocking_matrix"] = [
        [
            1 if random.random() < args.blocking_density else 0
            for _ in range(args.n_items)
        ]
        for _ in range(args.n_items)
    ]
else:
    departures = list(range(1, args.n_items + 1))
    random.shuffle(departures)
    data["departures"] = departures

if args.conflict_density is not None:
    data["conflict_matrix"] = [
        [
            1 if random.random() < args.conflict_density else 0
            for _ in range(args.n_items)
        ]
        for _ in range(args.n_items)
    ]
else:
    weights = list(range(1, args.n_items + 1))
    random.shuffle(weights)
    data["weights"] = weights

if args.output:
    with open(args.output, mode="w") as fp:
        json.dump(data, fp, separators=(",", ":"))
else:
    print(json.dumps(data, separators=(",", ":")))
