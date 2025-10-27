#!/usr/bin/env python3

"""
Script serves as a "load balancer" for samestr convert and samestr merge outputs.
Output sizes are quite heterogeneous, so some kind of balancing is required 
to avoid swarms of tiny < 1min jobs as well as allowing sample-rich clades
to be computed in time.

samestr convert outputs 1 file (matrix) per clade per sample -> optimisation by file count
samestr merge outputs 2 files (sample-names, matrix) per clade -> optimisation by file size

Method is a naive bucket-fill approach. For each file (pair) item check if it fits into one 
of a series of buckets. If not, add a new bucket containing the item.
"""

import argparse
import os




def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir")
	ap.add_argument("max_size", type=int)
	ap.add_argument("criterion", type=str, choices=("count", "size",), default="count")
	args = ap.parse_args()

	_, _, files = next(os.walk(args.input_dir))

	"t__SGB15164.SAMN17486269.mp4.npz"

	d = {}
	for f in files:
		if f.endswith("_clades.txt"):
			with open(f, "rt") as _in:
				for path in _in.read().rstrip().split("\n"):
					# /scratch/schudoma/WORK/mpow_ssdev/flow_test/e4/5394ca9f37fb1d77be77d211c6e5e2/sstr_convert/SAMN17485756.mp4/t__SGB1861.SAMN17485756.mp4.npz
					# /path/to/t__SGB14039.names.txt
					ff = path[path.rfind("/") + 1:]
					d.setdefault(ff[:ff.find(".")], []).append(path)

	batches = [[]]
	batch_sizes = [0]

	if args.criterion == "count":
		fbatch_sizes = {k: len(v) for k, v in d.items()}

	else:
		fbatch_sizes = {
			k: len(open(v[0] if v[0].endswith(".txt") else v[1]).readlines())
			for k, v in d.items()
		}

	for fbid, fbsize in fbatch_sizes.items():
		for i, (batch, bsize) in enumerate(zip(batches, batch_sizes)):
			if bsize + fbsize <= args.max_size:
				batch += d[fbid]
				batch_sizes[i] += fbsize
				break
		else:
			batches.append(d[fbid])
			batch_sizes.append(fbsize)
		
	
	for i, (batch, bsize) in enumerate(zip(batches, batch_sizes)):
		for f in batch:
			print(i, bsize, f, sep="\t")


if __name__ == "__main__":
	main()
