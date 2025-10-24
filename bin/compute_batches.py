#!/usr/bin/env python3

import argparse
import os
import pathlib


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir")
	ap.add_argument("max_size", type=int)
	ap.add_argument("criterion", type=str, choices=("count", "size",), default="count")
	args = ap.parse_args()

	input_dir = pathlib.Path(args.input_dir)
	pwd, dirs, files = next(os.walk(args.input_dir))

	"t__SGB15164.SAMN17486269.mp4.npz"

	d = {}
	for f in files:
		if f.endswith("_clades.txt"):
			with open(f, "rt") as _in:
				for path in _in.read().rstrip().split("\n"):
					# /scratch/schudoma/WORK/mpow_ssdev/flow_test/e4/5394ca9f37fb1d77be77d211c6e5e2/sstr_convert/SAMN17485756.mp4/t__SGB1861.SAMN17485756.mp4.npz
					# /path/to/t__SGB14039.names.txt
					ff = path[path.rfind("/") + 1:]
					# clade = ff[:ff.find(".")]
					d.setdefault(ff[:ff.find(".")], []).append(path)

	batches = [[]]
	if args.criterion == "count":
		fbatch_sizes = {k: len(v) for k, v in d.items()}
		
		for fbid, fbsize in fbatch_sizes.items():
			for batch in batches:
				if len(batch) + fbsize <= args.max_size:
					batch += d[fbid]
					break
			else:
				batches.append(d[fbid])
	else:
		batch_sizes = [0]
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
		
	
	for i, batch in enumerate(batches):
		for f in batch:
			# print(i, (input_dir / f).resolve(), sep="\t")
			print(i, f, sep="\t")
		





if __name__ == "__main__":
	main()