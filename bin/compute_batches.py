#!/usr/bin/env python3

import argparse
import os
import pathlib


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("input_dir")
	ap.add_argument("max_size", type=int)
	args = ap.parse_args()

	input_dir = pathlib.Path(args.input_dir)
	pwd, dirs, files = next(os.walk(args.input_dir))

	"t__SGB15164.SAMN17486269.mp4.npz"

	d = {}
	for f in files:
		with open(f, "rt") as _in:
			for path in _in.read().rstrip().split("\n"):
				# /scratch/schudoma/WORK/mpow_ssdev/flow_test/e4/5394ca9f37fb1d77be77d211c6e5e2/sstr_convert/SAMN17485756.mp4/t__SGB1861.SAMN17485756.mp4.npz
				ff = path[path.rfind("/") + 1:]
				# clade = ff[:ff.find(".")]
				d.setdefault(ff[:ff.find(".")], []).append(path)

	fbatch_sizes = {k: len(v) for k, v in d.items()}
	
	batches = [[]]
	for fbid, fbsize in fbatch_sizes.items():
		for batch in batches:
			if len(batch) + fbsize <= args.max_size:
				batch += d[fbid]
				break
		else:
			batches.append(d[fbid])
	
	for i, batch in enumerate(batches):
		for f in batch:
			# print(i, (input_dir / f).resolve(), sep="\t")
			print(i, f, sep="\t")





if __name__ == "__main__":
	main()