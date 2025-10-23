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
		d.setdefault(f[:f.find(".")], []).append(f)

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
			print(i, (input_dir / f).resolve(), sep="\t")




if __name__ == "__main__":
	main()