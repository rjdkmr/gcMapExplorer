#!/usr/bin/env python3
#
# This file is part of gcMapExplorer
# Copyright (C) 2016-2018  Rajendra Kumar, Ludvig Lizana, Per Stenberg
#
# gcMapExplorer is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gcMapExplorer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gcMapExplorer.  If not, see <http://www.gnu.org/licenses/>.
#
# =============================================================================
from argparse import ArgumentParser, ArgumentTypeError

from gcMapExplorer.lib.hic2gcmap import hic2gcmap
from gcMapExplorer.lib.hicparser import HicFileType, NormType


def resolution(string):
    if string == "finest":
        return string

    if string.endswith("kb"):
        string = string[:-2] + "000"

    try:
        value = int(string)
    except ValueError:
        raise ArgumentTypeError("Resolution must be a positive integer or integer ending with 'kb'")
    else:
        if value <= 0:
            raise ArgumentTypeError("Resolution must be a positive integer")
        return value


def norm(string):
    if string == "none":
        return None
    return NormType(string)


def main():
    from os import path
    import sys

    parser = ArgumentParser(prog="gcMapExplorer hic2gcmap", description="Convert hic files to gcmap",
                            allow_abbrev=False)
    parser.add_argument("input", type=HicFileType(), help="hic input file")
    parser.add_argument("output", type=str, nargs="?", help="output file or directory", default=".")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-c", "--chromosomes", type=str, nargs=2, metavar=("A", "B"), help="a pair of chromosomes A B")
    group.add_argument("-l", "--list", action="store_true", help="list all available chromosomes")
    parser.add_argument("--compression", type=str, choices=("lzf", "gzip"), default="lzf", metavar="C",
                        help="compression type, choose between lzf, gzip (default: lzf)")
    parser.add_argument("-r", "--resolution", type=resolution, default="finest", metavar="R",
                        help="the resolution, as an integer or as kb (default: finest)")
    parser.add_argument("-n", "--norm", type=str, choices=("VC", "VC_SQRT", "KR", "none"), metavar="N",
                        default="none",
                        help="the type of norm to use, choose between VC, VC_SQRT, KR, none (default: none)")
    parser.add_argument("--downsampling", type=str, choices=("sum", "mean", "max", "none"), default="sum", metavar="D",
                        help="the downsampling method to use, choose between sum, mean, max, none (default: sum)")

    args = parser.parse_args(args=sys.argv[sys.argv.index("hic2gcmap") + 1:])

    hic = args.input

    if args.list:
        print(", ".join(hic.chromosomes))
        return

    if path.isdir(args.output):
        filename = path.splitext(path.basename(args.input.name))[0]
        output_file = path.join(args.output, filename)
        if args.chromosomes:
            output_file += "_" + "_".join(args.chromosomes)
        if args.downsampling == "none" and args.resolution != "finest":
            output_file += "_" + str(args.resolution)[:-3] + "kb"
        if args.norm != "none":
            output_file += "_" + args.norm
        output_file += "_" + args.compression
        output_file += ".gcmap"
    else:
        output_file = args.output

    if args.chromosomes:
        chr1, chr2 = args.chromosomes
        try:
            hic2gcmap(hic, chr1, chr2, output_file, resolution=args.resolution, norm_type=norm(args.norm),
                      compression=args.compression,
                      downsampling=args.downsampling)
        except Exception as e:
            print("Error: {}".format(e))
    else:  # all chromosomes
        chromosome_names = filter(lambda pair: "All" not in pair, (hic.chromosome_names(r) for r in hic.records))
        for chr1, chr2 in chromosome_names:
            try:
                hic2gcmap(hic, chr1, chr2, output_file, resolution=args.resolution, norm_type=norm(args.norm),
                          compression=args.compression,
                          downsampling=args.downsampling)
            except Exception as e:
                print("Error: {}".format(e))


if __name__ == "__main__":
    main()
