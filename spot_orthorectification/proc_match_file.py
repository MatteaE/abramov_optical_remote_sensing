# This code processes the binary result of ipmatch into a text file.
from lib_parse_match_file import * # process() function for the binary matches file of ipmatch
import argparse

parser = argparse.ArgumentParser(description='Convert ipmatch output into a text file.')
parser.add_argument('in_file_path', metavar='in_file_path', type=str, nargs=1,
                    help='Path to input binary file with matches')
parser.add_argument('out_file_path', metavar='out_file_path', type=str, nargs=1,
                    help='Path to output text file with matches')
args = parser.parse_args()

matches_bin_path = os.path.abspath(args.in_file_path[0])
out_path = os.path.abspath(args.out_file_path[0])

process(matches_bin_path, out_path)
