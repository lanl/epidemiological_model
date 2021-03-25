import argparse

my_parser = argparse.ArgumentParser()
my_parser.add_argument('-v',
                       '--verbosity-level',
                       action='store',
                       type=int)

args = my_parser.parse_args()
print(args.verbosity_level)
