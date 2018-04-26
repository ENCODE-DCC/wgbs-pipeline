import glob
import argparse


def search(pattern):
    for file in glob.iglob(pattern, recursive=True):
        print(file)

def pattern_from_neighbour(nearest_file, pattern):
    subdirectory = nearest_file.rsplit('/', 1)[0]
    return subdirectory + '**/' + pattern

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print file path that matches pattern')
    parser.add_argument('--pattern', type=str, help='patter example:  subdirectory/**/*.gz', required=True)
    parser.add_argument('--nearest-neighbour', type=str, help='used for determining common parent directory')
    args = parser.parse_args()
    if args.nearest_neighbour:
        pattern = pattern_from_neighbour(args.nearest_neighbour, args.pattern)
    else:
        pattern = args.pattern
    search(pattern)
