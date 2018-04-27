import glob
import argparse
import os


def search(pattern):
    for file in glob.iglob(pattern, recursive=True):
        print(file)

def pattern_from_current_dir(directory, pattern, up=2):
    # In Cromwell directory structure, sibling tasks are
    # 2 parent directory away
    search_directory = directory.rsplit('/', up)[0]
    return search_directory + '/**/' + pattern


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print file path that matches pattern')
    parser.add_argument('--pattern', type=str, help='patter example:  subdirectory/**/*.gz', required=True)
    parser.add_argument('--nearness', type=int, default=2, help='Number of directories to climb before searching ')
    args = parser.parse_args()
    current_dir = os.getcwd()
    pattern = pattern_from_current_dir(current_dir, args.pattern, up=args.nearness)
    search(pattern)
