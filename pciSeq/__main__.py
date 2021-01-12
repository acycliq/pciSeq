import argparse
from pciSeq.app import pciSeq, stage_data


def main(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', action="store", help='Configuration file')
    args = parser.parse_args()

    """The main routine."""
    # 1. prepare the data
    stage_data(args.config.PREPROCESS)

    # 2. cell typing
    pciSeq(args.config.MOUSE)  # 'MOUSE' or 'HUMAN'


if __name__ == "__main__":
    main()