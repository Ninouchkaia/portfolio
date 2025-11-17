import argparse
from src.pipeline import run_signatures, run_all


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=str, default="all")
    args = parser.parse_args()

    if args.step == "signatures":
        run_signatures()
    elif args.step == "all":
        run_all()
    else:
        print("Unknown step:", args.step)


if __name__ == "__main__":
    main()
