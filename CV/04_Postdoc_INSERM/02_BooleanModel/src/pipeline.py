from pathlib import Path
from .signatures import compute_signatures
from .utils import Timer

DATA = Path("data")


def run_signatures():
    input_file = DATA / "input" / "Dorothea_TF_activity_scale.tsv"
    output_file = DATA / "output" / "scale_rescaled_output.txt"

    with Timer() as t:
        results = compute_signatures(input_file, output_file)

    print("M1_signature:", results["M1_signature"])
    print("M2_signature:", results["M2_signature"])
    print("NLC_signature:", results["NLC_signature"])
    print("NLC count:", len(results["NLC_signature"]))
    print("Time:", t.duration)


def run_all():
    run_signatures()
