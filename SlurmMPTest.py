import argparse
import random
import time
from joblib import Parallel, delayed

def run(file):
    sleep_time = random.randint(1, 10)
    time.sleep(sleep_time)
    print(f"Task for file {file} slept for {sleep_time} seconds and completed.")
    return file

if __name__ == "__main__":

    prog = "Slurm multi processing test script"
    usage = """
            python3.7 SlurmMPTest.py
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    args = parser.parse_args()

    files = [f"file_{i}.txt" for i in range(1, 51)]
    Parallel(n_jobs=args.threads)(delayed(run)(file) for file in files)
