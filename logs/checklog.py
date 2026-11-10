import os
import sys
import argparse


# Set directories for logging files
BASE_DIR = "/data/ZVAR-Forced-Phot/logs"


parser = argparse.ArgumentParser(description='Provide input log file.')
parser.add_argument(
    '-l',
    '--logfile',
    type=str,
    default=f'{BASE_DIR}/log_files/pipeline.log',
    help="Pipeline log filename."
)
args = parser.parse_args()
filename = args.logfile
if not os.path.isfile(filename):
    sys.exit(1)

fields = []
ccds = []
quads = []
filters = []
with open(filename) as file:

    for line in file.readlines():
        if "FAILED" in line:
            quadID = line.strip("\n").split(":")[-1].split(" ")[0]
        elif "MetricError" in line:
            quadID = line.strip("\n").split(":")[-1].split(" ")[0]
        else:
            continue

        fields.append(int(quadID.split("_")[0].strip("f")))
        ccds.append(int(quadID.split("_")[1].strip("c")))
        quads.append(int(quadID.split("_")[2].strip("q")))
        filters.append(quadID.split("_")[3].strip("z"))

Nfailed = len(fields)


# Save failed quadrants to file
savename = f"{BASE_DIR}/failed_quads.txt"
for i in range(Nfailed):
    line_text = f"{fields[i]} {ccds[i]} {quads[i]} {filters[i]}"
    with open(savename,'a') as file:
        file.write(f"{line_text}\n")
