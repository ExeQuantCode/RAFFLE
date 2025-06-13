# This script tidies up the n-body database files by organizing them into subdirectories based on the number of bodies.
# By running this script, then the following agox line in the command line:
#    agox analysis method_ratios/*
# you will get analysis of the success split up by the method ratios used in the generation of the structures.

import os
import re
import shutil
from math import gcd
from functools import reduce

src_dir = "method_ratio"
dst_base = "method_ratios"

# Regex to match the filename pattern
pattern = re.compile(
    r"^db(?P<seed>\d+)_grow(?P<grow>[\d\.]+)_min(?P<min>[\d\.]+)_rand(?P<rand>[\d\.]+)_void(?P<void>[\d\.]+)_walk(?P<walk>[\d\.]+)\.db$"
)

def float_str_to_int_ratios(values):
    # Convert all to floats
    floats = [float(v) for v in values]
    scale = 1000  # To avoid float precision issues
    ints = [round(f * scale) for f in floats]
    divisor = reduce(gcd, ints)
    return [i // divisor for i in ints]

os.makedirs(dst_base, exist_ok=True)

for fname in os.listdir(src_dir):
    match = pattern.match(fname)
    if match:
        parts = match.groupdict()
        ratios = float_str_to_int_ratios([
            parts["void"], parts["rand"], parts["walk"], parts["grow"], parts["min"]
        ])
        # Assign to the appropriate method
        v, r, w, g, m = ratios
        # exit()
        subdir_name = f"v{v}_r{r}_w{w}_g{g}_m{m}"
        fname_new = f"db{parts['seed']}.db" #_{subdir_name}.db"
        subdir_path = os.path.join(dst_base, subdir_name)
        os.makedirs(subdir_path, exist_ok=True)

        src_path = os.path.join(src_dir, fname)
        dst_path = os.path.join(subdir_path, fname_new)
        shutil.copy2(src_path, dst_path)
