# This script tidies up the n-body database files by organizing them into subdirectories based on the number of bodies.
# By running this script, then the following agox line in the command line:
#    agox analysis n-bodies/*
# you will get analysis of the n-body successes split up by the number of contributing n-body distribution functions.

import os
import re
import shutil
from math import gcd
from functools import reduce

src_dir = "n-body"
dst_base = "n-bodies"

# Regex to match the filename pattern
pattern = re.compile(
    r"^db(?P<seed>\d+)_dist_(?P<nbody>[\d\.]+)\.db$"
)

os.makedirs(dst_base, exist_ok=True)

for fname in os.listdir(src_dir):
    match = pattern.match(fname)
    if match:
        parts = match.groupdict()
        subdir_name = f"{parts['nbody']}-body"
        fname_new = f"db{parts['seed']}.db" #_{subdir_name}.db"
        subdir_path = os.path.join(dst_base, subdir_name)
        os.makedirs(subdir_path, exist_ok=True)

        src_path = os.path.join(src_dir, fname)
        dst_path = os.path.join(subdir_path, fname_new)
        shutil.copy2(src_path, dst_path)
