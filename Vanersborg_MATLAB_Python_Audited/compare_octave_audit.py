"""Compare all Octave partition outcomes and references with Python results."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat
p=argparse.ArgumentParser();p.add_argument('--python-results',type=Path,required=True);p.add_argument('--octave-results',type=Path,required=True);a=p.parse_args()
r=loadmat(a.octave_results/'octave_audit_results.mat')
f=pd.read_csv(a.python_results/'05_baseline_split_sensitivity.csv')
expected=f[['detected_of_7','correct_zone_of_7','joint_of_7']].to_numpy()
assert np.array_equal(r['splitCounts'],expected),'Partition counts differ'
refs=f[f.configuration=='Reduced QR physical sensors'].references
expected_refs=np.array([[int(g.removeprefix('SG')) for g in s.split(',')] for s in refs])
assert np.array_equal(r['splitRefs'],expected_refs),'Partition QR references differ'
print('PASS: all 60 configuration/partition rows and all 20 QR reference sets match Python.')
