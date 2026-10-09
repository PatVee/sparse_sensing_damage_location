"""Build raw-event and reference fixtures for the shared MATLAB/Octave core."""
import argparse, importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import savemat

p=argparse.ArgumentParser(); p.add_argument('--python-results',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('pipeline',Path(__file__).with_name('vanersborg_full_pipeline_from_zip.py')); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
c=np.load(a.python_results/'00_aligned_opening_cache.npz'); X=c['aligned']; gauges=c['gauges']; files=c['files']; damage=np.arange(list(files).index(m.CRACK_FILE),len(files))
ri,energy,_=m.select_qr_references(X,gauges,m.FIT_IDX); wi=np.setdiff1d(np.arange(len(gauges)),ri); xv,_=m.fit_virtual_transfer(X,ri,wi,m.FIT_IDX)
d={'aligned':X,'gauges':gauges,'refIdx':ri+1,'energy':energy,'virtual':xv,'damageIdx':damage+1,'temp':c['temp']}
for k,(xx,gg) in enumerate([(X,gauges),(X[:,:,ri],gauges[ri]),(xv,gauges)],1):
 names,sc=m.zone_scores(xx,gg,m.FIT_IDX); d[f'scores{k}']=sc
 d[f'truthColumn{k}']=names.index(m.TRUTH_ZONE)+1
rng=np.random.default_rng(818)
for channels in [1,2,5,12]:
 arr=rng.normal(size=(256,channels)); d[f'featureInput{channels}']=arr; d[f'featureExpected{channels}']=m.area_gram_dct8(arr)
savemat(a.out/'reference.mat',d)
rawpaths={f.name:f for f in (a.python_results/'00_raw_extracted').rglob('*.csv')}
for i,f in enumerate(files):
 frame=pd.read_csv(rawpaths[f],usecols=m.USECOLS)
 raw=frame[m.USECOLS].to_numpy(dtype=float)
 savemat(a.out/f'raw_{i+1:02}.mat',{'raw':raw})
 if (i+1)%8==0:print(f'Prepared {i+1}/64 raw numeric fixtures',flush=True)
print('Prepared fixtures in',a.out)
