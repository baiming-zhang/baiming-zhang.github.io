"""Shorten two geometry-panel titles; no native computation or metric changes."""
from pathlib import Path
import csv,json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import R,O,RAW,sha,save

assert not (O/'presentation_patch.json').exists()
before={f.name:sha(f) for f in O.glob('geometry.*')}
summary_hash=sha(O/'summary.json')
sites=list(csv.DictReader((RAW/'sites.csv').open(encoding='utf-8')));edges=list(csv.DictReader((RAW/'central_edges.csv').open(encoding='utf-8')))
plt.rcParams.update({'font.size':11,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(10,4.3),layout='constrained');s=[r for r in sites if r['case']=='short_face'];axes[0].scatter([float(r['x']) for r in s],[float(r['y']) for r in s],s=50,color='#3C5488')
for row in [r for r in edges if r['case']=='short_face']:
    x=np.array([float(row['ax']),float(row['bx'])]);y=np.array([float(row['ay']),float(row['by'])]);omitted=int(row['omitted_by_original_cutoff'])==1
    axes[0].plot(x,y,color='#E64B35' if omitted else '#3C5488',lw=2)
    axes[1].plot((x-.5)/8e-9,(y-.5)/8e-9,color='#E64B35' if omitted else '#3C5488',lw=3 if omitted else 1.5,marker='o',ms=5)
axes[0].set(xlabel='x / site spacing',ylabel='y / site spacing',title='Nine manufactured sites',aspect='equal',xlim=(-1.3,1.3),ylim=(-1.3,1.3))
axes[1].axhline(0,color='#999999',ls='--',lw=.7);axes[1].axvline(0,color='#999999',ls='--',lw=.7);axes[1].set(xlabel='(x - 0.5) / delta',ylabel='(y - 0.5) / delta',title='Corner zoom (delta = 8e-9)',aspect='equal',xlim=(-.8,.15),ylim=(-.8,.15))
for ext in ['png','pdf']:fig.savefig(O/('geometry.'+ext),dpi=200)
plt.close(fig)
assert sha(O/'summary.json')==summary_hash
shutil.copy2(__file__,O/'source-at-run/layout_fix.py')
save(O/'presentation_patch.json',dict(source_sha256=sha(Path(__file__)),before=before,after={f.name:sha(f) for f in O.glob('geometry.*')},summary_unchanged_sha256=summary_hash,new_native_executions=0,scope='Shortened clipped panel titles; same observed coordinates and explicit zoom, original report source retained'))
