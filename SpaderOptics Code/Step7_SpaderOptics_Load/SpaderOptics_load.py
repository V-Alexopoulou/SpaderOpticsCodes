"""
Load the complete SpaderOptics database (-v7.3 .mat file) with h5py.

Result (kept in memory, needs about 3 GB):
    db[material][author][field]   -> data of every source
    spread[material][field]       -> source_spread of every material
                                     (None if a material has no spread)

Run with:   python3 -i SpaderOptics_load.py [path/to/SpaderOptics.mat]
(-i keeps Python open afterwards, so db and spread can be used.)
"""

import h5py
import numpy as np

import sys
path = sys.argv[1] if len(sys.argv) > 1 else 'SpaderOptics.mat'   # path to SpaderOptics.mat


def to_str(f, ref):
    """MATLAB text (stored as numbers) -> Python string."""
    return ''.join(chr(int(c)) for c in np.array(f[ref]).flatten())


def get(f, ref):
    """Value of one field of one source."""
    obj = f[ref]
    if 'MATLAB_object_decode' in obj.attrs:      # MATLAB 'string' objects cannot be read with h5py
        return None
    arr = np.array(obj)
    if arr.dtype == np.uint16:                   # MATLAB text
        return ''.join(chr(int(c)) for c in arr.flatten())
    return arr.T                                 # reverse axes: same shape as in MATLAB


def read_any(f, obj):
    """One field of the source_spread struct (numbers, text, or list of texts)."""
    arr = obj[()]
    if obj.dtype.kind == 'O':                    # cell array, e.g. the list of authors
        return [to_str(f, r) for r in obj[:].flatten()]
    if obj.dtype == np.uint16:                   # MATLAB text
        return ''.join(chr(int(c)) for c in np.array(arr).flatten())
    return np.array(arr).T                       # numbers: same shape as in MATLAB


db = {}
spread = {}

with h5py.File(path, 'r') as f:
    root = f['SpaderOptics']

    print('source_spread in file:', 'source_spread' in root)

    mat_refs = root['material'][:].flatten()
    ds_refs = root['dataset'][:].flatten()
    spread_refs = root['source_spread'][:].flatten() if 'source_spread' in root else None

    for m, (mref, dref) in enumerate(zip(mat_refs, ds_refs)):

        material = to_str(f, mref)
        ds = f[dref]
        fields = list(ds.keys())
        n_src = ds['author'].size

        # ---- sources of this material ----
        db[material] = {}
        for j in range(n_src):
            rec = {}
            for fld in fields:
                try:
                    rec[fld] = get(f, ds[fld][:].flatten()[j])
                except Exception:
                    rec[fld] = None              # e.g. MATLAB 'string' objects
            db[material][rec['author']] = rec

        # ---- source_spread of this material ----
        spread[material] = None
        if spread_refs is not None:
            try:
                sp = f[spread_refs[m]]
                if isinstance(sp, h5py.Group) and len(sp.keys()) > 0:
                    spread[material] = {k: read_any(f, sp[k]) for k in sp.keys()}
            except Exception as err:
                print(f'   could not read source_spread of {material}: {err}')

        print(f'{m + 1}/{len(mat_refs)}  {material}: {n_src} sources, '
              f'spread: {"yes" if spread[material] else "no"}')

print()
print('Materials:', list(db.keys()))
print('Al sources:', list(db['Al'].keys()))
print('E_hemispherical shape (Al, Rakic):',
      db['Al']['Rakic']['calculated_E_hemispherical'].shape)          # expected (50, 8)
if spread['Al']:
    print('Al spread fields:', list(spread['Al'].keys()))
    print('Al spread min shape:', spread['Al']['min'].shape)           # expected (50, 8)
