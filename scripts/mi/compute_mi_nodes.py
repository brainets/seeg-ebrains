import os
import glob

import numpy as np
import xarray as xr

from frites.dataset import DatasetEphy
from frites.workflow import WfMi


if __name__ == '__main__':
    ###########################################################################
    # mi settings
    regr = 'difficulty'
    n_perm = 20
    inference = 'rfx'
    mcp = 'cluster'
    nb_min_suj = 10
    
    # file settings
    reference = 'sample_stim'  # {'sample_stim', 'sample_resp'}
    freq = "f50f150"           # {"f50f150", "f8f24"}
    smoothing = "sm0"

    # folders location
    analysis = '/hpc/brainets/data/db_ebrains/analysis'
    save_to = os.path.join(analysis, 'mi')
    ###########################################################################
    
    # __________________________________ I/O __________________________________
    # build the path where the epochs are saved
    path_epochs = os.path.join(
        analysis, 'epochs', reference, f"{freq}-{smoothing}", 'data')
    assert os.path.isdir(path_epochs)
    files = glob.glob(os.path.join(path_epochs, "*.nc"))

    # get mi_type
    mi_type = {
        'difficulty': 'cd', 'correct': 'cd', 'T_orientation': 'cd',
        'response_time': 'cc'
    }[regr]
    
    # define how to save the file
    if not os.path.isdir(save_to): os.makedirs(save_to)
    _save_as = (
        f"mi-{mi_type}-{regr}_{reference}-{freq}-{smoothing}_nperm-{n_perm}_"
        f"inference-{inference}_mcp-{mcp}.nc"
    )
    save_as = os.path.join(save_to, _save_as)
    print(save_as)

    # __________________________________ DATA _________________________________
    x = []
    for n_f, f in enumerate(files[0:5]):  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # show the file that is loaded
        _, fname = os.path.split(f)
        print(f"Loading {fname}", end='\r')
        
        # load the data
        _data = xr.load_dataarray(f)
        
        # make trials categorical (if needed)
        tr_before = _data[regr].data.copy()
        if mi_type == 'cd':
            trials = np.unique(_data[regr].data, return_inverse=True)[1]
            _data['trials'] = trials

        x += [_data.isel(roi=[r]) for r in range(len(_data['roi']))]
    
    # _______________________________-___ MI ____________-_____________________
    # define the dataset
    ds = DatasetEphy(x, y='trials', roi='roi', times='times',
                     nb_min_suj=nb_min_suj)

    # run the workflow
    wf = WfMi(mi_type=mi_type, inference=inference)
    mi, pv = wf.fit(ds, mcp=mcp, n_perm=n_perm)
    
    # save the results
    xr.Dataset({'mi': mi, 'pv': pv}).to_netcdf(save_as)
    