import pickle
import numpy as np
import scipy.io as sio

def convert_to_numpy(obj):
    """
    Recursively convert nested lists/tuples in obj to numpy arrays
    so that scipy.io.savemat can handle them.
    """
    if isinstance(obj, dict):
        return {k: convert_to_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        try:
            return np.array(obj)
        except Exception:
            return np.array([convert_to_numpy(el) for el in obj], dtype=object)
    elif isinstance(obj, np.ndarray):
        return obj
    else:
        return obj

# 1. Load the pickle file
with open('profs_196154_1600_py3.pkl', 'rb') as f:
    data = pickle.load(f)

# 2. Convert to numpy-friendly structure
data_mat = convert_to_numpy(data)

# 3. Save to MATLAB .mat
output_mat = 'profs_196154_1600_py3.mat'
sio.savemat(output_mat, {'data': data_mat}, oned_as='row')

print(f"Saved MATLAB file to: {output_mat}")
