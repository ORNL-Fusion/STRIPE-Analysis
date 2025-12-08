import pickle
import os

import numpy as np
import matplotlib.pyplot as plt

def inspect_dict(d, indent=0):
    """Recursively print keys and types of nested dicts."""
    prefix = "  " * indent
    for k, v in d.items():
        if isinstance(v, dict):
            print(f"{prefix}- {k}: dict")
            inspect_dict(v, indent + 1)
        elif isinstance(v, (list, tuple, np.ndarray)):
            shape = np.shape(v)
            print(f"{prefix}- {k}: {type(v).__name__}, shape={shape}")
        else:
            print(f"{prefix}- {k}: {type(v).__name__}")

def main():
    # Adjust this path if needed
    file_path = "profs_200882_1800_py3.pkl"

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    # Load the pickle file
    with open(file_path, "rb") as f:
        data = pickle.load(f)

    # List all top-level keys
    print("Top-level keys:")
    for key in data.keys():
        print(f" - {key}")

    # Inspect each top-level item
    print("\nDetails of each item:")
    for key, value in data.items():
        print(f"\nItem '{key}':")
        if isinstance(value, dict):
            inspect_dict(value, indent=1)
        else:
            print(f"  Type: {type(value).__name__}")

    # Example: Try plotting from a nested dict
    # For illustration, attempt to find and plot 'y' data
    # Adjust this to the actual nested keys once you see them printed
    sample_key = "nedatpsi"
    if sample_key in data and isinstance(data[sample_key], dict):
        subdict = data[sample_key]
        # Example: try to get the first array-like item
        for subkey, subval in subdict.items():
            if isinstance(subval, (list, tuple, np.ndarray)):
                y = np.array(subval)
                x = np.arange(len(y))
                plt.figure()
                plt.plot(x, y, marker="o")
                plt.xlabel("Index")
                plt.ylabel(f"{sample_key} / {subkey}")
                plt.title(f"Profile: {sample_key} -> {subkey}")
                plt.grid(True)
                plt.tight_layout()
                plt.show()
                break
        else:
            print(f"No array-like data found in '{sample_key}' to plot.")
    else:
        print(f"'{sample_key}' is not a dict or not present.")

if __name__ == "__main__":
    main()
