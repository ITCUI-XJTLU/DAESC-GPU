import os
import pandas as pd
import numpy as np

def load_example_data():
    """
    Load the example dataset from the package.
    
    Returns:
        numpy.ndarray: The loaded example data as a NumPy array.
    """
    package_dir = os.path.dirname(__file__)
    data_path = os.path.join(package_dir, "data", "data_example.parquet")
    
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Dataset not found at {data_path}")
    
    example_df = pd.read_parquet(data_path)
    example_data = example_df.to_numpy()
    return example_data.astype(np.float32)
