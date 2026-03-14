# tools/euclid_loader.py
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import os

class EuclidLoader:
    def __init__(self, data_path="data/clusters/"):
        self.data_path = data_path
        if not os.path.exists(data_path):
            os.makedirs(data_path)

    def load_cluster_fits(self, filename):
        """
        Loads a FITS file (Baryonic mass map or Lensing map).
        Returns the data array and the Physical Header.
        """
        path = os.path.join(self.data_path, filename)
        with fits.open(path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header)
        
        # Clean NaNs (Common in real telescope data)
        data = np.nan_to_num(data, nan=0.0)
        return data, header, wcs

    def map_to_physical_units(self, data, header):
        """
        Euclid data is often in 'Flux' or 'Surface Brightness'.
        This maps pixel values to kg/m^3 based on redshift and distance.
        """
        # Placeholder for real Mass-to-Light ratio (M/L) calibration
        # For the Bullet Cluster (1E 0657-558), z=0.296
        mass_conversion_factor = 1.0e-15 
        return data * mass_conversion_factor

# --- SCRIPT TO PULL SAMPLE DATA ---
if __name__ == "__main__":
    print("EUCLID LOADER: Initialized.")
    print("READY TO INGEST FITS DATA FOR CLUSTER 1E 0657-558.")