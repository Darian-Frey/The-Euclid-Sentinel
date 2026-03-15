# tools/fetch_euclid_data.py
import os
from astroquery.esa.euclid import Euclid

class EuclidMissionFetcher:
    def __init__(self, download_dir="data/euclid_raw"):
        self.download_dir = download_dir
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        print(f"[SENTINEL_LOG] ESA TAP Interface initialized: {download_dir}")

    def fetch_perseus_ero(self):
        print("[SENTINEL_LOG] Probing ESA Archive Schema...")
        
        # 1. Standardized query for public Euclid observations
        # We use the 'ivoa.obscore' table which is universal across ESA archives
        query = """
        SELECT TOP 1 observation_id, access_url, target_name
        FROM ivoa.obscore
        WHERE (target_name LIKE '%Perseus%' OR target_name LIKE '%NGC1275%')
        AND instrument_name = 'VIS'
        AND dataproduct_type = 'image'
        ORDER BY t_min DESC
        """
        
        try:
            print("[SENTINEL_LOG] Executing ADQL Search for Perseus VIS Stacks...")
            job = Euclid.launch_job(query)
            results = job.get_results()
            
            if len(results) == 0:
                print("[RECOVERY_NODE] No Perseus data found in obscore. Attempting schema list...")
                tables = Euclid.load_tables(only_names=True)
                print(f"[SENTINEL_LOG] Available Tables: {tables[:5]}")
                return None

            obs_id = results[0]['observation_id']
            print(f"[SENTINEL_LOG] Target Locked: {obs_id} (Target: {results[0]['target_name']})")
            
            # 2. Downloading the primary VIS Stack
            # We use the standardized product retrieval for VIS images
            print(f"[SENTINEL_LOG] Downloading FITS product via Datalink...")
            manifest = Euclid.get_products(observation_id=obs_id, 
                                          download_dir=self.download_dir,
                                          product_type='SCIENCE_IMAGE')
            return manifest
            
        except Exception as e:
            print(f"[RECOVERY_NODE] TAP Query Failed: {e}")
            return None

if __name__ == "__main__":
    fetcher = EuclidMissionFetcher()
    manifest = fetcher.fetch_perseus_ero()
    if manifest:
        print("\n[SUCCESS] Euclid Data Ingested. Ready for Mimetic-Euclid Calibration.")