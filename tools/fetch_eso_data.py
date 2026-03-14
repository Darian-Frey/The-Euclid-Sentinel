# tools/fetch_eso_data.py
import os
import time
from astroquery.mast import Observations

class EuclidDataFetcher:
    def __init__(self, download_dir="data/raw_fits"):
        self.download_dir = download_dir
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        print(f"[SENTINEL_LOG] Storage initialized at: {download_dir}")

    def fetch_bullet_cluster(self):
        print("[SENTINEL_LOG] STEP 1: Querying 1E 0657-558 records...")
        obs_table = Observations.query_object("1E 0657-558", radius=".02 deg")
        
        # Filter for HST (The most reliable lensing source for this cluster)
        hst_obs = obs_table[obs_table['obs_collection'] == 'HST']
        
        if len(hst_obs) == 0:
            print("[RECOVERY_NODE] No HST data found. Using first available record.")
            target = obs_table[0]
        else:
            # Pick the most recent science observation
            target = hst_obs[0]

        print(f"[SENTINEL_LOG] STEP 2: Fetching product list for ObsID: {target['obsid']}")
        products = Observations.get_product_list(target)
        
        print(f"[SENTINEL_LOG] STEP 3: Filtering for Science FITS (DRZ)...")
        # Explicitly filter the product table to avoid KeyError in download_products
        # We look for 'DRZ' in productSubGroupDescription or description
        mask = [('DRZ' in str(row['productSubGroupDescription']).upper()) for row in products]
        filtered_products = products[mask]

        if len(filtered_products) == 0:
            print("[RECOVERY_NODE] No DRZ found. Downloading first 2 science products.")
            filtered_products = products[:2]

        print(f"[SENTINEL_LOG] STEP 4: Initiating download of {len(filtered_products)} files...")
        manifest = Observations.download_products(filtered_products, 
                                                 download_dir=self.download_dir, 
                                                 curl_flag=False)
        return manifest

if __name__ == "__main__":
    fetcher = EuclidDataFetcher()
    manifest = fetcher.fetch_bullet_cluster()
    if manifest:
        print("\n[SUCCESS] Data ingested. Proceed to Mimetic Integration.")