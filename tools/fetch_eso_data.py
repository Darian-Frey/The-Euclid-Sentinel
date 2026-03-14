# tools/fetch_eso_data.py
import os
import time
from astroquery.mast import Observations

class EuclidDataFetcher:
    def __init__(self, download_dir="data/raw_fits"):
        self.download_dir = download_dir
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)
        
        # TARGET REGISTRY: Coordinates and preferred datasets
        self.registry = {
            "bullet": {"name": "1E 0657-558", "radius": ".02 deg"},
            "abell370": {"name": "Abell 370", "radius": ".03 deg"},
            "elgordo": {"name": "ACT-CL J0102-4915", "radius": ".04 deg"}
        }
        print(f"[SENTINEL_LOG] Survey Engine initialized. Storage: {download_dir}")

    def fetch_target(self, target_key="bullet"):
        if target_key not in self.registry:
            print(f"[ERROR] Target {target_key} not in registry.")
            return None

        target_info = self.registry[target_key]
        print(f"\n[SENTINEL_LOG] INITIATING FETCH: {target_info['name']}")
        
        # 1. Query MAST
        obs_table = Observations.query_object(target_info['name'], radius=target_info['radius'])
        
        # 2. Filter for HST (The Source of Truth for lensing)
        hst_obs = obs_table[obs_table['obs_collection'] == 'HST']
        
        if len(hst_obs) == 0:
            print(f"[RECOVERY_NODE] No HST data found for {target_key}. Aborting.")
            return None

        # Select the most substantial observation (often the one with the most products)
        target = hst_obs[0]
        print(f"[SENTINEL_LOG] ObsID Locked: {target['obsid']} ({target['project']})")

        # 3. Get Products
        products = Observations.get_product_list(target)
        
        # 4. Filter for DRZ (Drizzled Science Images)
        # We prioritize ACS/WFC for lensing analysis
        mask = [('DRZ' in str(row['productSubGroupDescription']).upper()) for row in products]
        filtered_products = products[mask]

        if len(filtered_products) == 0:
            print("[RECOVERY_NODE] No DRZ found. Falling back to science products.")
            filtered_products = products[products['productType'] == 'SCIENCE'][:1]
        else:
            # We only need the primary drizzled FITS to save your SSD space
            filtered_products = filtered_products[:1]

        print(f"[SENTINEL_LOG] Downloading calibrated FITS: {filtered_products[0]['productFilename']}")
        manifest = Observations.download_products(filtered_products, 
                                                 download_dir=self.download_dir, 
                                                 curl_flag=False)
        return manifest

if __name__ == "__main__":
    import sys
    fetcher = EuclidDataFetcher()
    
    # Check if user provided a target via CLI, else default to Abell 370
    target = sys.argv[1] if len(sys.argv) > 1 else "abell370"
    
    manifest = fetcher.fetch_target(target)
    if manifest:
        print(f"\n[SUCCESS] {target.upper()} ingested. Ready for Mimetic Integration.")