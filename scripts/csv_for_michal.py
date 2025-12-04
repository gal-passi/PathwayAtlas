import pandas as pd
import os
import sys

# Standardize path setup
try:
    back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(back_dir)
except NameError:
    print("Running in an environment where __file__ is not defined. Assuming project root is correctly configured.")

from utils import KeggApi






# Goal of this script:
#       get all for each cancer his description, and all the modules and pathways of him in the results:
#       ./results_and_graphs/scores/clinvar_reg_dis_ordered_prob-kl_divergence
#       and store them nicely in tables, to show Michal





# ---------------------------------------------------------
# 1. SETUP PATHS AND CONFIGURATION
# ---------------------------------------------------------

# Input Paths
SUMMARY_CSV_PATH = "/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/results_summary.csv"
MUTATIONS_DIR = "/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/scores/clinvar_reg_dis_ordered_prob-kl_divergence"

# Output Path
OUTPUT_FILE = "/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/pathway_analysis_combined_for_michal.xlsx"

# ---------------------------------------------------------
# 2. INITIALIZE KEGG API & DICTIONARIES
# ---------------------------------------------------------

# Global dictionary to hold both Pathways and Modules
global_knowledge_dict = {}


def fetch_modules_direct():
    """
    Fetches the list of KEGG modules directly from the REST API.
    Used to ensure module support (M00001, etc.) even if KeggApi class
    doesn't have a specific method for it.
    """
    print("Fetching KEGG Modules list (direct from rest.kegg.jp)...")
    try:
        # KEGG returns list as: "md:M00001 \t Description"
        url = "http://rest.kegg.jp/list/module"
        df = pd.read_csv(url, sep="\t", header=None, names=["id", "desc"])

        modules_map = {}
        for _, row in df.iterrows():
            full_id = str(row['id'])  # e.g. "md:M00001"
            desc = str(row['desc'])

            # Store exact match
            modules_map[full_id] = desc

            # Store stripped match (M00001)
            if full_id.startswith("md:"):
                modules_map[full_id.replace("md:", "")] = desc

        print(f"Fetched {len(modules_map)} module entries.")
        return modules_map
    except Exception as e:
        print(f"Warning: Could not fetch modules directly: {e}")
        return {}


try:
    # 1. Initialize API
    kegg = KeggApi()

    # 2. Fetch Pathways (hsa...)
    print("Fetching KEGG pathways dictionary...")
    pathways_dict = kegg.get_all_pathways()
    print("Fetching pathways done.")

    # 3. Fetch Modules (M00001...)
    # We try to use the class method if it exists, otherwise fallback to direct fetch
    try:
        if hasattr(kegg, 'get_all_modules'):
            modules_dict = kegg.get_all_modules()
        else:
            modules_dict = fetch_modules_direct()
    except Exception as e:
        print(f"Error checking for modules: {e}")
        modules_dict = fetch_modules_direct()

    # 4. Merge dictionaries
    # This creates one master lookup for both Pathways and Modules
    global_knowledge_dict = {**pathways_dict, **modules_dict}

except NameError:
    print("Warning: 'KeggApi' class not found. Creating an empty dictionary.")
    global_knowledge_dict = {}
except Exception as e:
    print(f"Error initializing KEGG data: {e}")
    global_knowledge_dict = {}


def get_description(item_id):
    """
    Looks up the info for a Pathway OR a Module.
    Handles variations: 'hsa001', 'path:hsa001', 'M00001', 'md:M00001'
    """
    item_id = str(item_id).strip()

    # 1. Try exact match
    if item_id in global_knowledge_dict:
        return global_knowledge_dict[item_id]

    # 2. Try common KEGG prefixes
    # Pathways often use 'path:' or 'pathway:'
    # Modules often use 'md:' or 'module:'
    prefixes = ["path:", "pathway:", "md:", "module:", "hsa:", "map:"]

    for prefix in prefixes:
        prefixed_id = f"{prefix}{item_id}"
        if prefixed_id in global_knowledge_dict:
            return global_knowledge_dict[prefixed_id]

    # 3. Handle cases where the ID in CSV has a prefix, but dict key doesn't
    # (e.g. input is "path:hsa001" but dict has "hsa001")
    if ":" in item_id:
        stripped_id = item_id.split(":")[-1]
        if stripped_id in global_knowledge_dict:
            return global_knowledge_dict[stripped_id]

    return "Unknown / Not Found"


# ---------------------------------------------------------
# 3. MAIN PROCESSING
# ---------------------------------------------------------

def main():
    # Load the summary file
    if not os.path.exists(SUMMARY_CSV_PATH):
        print(f"Error: Summary file not found at {SUMMARY_CSV_PATH}")
        return

    summary_df = pd.read_csv(SUMMARY_CSV_PATH)

    # Create an Excel Writer
    print(f"Creating Excel file: {OUTPUT_FILE}")

    # Note: Using 'xlsxwriter' as requested. If not installed, pip install XlsxWriter
    with pd.ExcelWriter(OUTPUT_FILE, engine='xlsxwriter') as writer:

        for index, row in summary_df.iterrows():
            short_name = row['cancer_short_name']
            cancer_type = row['cancer_type']

            # Construct path to the specific mutation file
            mutation_file_path = os.path.join(MUTATIONS_DIR, f"{short_name}_mutations.csv")

            if not os.path.exists(mutation_file_path):
                print(f"Skipping {short_name}: File not found ({mutation_file_path})")
                continue

            print(f"Processing: {cancer_type} ({short_name})...")

            # Read the mutation CSV
            df = pd.read_csv(mutation_file_path)

            # 1. Add pathway_info column
            # Use the unified get_description function
            df['pathway_info'] = df['pathway_name'].apply(get_description)

            # 2. Reorder columns
            # Ensure 'pathway_info' is right after 'pathway_name'
            cols = list(df.columns)
            if 'pathway_info' in cols:
                cols.remove('pathway_info')

            # Find index of pathway_name to insert after it
            try:
                target_index = cols.index('pathway_name') + 1
                cols.insert(target_index, 'pathway_info')
            except ValueError:
                # If pathway_name missing, just append
                cols.append('pathway_info')

            df = df[cols]

            # 3. Sort the data
            # Primary: q_value (Ascending / Lower is better)
            # Secondary: distance (Descending / Higher is better)
            if 'q_value' in df.columns and 'distance' in df.columns:
                df = df.sort_values(by=['q_value', 'distance'], ascending=[True, False])

            # 4. Write to Excel Sheet
            # Excel sheet names are limited to 31 chars. We must truncate if necessary.
            safe_sheet_name = str(cancer_type)[:31].replace(":", "").replace("/", "-")

            df.to_excel(writer, sheet_name=safe_sheet_name, index=False)

    print("Done! File saved successfully.")


if __name__ == "__main__":
    main()