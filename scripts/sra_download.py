# (Optional) For advanced SRA metadata handling
import pandas as pd
from Bio import Entrez

def fetch_sra_ids(gse_id: str):
    Entrez.email = "your@email.com"
    handle = Entrez.esearch(db="sra", term=f"{gse_id}[GSE]")
    return Entrez.read(handle)["IdList"]

# Example: Add to geo_fetch.py to auto-populate config["sra_ids"]
