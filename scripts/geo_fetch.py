import pandas as pd
from Bio import Entrez
import yaml

# Load GEO query terms from config
with open("config/geo_query.yaml") as f:
    geo_config = yaml.safe_load(f)

# Extract: Fetch GEO metadata
Entrez.email = geo_config["email"]
handle = Entrez.esearch(db="gds", term=geo_config["query"])
record = Entrez.read(handle)

# Transform: Clean data
df = pd.DataFrame({
    "GSE_ID": record["IdList"],
    "Title": [Entrez.read(Entrez.esummary(db="gds", id=id))[0]["title"] for id in record["IdList"]]
})

# Load: Save CSV
df.to_csv("results/geo/geo_metadata.csv", index=False)
