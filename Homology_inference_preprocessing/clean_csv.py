import pandas as pd

# Load CSV with correct delimiter and encoding
df = pd.read_csv("../Species_names_IDs_wo_spaces.csv", delimiter=";", encoding="utf-8-sig")

# Clean all spaces and non-breaking spaces
df = df.applymap(lambda x: str(x).replace('\xa0', '').replace(' ', '') if isinstance(x, str) else x)
df.columns = [col.replace('\xa0', '').replace(' ', '') for col in df.columns]

# Save cleaned file
df.to_csv("Species_names_IDs_cleaned.csv", index=False)
