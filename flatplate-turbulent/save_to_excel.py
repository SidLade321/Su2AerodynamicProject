import subprocess
import pandas as pd

command = 'SU2_CFD -d turb_SA_flatplate.cfg | grep -iE "skin|y_plus|pressure"'
result = subprocess.run(command, shell=True, capture_output=True, text=True)
lines = result.stdout.strip().split('\n')

parsed_data = []
for line in lines:
    if line.strip():
        parsed_data.append({"Terminal Output Log": line})

df = pd.DataFrame(parsed_data)
excel_file = "su2_simulation_data.xlsx"
df.to_excel(excel_file, index=False)
print(f"Success! Data exported to {excel_file}")
