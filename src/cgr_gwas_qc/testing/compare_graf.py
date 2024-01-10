import pandas as pd 
import matplotlib.pyplot as plt

df1 = pd.read_csv('~/GwasQcPipeline/with_modifications_run/sample_level/ancestry/graf_populations.txt', delimiter='\t')
df2 = pd.read_csv('~/GwasQcPipeline/with_modifications_run/sample_level/call_rate_2/testing.txt', delimiter='\t', skiprows=7)

merged = pd.merge(df1, df2, on='Sample', how='inner')

print(merged.columns)
plt.scatter(merged['E(%)_x'], merged['E(%)_y'], label='Variable 1')
plt.xlabel('Old European %')
plt.ylabel('New European %')
plt.title('Comparing Old vs New Ancestry Results')
plt.plot(merged['E(%)_x'], merged['E(%)_x'], color='red', linestyle='--', label='y=x Line')

plt.savefig("compare.png")
