import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} <csv_file> <dim>")
    sys.exit(1)

filename = sys.argv[1]
df = pd.read_csv(filename)
dim = int(sys.argv[2])

plt.figure(figsize=(6,6))
for _, row in df.iterrows():
    if row['death'] == float('inf'):
        plt.scatter(row['birth'], 5, marker='^', label=f"dim {int(row['dim'])}")
    else:
        plt.scatter(row['birth'], row['death'], label=f"dim {int(row['dim'])}")
plt.plot([0, df[['birth','death']].replace("inf", np.nan).max().max()],
         [0, df[['birth','death']].replace("inf", np.nan).max().max()],
         'k--')
plt.xlabel("Birth")
plt.ylabel("Death")
plt.title("Persistence Diagram")
plt.show()

plt.figure(figsize=(8,4))
y = 0
for _, row in df.iterrows():
    if row['dim'] != dim:
        continue
    if row['death'] == float('inf'):
        plt.hlines(y, row['birth'], df[['birth','death']].replace("inf", np.nan).max().max(), colors='r')
    else:
        plt.hlines(y, row['birth'], row['death'])
    y += 1
plt.xlabel("Filtration value")
plt.ylabel("Intervals")
plt.title("Persistence Barcode")
plt.show()
