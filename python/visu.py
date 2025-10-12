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
df = df[df['dim'] == dim]

df['length'] = df['death'] - df['birth']

finite_lengths = df.loc[df['death'] < float('inf'), 'length']
min_length = finite_lengths.min() if not finite_lengths.empty else 0

df = df[df['length'] > min_length + 1e-6]

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

if dim == 0:
    df0 = df.copy()
    max_death = df0['death'].replace(np.inf, 0).max() + 1
    df0['death_plot'] = np.where(np.isinf(df0['death']), max_death, df0['death'])

    events = []
    for _, row in df0.iterrows():
        events.append((row['birth'], +1))
        events.append((row['death_plot'], -1))

    events = []
    for _, row in df0.iterrows():
        if row['death'] < np.inf:
            events.append((row['death'], -1))

    events.sort()
    x_vals = [0]
    y_vals = [len(df0)]
    n_clusters = len(df0)
    for e, delta in events:
        n_clusters += delta
        x_vals.append(e)
        y_vals.append(n_clusters)

    plt.figure(figsize=(8,4))
    plt.step(x_vals, y_vals, where='post')
    plt.xlabel("Filtration value ε")
    plt.ylabel("Number of clusters")
    plt.title("Number of connected components vs ε (dim 0)")
    plt.show()
