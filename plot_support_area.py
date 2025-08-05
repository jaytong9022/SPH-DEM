import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load particle data
df = pd.read_csv('density.csv')  # Replace with the correct path if needed

# Parameters
target_id = 300
h = 0.024
h2 = h * h

# Check required columns
required_columns = {'id', 'x', 'y'}
if not required_columns.issubset(df.columns):
    raise ValueError(f"CSV file must contain columns: {required_columns}")

# Find target particle
target_particle = df[df['id'] == target_id]
if target_particle.empty:
    raise ValueError(f"Target particle with id={target_id} not found.")

x0, y0 = target_particle.iloc[0][['x', 'y']]

# Compute squared distances to all particles
df['r2'] = (df['x'] - x0)**2 + (df['y'] - y0)**2

# Find neighbors within support radius
neighbors = df[df['r2'] <= h2].copy()
neighbors['distance'] = np.sqrt(neighbors['r2'])

# Print ID and distance to command line
print(f"Neighbors within h={h} of particle id={target_id}:")
for _, row in neighbors.iterrows():
    pid = int(row['id'])
    dist = row['distance']
    print(f"  id={pid:<5} distance={dist:.6f}")

# Plot: only particles and support circle, no annotations
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(df['x'], df['y'], s=5, color='gray', label='All particles')
ax.scatter(neighbors['x'], neighbors['y'], s=20, color='blue', label='Neighbors')
ax.scatter([x0], [y0], s=50, color='red', label=f'Target id={target_id}')

# Draw support radius circle
circle = plt.Circle((x0, y0), h, color='red', fill=False, linestyle='--')
ax.add_patch(circle)

# Plot settings
ax.set_aspect('equal')
ax.set_title(f"Particles within support h={h} of id={target_id}")
ax.legend()
plt.tight_layout()
plt.show()
