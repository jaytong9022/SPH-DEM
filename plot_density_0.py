import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取 Fortran 输出的密度文件
df = pd.read_csv('density.csv')

# 绘制密度分布图
plt.figure(figsize=(6, 5))
sc = plt.scatter(df['x'], df['y'], c=df['rho'], cmap='viridis', s=10)
plt.colorbar(sc, label='Density (kg/m³)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('SPH Particle Density Field')
plt.axis('equal')
plt.tight_layout()
plt.show()
