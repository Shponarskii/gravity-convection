#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def load_field(filename):

    raw = []

    xs = []
    ys = []

    with open(filename, 'r') as f:
        for line in f:
            x, y, val = line.strip().split(',')
            x = float(x)
            y = float(y)
            val = float(val)

            raw.append((x, y, val))
            xs.append(x)
            ys.append(y)

    xs = np.array(xs)
    ys = np.array(ys)

    # --- определяем сетку ---
    x_unique = np.unique(xs)
    y_unique = np.unique(ys)

    M = len(x_unique) - 1
    N = len(y_unique) - 1

    Hx = x_unique.max() - x_unique.min()
    Hy = y_unique.max() - y_unique.min()

    hx = Hx / M if M > 0 else 1
    hy = Hy / N if N > 0 else 1

    # --- поле ---
    field = np.zeros((M + 1, N + 1))

    # mapping через ближайший индекс
    for x, y, val in raw:
        i = int(round((x - x_unique.min()) / hx))
        j = int(round((y - y_unique.min()) / hy))
        field[i, j] = val

    return field, x_unique.min(), x_unique.max(), y_unique.min(), y_unique.max(), M, N

base_path = "../Гравитационная конвекция/"
C, x0, x1, y0, y1, M, N = load_field(base_path + "C.csv")
T, _, _, _, _, _, _ = load_field(base_path + "T.csv")
Vx, _, _, _, _, _, _ = load_field(base_path + "Vx.csv")
Vy, _, _, _, _, _, _ = load_field(base_path + "Vy.csv")

x = np.linspace(x0, x1, M + 1)
y = np.linspace(y0, y1, N + 1)

X, Y = np.meshgrid(x, y, indexing='ij')
plt.close('all')

# Поле концентрации

fig1 = plt.figure(figsize=(12, 8))
ax1 = fig1.add_subplot(111, projection='3d')

surf1 = ax1.plot_surface(X, Y, C, cmap='jet')
fig1.colorbar(surf1, shrink=0.5, aspect=10)

ax1.set_title('Поле концентрации')

plt.show(block=False)

# Поле температуры

fig2 = plt.figure(figsize=(12, 8))
T_smooth = gaussian_filter(T, sigma=8.0)

ax2 = fig2.add_subplot(111, projection='3d')

surf2 = ax2.plot_surface(X, Y, T_smooth, cmap='jet', antialiased=True)
fig2.colorbar(surf2, shrink=0.5, aspect=10)

ax2.set_title('Поле температуры')

plt.show(block=False)

# Поле скоростей

fig3 = plt.figure(figsize=(12,8))
Vmag = np.sqrt(Vx**2 + Vy**2)
contour = plt.contourf(
    X,
    Y,
    Vmag,
    levels=50,
    cmap='jet'
)
cbar = plt.colorbar(contour)
cbar.set_label('Величина скорости')
stream = plt.streamplot(
    x, y,
    Vx.T,
    Vy.T,
    density=2,
    linewidth=1.0,
    color='black'
)
plt.axis("equal")
plt.title("Распределение скоростей")
plt.show()

# In[ ]:




