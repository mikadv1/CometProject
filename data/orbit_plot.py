import numpy as np
import matplotlib.pyplot as plt

# Загрузка данных из CSV (пропускаем заголовок)
data = np.loadtxt("orbit_output.csv", delimiter=",", skiprows=1)

x = data[:, 1]  # колонка X
y = data[:, 2]  # колонка Y

# Построение графика
plt.figure(figsize=(8, 8))
plt.plot(x, y, 'b-', linewidth=1, label='3I/ATLAS')
plt.plot(0, 0, 'o', color='orange', markersize=15, label='Солнце')
plt.plot(x[0], y[0], 'go', markersize=3, label='Начальное положение')
plt.plot(x[-1], y[-1], 'ro', markersize=3, label='Конечное положение')

plt.xlabel('X (а.е.)')
plt.ylabel('Y (а.е.)')
plt.title('Траектория кометы 3I/ATLAS (плоскость XY)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.axis('equal')

plt.savefig('orbit_xy.png', dpi=300)
plt.show()
