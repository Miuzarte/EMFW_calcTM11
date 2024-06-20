import numpy as np
import matplotlib.pyplot as plt

# 长宽 间隔
a = 20 * 1e-3
b = 10 * 1e-3
h = 1 * 1e-3

# 误差阈值 最大迭代次数
tolerance = 1e-5
max_iterations = 1e5

# 创建矩阵
nx = int(a / h) + 1
ny = int(b / h) + 1
Ez = np.zeros((nx, ny))

# 内部初始条件
Ez[1:-1, 1:-1] = 1


# 更新Ez
def update_Ez(Ez, h, kc):
    rows, cols = Ez.shape
    Ez_new = np.zeros_like(Ez)
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            Ez_new[i, j] = (
                Ez[i + 1, j] +
                Ez[i - 1, j] +
                Ez[i, j + 1] +
                Ez[i, j - 1]
                ) / (4 - kc * kc * h * h)

    return Ez_new


# 更新kc
def update_Kc(Ez, h):
    kc_squared = np.zeros_like(Ez[1:-1, 1:-1])
    for i in range(1, Ez.shape[0] - 1):
        for j in range(1, Ez.shape[1] - 1):
            term = -(
                Ez[i + 1, j] +
                Ez[i - 1, j] +
                Ez[i, j + 1] +
                Ez[i, j - 1] -
                4 * Ez[i, j]
                ) / (Ez[i, j] * h * h)
            kc_squared[i - 1, j - 1] = term

    return np.sqrt(np.mean(kc_squared))


# 初始kc, 要大于0
kc = 0.1
kc_old = kc + 2*tolerance

n = 0
while True:
    n = n + 1
    if n >= max_iterations:
        break

    Ez_old = Ez.copy()
    Ez = update_Ez(Ez, h, kc)

    if n % 10 == 0:
        kc_old = kc
        kc = update_Kc(Ez, h)

    EzDiff = np.linalg.norm(Ez - Ez_old)
    KcDiff = np.linalg.norm(kc - kc_old)
    if EzDiff < tolerance or KcDiff < tolerance:
        break

# 理论值
kc_theory = np.sqrt((1 * np.pi / a) ** 2 + (1 * np.pi / b) ** 2)
lambda_c_theory = 2 * np.pi / kc_theory
f_c_theory = 3e8 / lambda_c_theory

# 计算值
kc = update_Kc(Ez, h)
lambda_c = 2 * np.pi / kc
f_c = 3e8 / lambda_c

print("Ez:")
r, c = Ez.shape
for i in range(r):
    for j in range(c):
        print("{:.6f}".format(Ez[i, j]), end=" ")
    print()
print("kc: {:.6f}".format(kc))
print("截止波长: {:.6f} mm".format(lambda_c))
print("截止频率: {:.6f} Hz".format(f_c))
print("理论kc: {:.6f}".format(kc_theory))
print("理论截止波长: {:.6f} mm".format(lambda_c_theory))
print("理论截止频率: {:.6f} Hz".format(f_c_theory))

Ez = Ez.T
X, Y = np.meshgrid(np.linspace(0, a, nx), np.linspace(0, b, ny))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Ez, cmap='viridis')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
plt.title('TM11 Ez')
plt.show()
