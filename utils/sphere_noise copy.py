import pyshtools as pysh
import numpy as np
import trimesh
from scipy.special import sph_harm
import matplotlib.pyplot as plt
import cmath



def icosphere():
    sphere_cvt = trimesh.load('/media/box/Elements/Exp/HCPSOTS/sphere_cvt.obj')
    pts_xyz = np.array(sphere_cvt.vertices)  # radius = 1

    return pts_xyz
    # print(points_xyz.shape)
    # N = points_xyz.shape[0]

    # x = points_xyz[:, 0]
    # y = points_xyz[:, 1]
    # z = points_xyz[:, 2]

    # # convert to spherical coordinates
    # phi = np.arctan2(y, x)  # in radians
    # phi = np.mod(phi + 2 * np.pi, 2 * np.pi)   # [0, 2pi]
    # theta = np.arccos(z)   # in radians, [-pi/2, pi/2]

    # lonlat = np.stack([phi, theta], axis=1)  # (N, 2)

    # return lonlat

def random_on_sphere(N=3000):
    pts_xyz = np.random.uniform(-1, 1, (N, 3))
    pts_xyz = pts_xyz / np.linalg.norm(pts_xyz, axis=1, keepdims=True)

    return pts_xyz

def generate_fibonacci_points(N):
    # 黄金角（以弧度为单位）
    golden_angle = np.pi * (3 - np.sqrt(5))
    # 生成点
    pts_xyz = np.zeros((N, 3))
    for i in range(N):
        # 计算纬度和经度
        lat = np.arcsin(1 - 2 * (i + 0.5) / N)
        lon = golden_angle * i
        # 转换为笛卡尔坐标
        x = np.cos(lon) * np.cos(lat)
        y = np.sin(lon) * np.cos(lat)
        z = np.sin(lat)
        pts_xyz[i] = [x, y, z]
    return pts_xyz


# 将球面上的点转换为经纬度
def points_to_lat_lon(points):
    lats = np.degrees(np.arcsin(points[:, 2]))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    lons = np.mod(lons + 2 * np.pi, 2 * np.pi)   # [0, 2pi]
    return lats, lons

# 将数据转换为网格
def data_to_grid(lats, lons, nlat=800, nlon=1600):
    grid = np.zeros((nlat, nlon))
    for lat, lon in zip(lats, lons):
        ilat = int((lat + 90) / 180 * (nlat - 1))
        ilon = int(lon / 360 * (nlon - 1))
        grid[ilat, ilon] = 20.0
    return grid

def get_spectrum(pts_xyz, lmax):
    lats, lons = points_to_lat_lon(pts_xyz)
    grid = data_to_grid(lats, lons)
    grid = pysh.SHGrid.from_array(grid, grid='DH')

    # 计算球谐变换
    lmax = 256  # 最大球谐阶数
    coeffs = grid.expand(lmax_calc=lmax, normalization='schmidt')
    # 创建一个SHT对象
    clm = pysh.SHCoeffs.from_array(coeffs.coeffs, lmax=lmax, normalization='4pi', csphase=1)
    power_spectrum = clm.spectrum()

    return power_spectrum

# lonlat = icosphere()
pts_xyz_uniform = random_on_sphere(N=9000)
pts_xyz_icon = icosphere()
pts_xyz_Fibo = generate_fibonacci_points(N=10000)


# 计算球谐变换
lmax = 256  # 最大球谐阶数

power_spectrum_uniform = get_spectrum(pts_xyz_uniform, lmax)
power_spectrum_icon = get_spectrum(pts_xyz_icon, lmax)
power_spectrum_Fibo = get_spectrum(pts_xyz_Fibo, lmax)

degrees = np.arange(lmax+1)
plt.figure(figsize=(10, 7))
plt.plot(degrees, power_spectrum_uniform, linestyle='-', color='tomato', label='Uniform')
plt.plot(degrees, power_spectrum_icon, linestyle='-', color='deeppink', label='Icon')
plt.plot(degrees, power_spectrum_Fibo, linestyle='-', color='forestgreen', label='Fibo')
# plt.title('Angular Power Spectrum')
plt.legend()
plt.xlabel('Degree')
plt.ylabel('Power')
plt.grid(True, which='both', ls='-.')
plt.show()


# N = lonlat.shape[0]
# data = np.zeros((N, N))
# for pt in lonlat:
#     coords = int(pt[0]), int(pt[1])
#     data[coords] = 1e6*np.random.uniform()

# grid = pysh.SHGrid.from_array(data, grid='DH')

# # lmax 是球谐变换的最大阶数，它决定了频谱的分辨率
# lmax = 1096

# # 将网格数据转换为球谐系数
# coeffs = grid.expand(lmax_calc=lmax, normalization='schmidt')
# print(coeffs.coeffs.shape)


# # # 现在 coeffs 包含了球面上数据点的频谱
# # # coeffs.coeffs 是一个二维数组，其中每一行对应一个球谐阶数 l，每一列对应一个球谐阶数 m
# # # coeffs.coeffs[l, m] 就是对应于 l, m 的球谐系数

# # # # 打印频谱信息
# # # for l in range(lmax+1):
# # #     for m in range(-l, l+1):
# # #         print(f"l={l}, m={m}: {coeffs.coeffs[0,l, m+l]}")

# # # 用球谐变换计算球面上的数据点的复数球谐系数,并绘制出角功率频谱图

# # 计算角功率频谱
# power_spectrum = coeffs.spectrum()

# fig, ax = plt.subplots(figsize=(10, 6))
# plt.plot(np.arange(lmax + 1), power_spectrum, linestyle='-', color='b')
# plt.xlabel('Normalized Frequencies')
# plt.ylabel('Angular Power Spectrum')
# ax.set_ylim(bottom=0)
# # plt.title('Angular Power Spectrum')
# plt.grid(True, which='both', ls='-.')
# plt.show()