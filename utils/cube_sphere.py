import os
import numpy as np
import trimesh
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams

import pyvista as pv

out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'

sphere = pv.Sphere(radius=1.0, center=(0, 0, 0), direction=(0, 0, 1), theta_resolution=179, phi_resolution=80)

# mask = sphere.points[:, 2] > 0
# hemisphere = sphere.extract_points(mask, adjacent_cells=True)

# mask = (sphere.points[:, 0] >= 0) & (sphere.points[:, 1] >= 0) & (sphere.points[:, 2] >= 0)


octant = sphere.clip('x',origin=(0,0,0),invert=False).clip('y',origin=(0,0,0),invert=False).clip('z',origin=(0,0,0),invert=False)

png_file = out_dir + f'cube_sphere_octant.png'

plane1 = pv.Plane(center=(0.5, 0.5, 1.0), direction=(0, 0, 1), i_size=1, j_size=1)
plane2 = pv.Plane(center=(1.0, 0.5, 0.5), direction=(1, 0, 0), i_size=1, j_size=1)
plane3 = pv.Plane(center=(0.5, 1.0, 0.5), direction=(0, 1, 0), i_size=1, j_size=1)

# p = pv.Plotter(off_screen=True)
p = pv.Plotter(window_size=(1080, 900))
p.add_mesh(octant, color=(1.0, 0.376, 0.91), smooth_shading=True, opacity=0.5)
p.add_mesh(plane1, color=(0.275966, 0.792482, 0.8), opacity=0.5, smooth_shading=True)
p.add_mesh(plane2, color=(0.275966, 0.792482, 0.8), opacity=0.5, smooth_shading=True)
p.add_mesh(plane3, color=(0.275966, 0.792482, 0.8), opacity=0.5, smooth_shading=True)
p.show_axes()


# one octant of the cube and sphere
vtx = np.array([[1.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0]])
edges = np.array([[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [0, 6]])


corner = np.ones(3) / np.sqrt(3)
arc_mid1 = np.array([0.0, 1/np.sqrt(2), 1/np.sqrt(2)])
arc_mid2 = np.array([1/np.sqrt(2), 0.0, 1/np.sqrt(2)])
arc_mid3 = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0.0])

def great_circle(p1: np.ndarray, p2: np.ndarray, num_pts: int):
    n = np.cross(p1, p2)
    u = p1 / np.linalg.norm(p1)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    theta = np.arccos(np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2)))
    t = np.linspace(0, theta, num_pts)
    path = [(np.cos(t_i) * u + np.sin(t_i) * v) for t_i in t]
    
    return np.array(path)


arc_pts1 = great_circle(corner, arc_mid1, num_pts=20)
arc_edges1 = np.stack([[i, i+1] for i in range(len(arc_pts1)-1)], axis=0)

arc_pts2 = great_circle(corner, arc_mid2, num_pts=20)
arc_edges2 = np.stack([[i, i+1] for i in range(len(arc_pts2)-1)], axis=0)

arc_pts3 = great_circle(corner, arc_mid3, num_pts=20)
arc_edges3 = np.stack([[i, i+1] for i in range(len(arc_pts3)-1)], axis=0)


p.add_mesh(pv.PolyData(corner), color='red', point_size=10, smooth_shading=True, render_points_as_spheres=True)
p.add_mesh(pv.PolyData(np.array([1.0, 1.0, 1.0])), color='red', point_size=10, smooth_shading=True, render_points_as_spheres=True)
p.add_mesh(pv.Spline(arc_pts1, 100), color=(0.0, 0.0, 0.0), line_width=2.0, smooth_shading=True)
p.add_mesh(pv.Spline(arc_pts2, 100), color=(0.0, 0.0, 0.0), line_width=2.0, smooth_shading=True)
p.add_mesh(pv.Spline(arc_pts3, 100), color=(0.0, 0.0, 0.0), line_width=2.0, smooth_shading=True)


# 定义一个回调函数，在关闭窗口时保存截图
def save_screenshot():
    p.screenshot(png_file)
    print("Screenshot saved")

# 设置关闭窗口时的回调函数
p.add_key_event('s', save_screenshot)  # 按 's' 键保存截图

p.show()
