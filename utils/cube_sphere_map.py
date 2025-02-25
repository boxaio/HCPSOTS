import os
import numpy as np
import trimesh
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams

import pyvista as pv

from hilbert_curve import hilbert_curve, curve_func

out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'


hc_pts = np.array(hilbert_curve(0, 0, 1, 0, 0, 1, n=3))
c_pts = curve_func(hc_pts)


def pq_to_xy(pq:np.ndarray):
    # note that p > q
    p, q = pq[:,0], pq[:,1]
    S = np.sin(np.pi*p*p/6)
    T = q/p
    # note that x > y
    x = np.sqrt(S/(1-S))
    y = T*x*np.sqrt((1+x*x)/(1+2*x*x-T*T*x*x))

    return x, y

    # return np.stack((x, y)).transpose()

def xy_to_XYZ(x:np.ndarray, y:np.ndarray):
    Z = 1/np.sqrt(1+x*x+y*y)
    X = x*Z
    Y = y*Z

    return np.stack((X, Y, Z)).transpose()


# print(np.linalg.norm(XYZ, axis=-1))

def grid_points_in_right_triangle(step):
    x = np.arange(0.0, 1.0+step, step)
    y = np.arange(0.0, 1.0+step, step)
    
    grid_points = []
    for xi in x:
        for yi in y:
            if xi >= yi:
                if not (xi == 0 and yi == 0):
                    grid_points.append((xi, yi))
                    
    return np.array(grid_points)


def cube_face_map_to_sphere(pq_grid_pts: np.ndarray):
    assert pq_grid_pts.shape[1] == 2
    # (p, q) to (x, y)
    mask1 = pq_grid_pts[:,0] >= pq_grid_pts[:,1]
    mask2 = pq_grid_pts[:,0] < pq_grid_pts[:,1]

    x1, y1 = pq_to_xy(pq_grid_pts[mask1])
    y2, x2 = pq_to_xy(np.fliplr(pq_grid_pts[mask2]))
    
    # (x, y) to (X, Y, Z)
    sphere_pts = np.zeros((len(pq_grid_pts), 3))
    sphere_pts[mask1] = xy_to_XYZ(x1, y1)
    sphere_pts[mask2] = xy_to_XYZ(x2, y2)

    return sphere_pts


def xy_to_pq(xy: np.ndarray):
    # note that x > y
    assert xy.shape[1] == 2
    x, y = xy[:,0], xy[:,1]
    p = np.sqrt(6*np.arcsin(x*x/(1+x*x))/np.pi)
    q = p*y*np.sqrt((1+2*x*x)/(1+x*x+y*y))/x

    return p, q

def sphere_map_to_cube_face(sphere_pts: np.ndarray):
    assert sphere_pts.shape[1] == 3
    # (X, Y, Z) to (x, y)
    x = sphere_pts[:,0] / sphere_pts[:,2]
    y = sphere_pts[:,1] / sphere_pts[:,2]
    xy = np.stack((x, y)).transpose()

    mask1 = x >= y
    mask2 = x < y

    # (x, y) tp (p, q)
    p1, q1 = xy_to_pq(xy[mask1])
    q2, p2 = xy_to_pq(np.fliplr(xy[mask2]))

    pq_pts = np.stack((np.hstack((p1, p2)), np.hstack((q1, q2)))).transpose()

    return pq_pts


def great_circle(p1: np.ndarray, p2: np.ndarray, num_pts: int=1, return_distance: bool=False):
    assert np.abs(np.linalg.norm(p1)-1.0)<1e-6 and np.abs(np.linalg.norm(p2)-1.0)<1e-6
    n = np.cross(p1, p2)
    u = p1 / np.linalg.norm(p1)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    theta = np.arccos(np.dot(p1, p2))
    # print(np.dot(p1, p2), p1, p2)
    t = np.linspace(0, theta, num_pts)
    path = [(np.cos(t_i) * u + np.sin(t_i) * v) for t_i in t]
    
    if return_distance:
        # distance = np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))
        distance = theta
        return distance, np.array(path)
    else:
        return np.array(path)



if __name__ == '__main__':
    
    hc_pts = np.array(hilbert_curve(0, 0, 1, 0, 0, 1, n=4))
    hc_pts = curve_func(hc_pts)
    hc_length = np.linalg.norm(np.diff(hc_pts, axis=0), axis=-1).sum()

    print("Hilbert curve length", hc_length)

    # step = 0.1
    # pq_grid_pts = [(p, q) for p in np.arange(0.0, 1.0+step, step) for q in np.arange(0.0, 1.0+step, step)]
    # pq_grid_pts = np.array(pq_grid_pts)

    pq_grid_pts = hc_pts

    sphere_pts = cube_face_map_to_sphere(pq_grid_pts)
    # shc_length = np.sum(np.array([great_circle(p1, p2, return_distance=True)[0] for p1, p2 in zip(sphere_pts[0:-1], sphere_pts[1:])]))
    shc_length = np.linalg.norm(np.diff(sphere_pts, axis=0), axis=-1).sum()
    print("Spherical Hilbert curve length", shc_length)

    pq_grid_pts = np.hstack((pq_grid_pts, np.ones([len(pq_grid_pts), 1])))

    # # random points on the sphere octant
    # lons = np.linspace(0.05, np.pi/4, 9)
    # lats = np.linspace(np.pi/4, np.pi/2, 9)
    # sphere_pts = [[np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)] for lon in lons for lat in lats]
    # sphere_pts = np.array(sphere_pts)

    # pq_grid_pts = sphere_map_to_cube_face(sphere_pts)
    # pq_grid_pts = np.hstack((pq_grid_pts, np.ones([len(pq_grid_pts), 1])))


    sphere = pv.Sphere(radius=1.0, center=(0, 0, 0), direction=(0, 0, 1), theta_resolution=179, phi_resolution=80)
    octant = sphere.clip('x',origin=(0,0,0),invert=False).clip('y',origin=(0,0,0),invert=False).clip('z',origin=(0,0,0),invert=False)

    png_file = out_dir + f'HC_cube_sphere_octant.png'

    plane1 = pv.Plane(center=(0.5, 0.5, 1.0), direction=(0, 0, 1), i_size=1, j_size=1)
    plane2 = pv.Plane(center=(1.0, 0.5, 0.5), direction=(1, 0, 0), i_size=1, j_size=1)
    plane3 = pv.Plane(center=(0.5, 1.0, 0.5), direction=(0, 1, 0), i_size=1, j_size=1)

    p = pv.Plotter(window_size=(1340, 1100))
    p.add_mesh(octant, color=(1.0, 0.376, 0.91), smooth_shading=True, opacity=0.5)
    p.add_mesh(plane1, color=(0.89019608, 0.92941176, 0.92941176), opacity=0.3, smooth_shading=True)
    p.add_mesh(plane2, color=(0.89019608, 0.92941176, 0.92941176), opacity=0.3, smooth_shading=True)
    p.add_mesh(plane3, color=(0.89019608, 0.92941176, 0.92941176), opacity=0.3, smooth_shading=True)
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


    corner = np.ones(3)/np.sqrt(3)
    arc_mid1 = np.array([0.0, 1/np.sqrt(2), 1/np.sqrt(2)])
    arc_mid2 = np.array([1/np.sqrt(2), 0.0, 1/np.sqrt(2)])
    arc_mid3 = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0.0])



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

    polydata_plane = pv.PolyData()
    polydata_plane.points = pq_grid_pts
    polydata_plane.lines = np.hstack(([len(pq_grid_pts),] + list(range(len(pq_grid_pts)))))

    polydata_sphere = pv.PolyData()
    polydata_sphere.points = sphere_pts
    polydata_sphere.lines = np.hstack(([len(sphere_pts),] + list(range(len(sphere_pts)))))

    p.add_mesh(polydata_plane, color='blue', line_width=2, smooth_shading=True, render_points_as_spheres=True)
    p.add_mesh(polydata_sphere, color='red', line_width=2, smooth_shading=True, render_points_as_spheres=True)

    # 定义一个回调函数，在关闭窗口时保存截图
    def save_screenshot():
        p.screenshot(png_file)
        print("Screenshot saved")

    # 设置关闭窗口时的回调函数
    p.add_key_event('s', save_screenshot)  # 按 's' 键保存截图

    p.show()

