import os
import pyshtools as pysh
import numpy as np
import trimesh
import json
from scipy.special import sph_harm
import matplotlib.pyplot as plt
from matplotlib import rcParams

import polyscope as ps


config = {
    # 'figure.figsize': (8, 6),
    'lines.linewidth': 2.0, 
    # "font.family": 'Arial',
    'font.size': 24,
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    "axes.unicode_minus": False,
    'axes.linewidth': 2.0
}
rcParams.update(config)

out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'


mesh_names = ['spot', 'armadillo', '551074_sf', '275400_sf_simp']


for meshname in mesh_names:
    mesh_file = f'/media/box/Elements/Exp/HCPSOTS/data/{meshname}.obj'
    mesh = trimesh.load(mesh_file)

    # mesh_pts_file = f'/media/box/Elements/Exp/HCPSOTS/out/HCPSOTS_mesh_sampling_{meshname}.pts'
    mesh_pts_file = f'/media/box/Elements/Exp/HCPSOTS/out/NESOTS_mesh_sampling_{meshname}.pts'
    mesh_pts = np.loadtxt(mesh_pts_file)

    ps.init()
    ps.set_ground_plane_mode("shadow_only")  
    # ps.set_ground_plane_mode("tile")  
    ps.set_shadow_darkness(0.9)

    ps.register_surface_mesh(
        "mesh", 
        mesh.vertices, 
        mesh.faces,
        color=(0.89019608, 0.92941176, 0.92941176),
        enabled=True,
        edge_color=(0.0, 0.0, 0.0),
        edge_width=1.0,
        smooth_shade=True,
        transparency=1.0
        )

    ps.register_point_cloud(
        "points", 
        mesh_pts,
        color=(1.0, 0.0, 0.85490196),
        enabled=True,
        radius=0.0029,
        transparency=0.9
        )

    ps.show()

    # ps.screenshot(out_dir+f'show_HCPSOTS_mesh_sampling_{meshname}.jpg')
    ps.screenshot(out_dir+f'show_NESOTS_mesh_sampling_{meshname}.jpg')

