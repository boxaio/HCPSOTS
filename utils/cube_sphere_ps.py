import os
import pyshtools as pysh
import numpy as np
import trimesh
import json
from scipy.special import sph_harm
import matplotlib.pyplot as plt
from matplotlib import rcParams

import polyscope as ps

out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'

# unit sphere
sphere = trimesh.primitives.Sphere(radius=1.0, subdivisions=5)

# define the surrounding cube
cube = trimesh.primitives.Box(extents=[2.0, 2.0, 2.0])

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

ps.init()
ps.set_ground_plane_mode("shadow_only")  
ps.set_shadow_darkness(0.)


ps.register_surface_mesh("sphere", sphere.vertices, sphere.faces, color=(1.0, 0.376, 0.91), enabled=True,
    # edge_color=(0.0, 0.0, 0.0),
    # edge_width=1.0,
    smooth_shade=True, transparency=0.26
    )
ps.register_surface_mesh("cube", cube.vertices, cube.faces, color=(0.275966, 0.792482, 0.8), enabled=True,
    # edge_color=(0.0, 0.0, 0.0),
    # edge_width=2.0,
    smooth_shade=True, transparency=0.52
    )
ps.set_transparency_mode('simple')
ps.register_curve_network("octant", vtx, edges, color=(0.0, 0.0, 0.0), radius=0.0005, enabled=True, material='flat')
ps.register_curve_network("arc1", arc_pts1, arc_edges1, color=(0.0, 0.0, 0.0), radius=0.0005, enabled=True, material='flat')
ps.register_curve_network("arc2", arc_pts2, arc_edges2, color=(0.0, 0.0, 0.0), radius=0.0005, enabled=True, material='flat')
ps.register_curve_network("arc3", arc_pts3, arc_edges3, color=(0.0, 0.0, 0.0), radius=0.0005, enabled=True, material='flat')
ps.register_point_cloud("corner", corner[None], color=(1.0, 0.0, 0.0), radius=0.003, point_render_mode='sphere')

ps_plane0 = ps.add_scene_slice_plane()
# ps_plane0.set_draw_plane(True) # render the semi-transparent gridded plane
ps_plane0.set_pose(plane_position=np.array([0,0,0]), plane_normal=(1.0, 0.0, 0.0))
# ps_mesh0.set_ignore_slice_plane(ps_plane0, True)

ps_plane1 = ps.add_scene_slice_plane()
ps_plane1.set_pose(plane_position=np.array([0,0,0]), plane_normal=(0.0, 1.0, 0.0))
# ps_mesh1.set_ignore_slice_plane(ps_plane1, True)

ps_plane2 = ps.add_scene_slice_plane()
ps_plane2.set_pose(plane_position=np.array([0,0,0]), plane_normal=(0.0, 0.0, 1.0))

camera_json = out_dir + f"octant_view_mat.json"
if os.path.exists(camera_json):
    os.remove(camera_json)

ps.show()

camera = json.loads(ps.get_view_as_json())
# Save camera parameters string into JSON
with open(camera_json, 'w') as json_file:
    json.dump(camera, json_file, sort_keys=True, indent=4)

# Load the camera parameters
with open(camera_json) as json_file:
    camera = json_file.read()

ps.set_view_from_json(camera)

ps.screenshot(out_dir+f'cube_sphere_octant.jpg')


