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


base_path = '/media/box/Elements/Exp/sphere-code/results/'
tmp = '2024-12-31/datafiles'
N = 3000
out_dir = '/media/box/Elements/Exp/HCPSOTS/out/'

methods = ['regular', 'whitenoise', 'dartthrowing', 'stratified', 'poissondisk', 'NESOTS', 'HCPSOTS']
labels = ['regular', 'white noise', 'Dart Throwing', 'Stratified', 'Poisson Disk', 'NESOTS', 'HCPSOTS']
colors = ['lightsteelblue', 'slategray', 'purple', 'gold', 'forestgreen', 'cyan', 'firebrick']
alpha = [0.5, 0.9, 0.5, 0.5, 0.5, 0.5, 1.0]

spectras = []
samples = []
for m in methods:
    if m == 'NESOTS':    
        spect_file = base_path + f"powspec-sphere-{m}-n{N}.txt"
        pts_file = f"/media/box/Elements/Exp/NESOTS_debug/out/NESOTS_spherical_bluenoise.pts"
    elif m == 'HCPSOTS':
        spect_file = base_path + f"powspec-sphere-{m}-n{N}.txt"
        pts_file = f"/media/box/Elements/Exp/HCPSOTS/out/HCPSOTS_spherical_bluenoise.pts"
    else:
        spect_file = base_path + f"powspec-sphere-{m}-n{N}/{tmp}/powspec-sphere-{m}-n{N}.txt"
        pts_file = base_path + f"powspec-sphere-{m}-n{N}/{tmp}/samples-{m}-n{N}.txt"
    spect = np.loadtxt(spect_file)
    pts = np.loadtxt(pts_file)
    spectras.append(spect)
    samples.append(pts)


fig, ax = plt.subplots(figsize=(14, 9))
for i in range(len(spectras)):
    plt.plot(spectras[i][1:,0], spectras[i][1:,1], linestyle='-', color=colors[i], alpha=alpha[i],label=labels[i])
plt.xlabel('Normalized Frequencies')
plt.ylabel('Angular Power Spectrum')
ax.set_xlim(left=0)
ax.set_ylim(bottom=0, top=3.0)
plt.grid(True, which='both', ls='-.')
plt.legend(fontsize=18, loc='upper left', frameon=False)
plt.savefig(out_dir+f"powspec-sphere-all-n{N}.png", dpi=600)


sphere_cvt = trimesh.load('/media/box/Elements/Exp/HCPSOTS/sphere_cvt.obj')

for i in range(len(samples)):
    ps.init()
    ps.set_ground_plane_mode("shadow_only")  
    # ps.set_ground_plane_mode("tile")  
    ps.set_shadow_darkness(0.9)

    ps.set_screenshot_extension(".png")

    ps.register_surface_mesh(
        "sphere", 
        sphere_cvt.vertices, 
        sphere_cvt.faces,
        color=(0.89019608, 0.92941176, 0.92941176),
        # enabled=True,
        # edge_color=(0.0, 0.0, 0.0),
        # edge_width=1.0,
        # smooth_shade=True,
        transparency=1.0
        )

    ps.register_point_cloud(
        "points {:}".format(methods[i]), 
        samples[i],
        color=(1.0, 0.0, 0.85490196),
        enabled=True,
        radius=0.0035,
        transparency=0.9
        )
    # ps.reset_camera_to_home_view()

    # Look for camera parameters JSON and use it, if it is available and not to be ignored;
    # Otherwise, open GUI so that the user could choose the best camera view
    camera_json = f"/media/box/Elements/Exp/HCPSOTS/utils/view_sphere_{methods[i]}.json"
    if os.path.exists(camera_json):
        os.remove(camera_json)

    # View the point cloud and mesh we just registered in the 3D UI
    ps.show()

    camera = json.loads(ps.get_view_as_json())

    # Save camera parameters string into JSON
    with open(camera_json, 'w') as json_file:
        json.dump(camera, json_file, sort_keys=True, indent=4)

    # Load the camera parameters
    with open(camera_json) as json_file:
        camera = json_file.read()

    ps.set_view_from_json(camera)

    ps.screenshot(out_dir+f'show_{methods[i]}.jpg')
    ps.remove_all_structures()

