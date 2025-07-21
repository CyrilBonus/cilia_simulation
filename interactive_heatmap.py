import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap

yellow_red = LinearSegmentedColormap.from_list('yellow_red', ['yellow', 'red'])

# Load .dat files
dat_files = sorted(
    [f for f in glob.glob("*.dat") if re.match(r"^\d+\.dat$", f)],
    key=lambda x: int(x.split('.')[0])
)
records = []
for idx, f in enumerate(dat_files):
    with open(f, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('//') or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) == 5:
                x = int(parts[1])
                y = int(parts[2])
                phi = float(parts[3])
                records.append({'frame': idx, 'x': x, 'y': y, 'phi': phi})

df = pd.DataFrame(records)
frames = sorted(df['frame'].unique())

# Interpolation grid
X_vals = sorted(df['x'].unique())
Y_vals = sorted(df['y'].unique())
grid_x, grid_y = np.meshgrid(
    np.linspace(min(X_vals), max(X_vals), 100),
    np.linspace(min(Y_vals), max(Y_vals), 100)
)

# Fixed phi scale
phi_min = 0
phi_max = 6

# Plot setup for 3D
fig3d = plt.figure(figsize=(10, 7))
ax3d = fig3d.add_subplot(projection='3d')
plt.subplots_adjust(bottom=0.18)

# Plot setup for 2D heatmap
fig2d, ax2d = plt.subplots(figsize=(7, 7))
plt.subplots_adjust(bottom=0.18)

# Colorbar for 3D
mappable = plt.cm.ScalarMappable(cmap='plasma')
mappable.set_array([phi_min, phi_max])
mappable.set_clim(phi_min, phi_max)
plt.colorbar(mappable, ax=ax3d, shrink=0.6, pad=0.1, label='φ')

# Slider
ax_slider = fig3d.add_axes([0.2, 0.05, 0.6, 0.04])
slider = Slider(ax=ax_slider, label='Frame', valmin=0, valmax=len(frames) - 1, valinit=0, valstep=1, color='orange')

im = None
colorbar2d = [None]

def get_phi_grid(frame_df, X_vals, Y_vals):
    phi_grid = np.full((len(Y_vals), len(X_vals)), np.nan)
    x_to_idx = {x: i for i, x in enumerate(X_vals)}
    y_to_idx = {y: i for i, y in enumerate(Y_vals)}
    for _, row in frame_df.iterrows():
        xi = x_to_idx[row['x']]
        yi = y_to_idx[row['y']]
        phi_grid[yi, xi] = row['phi']
    return phi_grid

def update(val):
    ax3d.cla()
    idx = int(slider.val)
    frame_df = df[df['frame'] == frames[idx]]

    # 3D Surface plot
    xs = frame_df['x'].values
    ys = frame_df['y'].values
    zs = frame_df['phi'].values
    grid_z = griddata((xs, ys), zs, (grid_x, grid_y), method='linear')
    grid_z = np.clip(grid_z, phi_min, phi_max)
    ax3d.plot_surface(grid_x, grid_y, grid_z, cmap='plasma', vmin=phi_min, vmax=phi_max, linewidth=0, antialiased=True)
    ax3d.set_title(f'3D Surface - Frame = {idx} ({dat_files[idx]})')
    ax3d.set_xlabel('x')
    ax3d.set_ylabel('y')
    ax3d.set_zlabel('φ')
    ax3d.set_xlim(min(X_vals), max(X_vals))
    ax3d.set_ylim(min(Y_vals), max(Y_vals))
    ax3d.set_zlim(phi_min, phi_max)

    # 2D Heatmap
    phi_grid = get_phi_grid(frame_df, X_vals, Y_vals)
    global im
    if im is None:
        im = ax2d.imshow(
            np.ma.masked_where(np.isnan(phi_grid) | (phi_grid <= 0), phi_grid),
            cmap=yellow_red,
            origin='lower',
            extent=(min(X_vals), max(X_vals), min(Y_vals), max(Y_vals)),
            aspect='equal',
            vmin=0, vmax=phi_max
        )
        ax2d.set_title(f'2D Heatmap - Frame = {idx} ({dat_files[idx]})')
        ax2d.set_xlabel('x')
        ax2d.set_ylabel('y')
        
        colorbar2d[0] = fig2d.colorbar(im, ax=ax2d, label='φ', shrink=0.8)
    else:
        im.set_data(np.ma.masked_where(np.isnan(phi_grid) | (phi_grid <= 0), phi_grid))
        ax2d.set_title(f'2D Heatmap - Frame = {idx} ({dat_files[idx]})')
    fig2d.canvas.draw_idle()
    fig3d.canvas.draw_idle()

slider.on_changed(update)
update(0)
plt.show()
