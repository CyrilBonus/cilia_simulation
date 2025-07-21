import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector

# Load file with no header and treat each line as a single string
with open('cluster_centers.csv', 'r') as f:
    lines = f.readlines()

x_vals = []
y_vals = []
for line in lines:
    line = line.strip().strip('"')
    x_str, y_str = line.split(',')
    x_vals.append(float(x_str))
    y_vals.append(float(y_str))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

ax1.scatter(x_vals, y_vals, color='blue')
ax1.set_title('Left Plot')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)

ax2.scatter(x_vals, y_vals, color='blue')
ax2.set_title('Right Plot with Selector')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.grid(True)

def onselect(eclick, erelease):
    xmin, xmax = sorted([eclick.xdata, erelease.xdata])
    ymin, ymax = sorted([eclick.ydata, erelease.ydata])

    selected_points = [(x, y) for x, y in zip(x_vals, y_vals)
                       if xmin <= x <= xmax and ymin <= y <= ymax]

    print(f"Selected {len(selected_points)} points.")
    if selected_points:
        df_selected = pd.DataFrame(selected_points, columns=['x', 'y'])
        df_selected.to_csv('selected_points.csv', index=False, header=False)
        print("Saved selected points to 'selected_points.csv'.")

selector = RectangleSelector(
    ax2,
    onselect,
    useblit=False,
    button=[1],
    minspanx=0, minspany=0,
    spancoords='data',
    interactive=True,
    props=dict(edgecolor='black', facecolor='none', linewidth=1),
    handle_props=dict(marker='s', markersize=5, markerfacecolor='black', markeredgecolor='black')
)

plt.tight_layout()
plt.show()
