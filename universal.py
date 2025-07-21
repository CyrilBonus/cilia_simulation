import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
from matplotlib.widgets import Button
from tkinter import Tk
from tkinter.filedialog import askopenfilename, asksaveasfilename

Tk().withdraw()

file_path = askopenfilename(title="Select a data file")

if not file_path:
    print("No file selected.")
    exit()

x_vals = []
y_vals = []
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip().strip('"')
        if not line:
            continue
        try:
            x_str, y_str = line.split(',')
            x_vals.append(float(x_str))
            y_vals.append(float(y_str))
        except ValueError:
            print(f"Skipping invalid line: {line}")

selected_points = []

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

for ax in (ax1, ax2):
    ax.scatter(x_vals, y_vals, color='blue')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True)

ax1.set_title('Left Plot')
ax2.set_title('Right Plot with Selector')

def onselect(eclick, erelease):
    global selected_points
    xmin, xmax = sorted([eclick.xdata, erelease.xdata])
    ymin, ymax = sorted([eclick.ydata, erelease.ydata])
    selected_points = [(x, y) for x, y in zip(x_vals, y_vals)
                       if xmin <= x <= xmax and ymin <= y <= ymax]
    print(f"Selected {len(selected_points)} points.")

def save_selected(event):
    if not selected_points:
        print("No points selected.")
        return
    save_path = asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        title="Save selected points as"
    )
    if save_path:
        df_selected = pd.DataFrame(selected_points, columns=['x', 'y'])
        df_selected.to_csv(save_path, index=False, header=False)
        print(f"Saved {len(selected_points)} points to '{save_path}'.")

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

save_ax = plt.axes([0.45, 0.01, 0.1, 0.05])
save_button = Button(save_ax, 'Save', color='white', hovercolor='lightblue')
save_button.on_clicked(save_selected)

fig.set_constrained_layout(True)
plt.show()
