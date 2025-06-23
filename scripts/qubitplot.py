from matplotlib import pyplot as plt
import numpy as np

# Generate 100 random data points along 3 dimensions
#x, y, scale = np.random.randn(3, 100)
groups=("28 hpi", "32 hpi", "36 hpi")
qubit_readings = {
    "C1": [1.42, 10.99, 2.33],
    "C2": [3.67, 10.03, 16.56],
    "K1": [5.71, 15.84, 15.6],
    "K2": [1.46, 12.48, 4.58]
    
}

colors = {
    "C1": "steelblue",
    "C2": "cornflowerblue",
    "K1": "firebrick",
    "K2": "indianred"
}

x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in qubit_readings.items():
    offset = width * multiplier

    this_color = colors[attribute]
    print(this_color)
    
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=this_color)
    # ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('RNA yield (µg)')
# ax.set_title('RNA yield')
ax.set_xticks(x + 0.3, groups)
ax.set_yticks(range(0, 20, 2))
ax.legend(loc='upper left')
# ax.set_ylim(0, 20)
print(x + width)
print(groups)
plt.axhline(1, linewidth=1, linestyle="--", color='black')

plt.show()