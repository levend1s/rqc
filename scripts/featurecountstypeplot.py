import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


groups = [
    [3496956,	3206,	593946,	14377,	1439,	3303],
    [3283650,	2339,	823046,	23185,	1396,	2534],
    [3393375,	2095,	1628387,	27022,	1175,	4349]
]

groups = [
    [3496956,	3283650,	3393375],
    [593946,	823046	,1628387],
    [1439,	1396	,1175],
    [14377,	23185,	27022],
    [3303,	2534	,4349],
    [3206,	2339,	2095]
]


# group_labels = ['28 hpi', '32 hpi', '36 hpi']
group_labels = ['mRNA', 'rRNA', 'tRNA', 'snoRNA', 'snRNA', 'pseudogene']


# Convert data to pandas DataFrame.
df = pd.DataFrame(groups, index=group_labels).T
# df.columns = ["C1", "C2", "K1", "K2"]

# Plot.
pd.concat(
    [
        df.iloc[0].rename('28 hpi'), 
        df.iloc[1].rename('32 hpi'), 
        df.iloc[2].rename('36 hpi'),
    ],
    axis=1,
).plot.bar(color=['crimson', 'darkslategray', 'darkgoldenrod'])

plt.xticks(rotation=0)
plt.ylabel("transcript count")
# plt.title("total mapped transcript types by timepoint")
# plt.legend(loc="upper right", ncol=4)
ax = plt.gca()
ax.set_yscale('log')
# ax.set_ylim(0)

groups = [
    [89730	,66807,	172314],
    [116	,0,	0],
    [0	,1122,	3260]
]

group_labels = ['mRNA', 'pseudogene', 'snoRNA']


# Convert data to pandas DataFrame.
df = pd.DataFrame(groups, index=group_labels).T
# df.columns = ["C1", "C2", "K1", "K2"]

# Plot.
pd.concat(
    [
        df.iloc[0].rename('28 hpi'), 
        df.iloc[1].rename('32 hpi'), 
        df.iloc[2].rename('36 hpi'),
    ],
    axis=1,
).plot.bar(color=['crimson', 'darkslategray', 'darkgoldenrod'])

plt.xticks(rotation=0)
plt.ylabel("transcript count")
# plt.title("differentially methylated transcript types by timepoint")
# plt.legend(loc="upper right", ncol=4)
ax = plt.gca()
ax.set_yscale('log')
# ax.set_ylim(0)

fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

data = [
    10173981,
    3045379,
    64584,
    10186,
    7640,
    4010,
]
ingredients = [
    "mRNA",
    "rRNA",
    "snoRNA",
    "snRNA",
    "pseudogene",
    "tRNA"
]


def func(pct, allvals):
    absolute = int(np.round(pct/100.*np.sum(allvals)))
    return f"{pct:.1f}%\n({absolute:d})"


wedges, texts = ax.pie(data)

percent = 100. * np.array(data) / sum(data)

print(sum(data))

labels = ['{0} - {1:1.2f}%'.format(i,j) for i,j in zip(ingredients, percent)]

print(labels)

ax.legend(wedges, labels,
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1),
          fontsize=20)

# plt.setp(autotexts, size=8, weight="bold")

# ax.set_title("Total feature counts, all samples")


plt.show()