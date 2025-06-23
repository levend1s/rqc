import pandas as pd
from matplotlib import pyplot as plt
import numpy as np


groups = [
    [1.196645609,	1.399813891,	0.971040555,	1.008199288],
    [1.243606457,	1.071557157,	0.891821486,	0.986118786],
    [1.187868377,	1.129854038,	1.007056742,	1.079522563]
]

# groups = [
#     [4.434065345,	5.363682623,	4.36245526,	4.506259078],
#     [4.79268054,	4.076326369,	4.08474427,	4.405089686],
#     [4.52760482,	4.239400869,	4.451438246,	4.677009644]
# ]
group_labels = ['28 hpi', '32 hpi', '36 hpi']

# Convert data to pandas DataFrame.


df = pd.DataFrame(np.array(groups), index=group_labels).T
# df.columns = ["C1", "C2", "K1", "K2"]

c = {
    "C1": "steelblue",
    "C2": "cornflowerblue",
    "K1": "firebrick",
    "K2": "indianred"
}

# Plot.
pd.concat(
    [
        df.iloc[0].rename('C1'), 
        df.iloc[1].rename('C2'), 
        df.iloc[2].rename('K1'),
        df.iloc[3].rename('K2')
    ],
    axis=1,
).plot.bar(color=['steelblue', 'cornflowerblue', 'firebrick', 'indianred'])

plt.xticks(rotation=0)
plt.ylabel("relative m6A/A abundance (%)")
# plt.title("relative m6A/A abundance (m6A probability >= 0.9)")
plt.legend(loc="upper right", ncol=4)

plt.savefig('filename.png', dpi=300)
plt.show()