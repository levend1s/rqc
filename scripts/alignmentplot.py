from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# Generate 100 random data points along 3 dimensions
#x, y, scale = np.random.randn(3, 100)
groups=("28 hpi", "32 hpi", "36 hpi")

colors = {
    "P. falciparum": "steelblue",
    "S. cerevisae": "orange",
    "H. sapien": "green",
    "Unaligned": "gray"
}

d_unfiltered = {
    "P. falciparum": [62.59,	56.30,	58.53,	72.04,	65.45,	50.61,	68.47,	67.08,	79.52,	67.09,	64.70,	70.47],
    "S. cerevisiae": [13.42,	29.24,	18.37,	18.46,	17.84,	29.26,	14.54,	20.52,	12.96,	14.25,	8.51,	6.64],
    "H. sapien": [0.99,	2.51,	1.50,	1.15,	1.33,	0.96,	1.40,	1.22,	0.82,	0.83,	0.98,	0.59],
    "Unaligned": [22.99,	11.95,	21.60,	8.35,	15.37,	19.16,	15.59,	11.18,	6.70,	17.83,	25.81,	22.29]
}

d_filtered = {
    "P. falciparum": [37.71,  24.33 , 35.96 , 38.39, 30.00 , 25.47 , 37.90 , 32.09, 45.51  ,35.56 , 35.13 , 28.87],
    "S. cerevisiae": [5.36 , 19.38 , 10.96  , 7.41, 8.93 , 21.81 ,  5.03 , 11.26, 1.99  , 6.76  , 4.51  , 1.70],
    "H. sapien": [0.07 ,  0.09  , 0.06 ,  0.04, 0.03  , 0.03   ,0.04  , 0.0, 0.07  , 0.08  , 0.07,   0.06],
    "MAPQ < 50": [56.85 , 56.21 , 53.02 , 54.16, 61.03 , 52.69,  57.04  ,56.60, 52.43 , 57.60 , 60.29 , 69.36]
}

orgs = ["P. falciparum", "S. cerevisiae", "H. sapien", "Unaligned"]

sample = ["C1", "C2", "K1", "K2", "C1", "C2", "K1", "K2", "C1", "C2", "K1", "K2"]
timepoint = ["28 hpi", "28 hpi", "28 hpi", "28 hpi", "32 hpi", "32 hpi", "32 hpi", "32 hpi", "36 hpi", "36 hpi", "36 hpi", "36 hpi"]
sample_name = ["28C1", "28C2", "28K1", "28K2", "32C1", "32C2", "32K1", "32K2", "36C1", "36C2", "36K1", "36K2"]

d = d_filtered
df = pd.DataFrame(d)
df['time point'] = timepoint
df['sample'] = sample
df['sample_name'] = sample_name



print(df)

df.plot.bar(x="sample_name", stacked=True)
# df.pivot_table(index='sample', columns='time point').plot(kind = 'bar', stacked = True)

plt.ylabel("% aligned")
plt.xticks(rotation=45)
plt.xlabel("")
plt.legend(loc="upper right")
# plt.title("Sample alignment to reference genomes")
plt.show()