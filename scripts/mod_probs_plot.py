from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

unaligned = "/Users/joshualevendis/OfflineDocuments/RNA/modkit/modkit_sampleprobs_unaligned/28C1/probabilities.tsv"
pfal_mapq50 = "/Users/joshualevendis/OfflineDocuments/RNA/modkit/modkit_sampleprobs_pfal_mapq50/28C1/probabilities.tsv"

df = pd.read_csv(unaligned, sep="\t")
matches = df.loc[df['code'] == 'a']

counts = matches['count'].to_numpy()
probs = matches['range_start'].to_numpy()

print(len(probs))

label='unaligned, unfiltered'
x = np.arange(0, 128)
plt.bar(x, counts, color="darkgoldenrod", label=label)
plt.xlim(xmin=0, xmax=128)

xticks=np.arange(0, 129, 25.6).round(1)

print(xticks)
xticks_labels=np.arange(0.5, 1.05, 0.1).round(1)

print(xticks_labels)

plt.xticks(ticks=xticks, labels=xticks_labels)
plt.ylabel("count")
plt.xlabel("modification probability")
plt.legend(loc="upper right")

# plt.title("m6A probability distribution (P. falciparum MAPQ >= 50)")

plt.show()
