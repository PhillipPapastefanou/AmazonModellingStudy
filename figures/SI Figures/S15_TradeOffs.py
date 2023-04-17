import  pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

data = pd.read_csv('../woody_trait.0827.txt', sep='\t')
dataIso = pd.read_csv('../Parameters_v1.7.3.csv', sep=',')

import string

lowerletters = string.ascii_lowercase[:26]

fig = plt.figure(figsize=(4, 9))

psi50 = data['P50']
id = 0

for name in ['LMA', 'Ks']:
    slice = data[name]

    nanpos = np.isnan(psi50)
    psi503 = psi50[~nanpos]
    slice = slice[~nanpos]

    nanpos = np.isnan(slice)
    psi503 = psi503[~nanpos]
    slice = slice[~nanpos]

    slope, intercept, r_value, p_value, std_err = stats.linregress(psi503, slice)
    psi50l = np.arange(np.min(psi503), np.max(psi503), 0.1)

    ax = fig.add_subplot(3, 1, id +1)
    ax.scatter(psi503, slice, s = 2)
    y = slope * psi50l + intercept
    ax.plot(psi50l, y, color = 'black')
    ax.set_ylabel("$\log{" + name + "}$")
    ax.set_xlabel("$\log{ ( -\psi_{50}) }$")
    ax.text(0.05, 0.18, r"$\log{" + name + "} = " + str(np.round(slope, 2)) + "\cdot \log{(-\psi_{50})} +" + str(np.round(intercept,2)) + "$", horizontalalignment='left',
            verticalalignment='center', c="black",
            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))
    ax.text(0.05, 0.08, r"$p  < 0.001$", horizontalalignment='left',
            verticalalignment='center', c="black",
            transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))
    ax.text(0.03, 0.92, lowerletters[id] + ")", horizontalalignment='left',
               verticalalignment='center',
               transform=ax.transAxes, size=10,
               bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))
    id += 1

p = dataIso['psi50s'].values
iso = dataIso['iso'].values
ax = fig.add_subplot(3, 1, id +1)
ax.plot(p, iso)
ax.set_ylabel("isohydricity $\lambda$")
ax.set_xlabel("$\psi_{50}$")
ax.text(0.05, 0.08, r"$\lambda = 0.65 + 0.15 \cdot \psi_{50}$", horizontalalignment='left',
        verticalalignment='center', c="black",
        transform=ax.transAxes, size=8, bbox=dict(facecolor='white', alpha=1.0, edgecolor='white'))
ax.text(0.03, 0.92, lowerletters[2] + ")", horizontalalignment='left',
        verticalalignment='center',
        transform=ax.transAxes, size=10,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='white'))

plt.subplots_adjust(wspace=0.4, left= 0.2)
plt.savefig('Figure_Empirical_TradeOffs.png', dpi = 300, bbox_inches='tight',
                    pad_inches=0)




