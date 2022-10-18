import numpy as np
import matplotlib.pyplot as plt
from correlate import CCG

# Generate spikes
n_unit = 3
duration = 3600
max_fr = 20
fr = np.random.rand(n_unit) * max_fr
iti = -np.log(np.random.rand(n_unit, max_fr * duration)) / fr[:, np.newaxis]
iti += np.random.rand(n_unit, max_fr * duration) * 0.002 + 0.001
sps = np.cumsum(iti, axis=1)

# Inject correlated spikes from unit 0 to unit 1
sp = sps[0].copy()
random_sp = sp[np.random.randint(low=0, high=len(sp), size=len(sp) // 2)] + (
    0.002 * np.random.randn(len(sp) // 2)) + 0.005

spikes = np.concatenate(sps)
indices = np.concatenate([i * np.ones_like(s) for i, s in enumerate(sps)])
spikes = np.concatenate([spikes, random_sp])
indices = np.concatenate([indices, np.ones_like(random_sp)])

sort_idx = np.argsort(spikes)
spikes = spikes[sort_idx]
indices = indices[sort_idx]

indices = indices[spikes <= duration]
spikes = spikes[spikes <= duration]

# Run CCG
fs = 3e4
ccg = CCG((fs * spikes).astype(int), indices.astype(int), fs=fs)

# Plot
plt.clf()
plt.ion()
for i in range(n_unit):
    for j in range(n_unit):
        k = n_unit * i + j + 1
        plt.subplot(n_unit, n_unit, k)
        plt.bar(np.arange(-25, 25), ccg[i][j])
        plt.title(f'{i} -> {j}')
plt.show()
