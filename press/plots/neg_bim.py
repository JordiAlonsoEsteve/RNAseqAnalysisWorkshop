import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom

def nb_pmf(k, mu, alpha):
    r = 1.0 / alpha
    p = 1.0 / (1.0 + alpha * mu)
    return nbinom.pmf(k, r, p)

alpha = 0.2
k_obs = 20
mus = {9: "tab:purple", 15: "tab:blue", 20: "tab:green",
       25: "tab:red", 31: "tab:orange"}

ref = nb_pmf(k_obs, 20, alpha)   # height of the mu = 20 curve at k = 20

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

# Left: PMFs in k, slice at k = k_obs
k = np.arange(0, 71)
ax = axes[0]
for m, c in mus.items():
    ax.plot(k, nb_pmf(k, m, alpha), "-o", ms=3, color=c, lw=1.4, label=rf"$\mu={m}$")
    ax.plot(k_obs, nb_pmf(k_obs, m, alpha), "o", color=c, ms=10,
            markeredgecolor="k", zorder=5)
ax.axvline(k_obs, color="k", ls="--", lw=1)
ax.axhline(ref, color="gray", lw=1.5, alpha=0.4)   # soft line
ax.set_xlabel("$k$")
ax.set_ylabel(r"$P(K=k\mid\mu,\alpha)$")
ax.set_title(rf"PMFs in $k$ ($\alpha={alpha}$); dots at $k={k_obs}$")
ax.legend()

# Right: likelihood in mu for fixed k_obs
mu_grid = np.linspace(2, 60, 1000)
lik = nb_pmf(k_obs, mu_grid, alpha)
ax = axes[1]
ax.plot(mu_grid, lik, "k-", lw=1.5)
for m, c in mus.items():
    ax.plot(m, nb_pmf(k_obs, m, alpha), "o", color=c, ms=10,
            markeredgecolor="k", zorder=5, label=rf"$\mu={m}$")
ax.axvline(k_obs, color="k", ls="--", lw=1)
ax.axhline(ref, color="gray", lw=1.5, alpha=0.4)   # soft line
ax.set_xlabel(r"$\mu$")
ax.set_ylabel(rf"$L(\mu\mid k={k_obs},\alpha)$")
ax.set_title(rf"Likelihood in $\mu$: maximum at $\hat\mu={mu_grid[lik.argmax()]:.1f}$")
ax.legend()

plt.tight_layout()
plt.show()