#!/usr/bin/env python3
"""Figures for "Short-interval Moebius sums, Chowla correlations, and polynomial
Moebius sums in characteristic 2".  Requires numpy and matplotlib.

    python make_figures.py        # writes ../figures/fig_planes.pdf,
                                  #        ../figures/fig_chowla.pdf,
                                  #        ../figures/fig_polysums.pdf

Everything is computed from tau(H) = sum_{j odd} h_j^2 p_{2j-1}(H) and the parity
law mu(H^2 + a t) = (-1)^{k eps_d} chi(tau(H)/a); the parity law itself is checked
by direct factorisation in verify_mobius_chowla.py and extra_checks.py.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

MOD = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011,
       7: 0b10000011, 8: 0b100011101, 9: 0b1000010001, 10: 0b10000001001}

# palette (validated categorical order) and inks
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
INK, INK2, GRID, PALE = "#0b0b0b", "#52514e", "#dddcd6", "#f0efec"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
})


class GF:
    def __init__(self, k):
        self.k, self.q = k, 1 << k
        q, mod = self.q, MOD[k]
        self.exp = np.zeros(2 * q, dtype=np.int64)
        self.log = np.zeros(q, dtype=np.int64)
        x = 1
        for i in range(q - 1):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & q:
                x ^= mod
        assert x == 1 or k == 1, "modulus not primitive"
        self.exp[q - 1:2 * (q - 1)] = self.exp[:q - 1]
        tr = np.zeros(q, dtype=np.int64)
        for a in range(q):
            s, y = 0, a
            for _ in range(k):
                s ^= y
                y = self.mul(y, y)
            tr[a] = s
        self.chi = 1 - 2 * tr              # chi(x) = (-1)^Tr(x)

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return int(self.exp[self.log[a] + self.log[b]])

    def inv(self, a):
        return int(self.exp[(self.q - 1 - self.log[a]) % (self.q - 1)])

    def tau(self, H):
        """tau(H) for monic H (coefficients low to high), via Newton's identities."""
        d = len(H) - 1
        e = [0] + [H[d - i] for i in range(1, d + 1)]
        p = [0] * (2 * d)
        for m in range(1, 2 * d):
            s = 0
            for i in range(1, min(m - 1, d) + 1):
                s ^= self.mul(e[i], p[m - i])
            if m <= d and m % 2 == 1:
                s ^= e[m]
            p[m] = s
        t = 0
        for j in range(1, d + 1, 2):
            t ^= self.mul(self.mul(H[j], H[j]), p[2 * j - 1])
        return t

    def tau_plane(self, Hf):
        """q x q array T[b1, b0] = tau(H_f + b1 t + b0)."""
        q = self.q
        T = np.zeros((q, q), dtype=np.int64)
        H = list(Hf)
        for b1 in range(q):
            for b0 in range(q):
                H[0], H[1] = b0, b1
                T[b1, b0] = self.tau(H)
        return T

    def chi_div(self, T, a):
        """chi(T / a) elementwise, a != 0."""
        out = np.ones(T.shape, dtype=np.int64)
        nz = T != 0
        out[nz] = self.chi[self.exp[(self.log[T[nz]] - self.log[a]) % (self.q - 1)]]
        return out


def eps(d):
    return 1 if d % 4 in (2, 3) else 0


def Hpoly(d, mid):
    H = [0] * (d + 1)
    H[d] = 1
    for j, c in mid.items():
        H[j] = c
    return H


HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)


# ---------------------------------------------------------------------------
# Figure 1: the plane H_f + P_{<=1} over F_16
# ---------------------------------------------------------------------------
def fig_planes():
    K = GF(4)
    q = K.q
    cases = [
        (3, {2: 2}, r"$n=6$:  $H_f=t^3+\omega t^2$",
         r"bi-affine, $\Delta=0$"),
        (5, {2: 1}, r"$n=10$:  $H_f=t^5+t^2$",
         r"bi-affine, $\Delta\neq0$"),
        (4, {3: 1}, r"$n=8$:  $H_f=t^4+t^3$",
         r"$\mathbb{F}_2$-affine, $T_{0,1}\neq0$"),
        (6, {5: 1, 3: 1, 2: 1}, r"$n=12$:  $H_f=t^6+t^5+t^3+t^2$",
         r"$\mathbb{F}_2$-affine, $T_{0,1}=0$"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.15))
    cm_zero = ListedColormap([PALE, BLUE])
    cm_mu = ListedColormap([ORANGE, BLUE])
    for col, (d, mid, title, kind) in enumerate(cases):
        T = K.tau_plane(Hpoly(d, mid))
        N = int((T == 0).sum())
        s = (-1) ** (K.k * eps(d)) * q * (N - q)
        mu = (-1) ** (K.k * eps(d)) * K.chi_div(T, 2)        # a = omega (the element 2)
        edges = np.arange(q + 1) - 0.5
        for row, (img, cm) in enumerate([((T == 0).astype(int), cm_zero), ((mu > 0).astype(int), cm_mu)]):
            ax = axes[row, col]
            ax.pcolormesh(edges, edges, img, cmap=cm, vmin=0, vmax=1, edgecolors="white", linewidth=0.5)
            ax.set_aspect("equal")
            ax.set_xticks([0, 5, 10, 15])
            ax.set_yticks([0, 5, 10, 15])
            ax.tick_params(length=2, labelsize=7)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if row == 1:
                ax.set_xlabel(r"$b_0$")
            if col == 0:
                ax.set_ylabel(r"$b_1$")
        axes[0, col].set_title(title + "\n" + kind, fontsize=8, color=INK, linespacing=1.35)
        axes[0, col].text(0.5, -0.17, rf"$N(\tau_f)={N}$,  $\sum\mu={s}$", transform=axes[0, col].transAxes,
                          ha="center", va="top", fontsize=8, color=INK2)
    fig.text(0.014, 0.72, r"zero locus of $\tau_f$", rotation=90, va="center", ha="center", fontsize=8.5)
    fig.text(0.014, 0.29, r"$\mu\,(H^2+\omega t)$", rotation=90, va="center", ha="center", fontsize=8.5)
    # legends
    hz = [plt.Rectangle((0, 0), 1, 1, fc=BLUE), plt.Rectangle((0, 0), 1, 1, fc=PALE, ec=GRID, lw=0.5)]
    hm = [plt.Rectangle((0, 0), 1, 1, fc=BLUE), plt.Rectangle((0, 0), 1, 1, fc=ORANGE)]
    fig.legend(hz, [r"$\tau_f=0$", r"$\tau_f\neq0$"], loc="center right", bbox_to_anchor=(1.0, 0.74),
               frameon=False, fontsize=8, handlelength=1, handleheight=1)
    fig.legend(hm, [r"$\mu=+1$", r"$\mu=-1$"], loc="center right", bbox_to_anchor=(1.0, 0.29),
               frameon=False, fontsize=8, handlelength=1, handleheight=1)
    fig.subplots_adjust(left=0.095, right=0.885, top=0.9, bottom=0.1, wspace=0.3, hspace=0.42)
    fig.savefig(os.path.join(FIG, "fig_planes.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: Chowla correlations on I_2(f)
# ---------------------------------------------------------------------------
def correlation(K, Hf, shifts):
    """C_f for shifts (beta0, beta1, a_i): sum over (b0,b1,a) of prod mu(F + alpha_i)."""
    q, d = K.q, len(Hf) - 1
    T = K.tau_plane(Hf)
    b = np.arange(q)
    sgn = (-1) ** (K.k * eps(d) * len(shifts))
    total = 0
    A = set(s[2] for s in shifts)
    for a in range(q):
        if a in A:
            continue
        prod = np.ones((q, q), dtype=np.int64)
        for (be0, be1, ai) in shifts:
            Ts = T[np.ix_(b ^ be1, b ^ be0)]
            prod *= K.chi_div(Ts, a ^ ai)
        total += int(prod.sum())
    return sgn * total


def fig_chowla():
    ks = list(range(2, 8))
    H12 = Hpoly(6, {3: 1})        # n = 12, deg U_f = 1: T01 = 0, tau_f != 0
    H8 = Hpoly(4, {3: 1})         # n = 8, T01 = 1
    series = [
        (r"2 points, shifts $0,\,1$  ($n=12$, $T_{0,1}=0$)", H12, [(0, 0, 0), (1, 0, 0)], BLUE, "o", True),
        (r"4 points, parallelogram $0,\,1,\,t^2,\,t^2+1$  ($n=8$)", H8,
         [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)], AQUA, "D", True),
        (r"2 points, shifts $0,\,t$  ($n=6$, $f=t^6$)", Hpoly(3, {}), [(0, 0, 0), (0, 0, 1)], ORANGE, "s", False),
        (r"2 points, shifts $0,\,t^2+t+1$  ($n=10$, $f=t^{10}+t^4$)", Hpoly(5, {2: 1}),
         [(0, 0, 0), (1, 1, 1)], YELLOW, "^", False),
    ]
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    qs = np.array([2 ** k for k in ks], dtype=float)
    ax.plot(qs, 2 / np.sqrt(qs), color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(qs[-1] * 1.08, 2 / np.sqrt(qs[-1]), r"$2q^{-1/2}$", color=INK2, fontsize=8, va="center")
    data = {}
    for label, Hf, shifts, color, marker, rigid in series:
        vals = []
        for k in ks:
            K = GF(k)
            vals.append(abs(correlation(K, Hf, shifts)) / K.q ** 3)
        vals = np.array(vals)
        data[label] = vals
        ok = vals > 0
        ax.plot(qs[ok], vals[ok], color=color, lw=1.6, marker=marker, ms=7 if marker == "o" else 5,
                mfc="white" if marker == "o" else color, mec=color, mew=1.4, label=label,
                zorder=3 if rigid else 2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(qs)
    ax.set_xticklabels([str(int(x)) for x in qs])
    ax.set_yticks([2.0 ** -e for e in range(0, 13, 2)])
    ax.set_yticklabels([r"$1$"] + [rf"$2^{{-{e}}}$" for e in range(2, 13, 2)])
    ax.set_ylim(2.0 ** -12.5, 1.8)
    ax.set_xlim(qs[0] / 1.2, qs[-1] * 1.9)
    ax.set_xlabel(r"$q=2^k$")
    ax.set_ylabel(r"$|\mathcal{C}_f|\,/\,q^3$")
    ax.grid(True, which="major", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.text(qs[-1] * 1.08, 1.0, r"$1-1/q$", color=INK2, fontsize=8, va="center")
    ax.legend(frameon=False, fontsize=7.5, loc="lower left", handlelength=2.2, labelspacing=0.35)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.14)
    fig.savefig(os.path.join(FIG, "fig_chowla.pdf"))
    plt.close(fig)
    return ks, data


# ---------------------------------------------------------------------------
# Figure 3: polynomial Moebius sums, closed forms of Theorem 1.4(a)
# ---------------------------------------------------------------------------
def fig_polysums():
    ks = list(range(1, 11))
    qs = np.array([2 ** k for k in ks], dtype=float)
    S3 = []
    for k in ks:
        K = GF(k)
        cubes = np.array([K.mul(K.mul(b, b), b) for b in range(K.q)])
        S3.append(int(K.chi[cubes].sum()))            # a = 1
    S3 = np.array(S3, dtype=float)
    sign = np.array([(-1) ** k for k in ks], dtype=float)
    S = {
        2: sign * qs * S3,
        3: sign * qs ** 2,
        4: qs ** 3,
        5: qs ** 3 * S3,
        6: sign * qs ** 4 + qs ** 3 * (qs - 1),
    }
    colors = {2: BLUE, 3: ORANGE, 4: AQUA, 5: YELLOW, 6: MAGENTA}
    markers = {2: "o", 3: "s", 4: "D", 5: "^", 6: "v"}
    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    for sl, lab in [(-1, r"$q^{-1}$"), (-2, r"$q^{-2}$"), (-3, r"$q^{-3}$")]:
        ax.plot(qs, qs ** sl, color=GRID, lw=0.9, zorder=1)
        ax.text(qs[-1] * 1.12, qs[-1] ** sl, lab, color=INK2, fontsize=8, va="center")
    for d in (2, 3, 4, 5, 6):
        v = np.abs(S[d]) / qs ** d
        ok = v > 0
        ax.plot(qs[ok], v[ok], color=colors[d], lw=2.6 if d == 3 else 1.4, marker=markers[d],
                ms=8 if d == 3 else 4.5, mew=1.4,
                mec=colors[d], mfc="white" if d == 3 else colors[d], label=rf"$d={d}$", zorder=2 if d == 3 else 3,
                ls="-" if d in (3, 4) else (0, (1, 1.5)))
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(qs)
    ax.set_xticklabels([rf"$2^{{{k}}}$" for k in ks])
    ax.set_xlim(qs[0] / 1.2, qs[-1] * 2.6)
    ax.set_ylim(2.0 ** -31, 2.5)
    ax.set_yticks([2.0 ** -e for e in range(0, 31, 5)])
    ax.set_yticklabels([r"$1$"] + [rf"$2^{{-{e}}}$" for e in range(5, 31, 5)])
    ax.set_xlabel(r"$q=2^k$")
    ax.set_ylabel(r"$\left|\sum_{H\in\mathcal{M}_d}\mu(H^2+t)\right|\,/\,q^{d}$")
    ax.grid(True, which="major", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower left", ncol=1, handlelength=2.4)
    fig.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.14)
    fig.savefig(os.path.join(FIG, "fig_polysums.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_planes()
    ks, data = fig_chowla()
    for lab, v in data.items():
        print(lab, [round(float(x), 4) for x in v])
    fig_polysums()
    print("figures written to", os.path.abspath(FIG))
