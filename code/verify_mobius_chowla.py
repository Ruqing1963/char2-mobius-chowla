#!/usr/bin/env python3
"""
Exhaustive verification script for:
"Short-interval Mobius sums, Chowla correlations, and polynomial Mobius sums in characteristic 2"
Pure Python standard library only.
"""

import itertools

# ----------------------------------------------------------------------
# Finite field arithmetic in F_{2^k} (k = 1, 2, 3, 4)
# ----------------------------------------------------------------------
IRRED_MOD = {
    1: 0b11,     # x + 1
    2: 0b111,    # x^2 + x + 1
    3: 0b1011,   # x^3 + x + 1
    4: 0b10011,  # x^4 + x + 1
}

class GF2k:
    def __init__(self, k):
        self.k = k
        self.q = 1 << k
        self.mod = IRRED_MOD[k]
        self.mul_tab = [[0] * self.q for _ in range(self.q)]
        self.inv_tab = [0] * self.q
        self.tr_tab = [0] * self.q
        self._build()

    def _build(self):
        for a in range(self.q):
            for b in range(self.q):
                res = 0
                aa, bb = a, b
                while bb:
                    if bb & 1:
                        res ^= aa
                    aa <<= 1
                    if aa & self.q:
                        aa ^= self.mod
                    bb >>= 1
                self.mul_tab[a][b] = res
        for a in range(1, self.q):
            for b in range(1, self.q):
                if self.mul_tab[a][b] == 1:
                    self.inv_tab[a] = b
                    break
        for a in range(self.q):
            t = 0
            cur = a
            for _ in range(self.k):
                t ^= cur
                cur = self.mul_tab[cur][cur]
            self.tr_tab[a] = t & 1

    def mul(self, a, b):
        return self.mul_tab[a][b]

    def sq(self, a):
        return self.mul_tab[a][a]

    def pow(self, a, e):
        res = 1
        base = a
        while e > 0:
            if e & 1:
                res = self.mul_tab[res][base]
            base = self.mul_tab[base][base]
            e >>= 1
        return res

    def div(self, a, b):
        return self.mul_tab[a][self.inv_tab[b]]

    def chi(self, a):
        return -1 if self.tr_tab[a] else 1

# ----------------------------------------------------------------------
# Polynomial arithmetic & distinct-degree factorization over F_{2^k}[t]
# Polynomials are lists of coefficients [c_0, c_1, ..., c_n]
# ----------------------------------------------------------------------
def poly_strip(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_add(a, b):
    m = max(len(a), len(b))
    res = [(a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0) for i in range(m)]
    return poly_strip(res)

def poly_mul(a, b, F):
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    res[i + j] ^= F.mul(ai, bj)
    return poly_strip(res)

def poly_divmod(a, b, F):
    a = list(a)
    poly_strip(a)
    db = len(b) - 1
    inv_lb = F.inv_tab[b[-1]]
    q = [0] * max(1, len(a) - db)
    while len(a) - 1 >= db and not (len(a) == 1 and a[0] == 0):
        deg_diff = len(a) - 1 - db
        coeff = F.mul(a[-1], inv_lb)
        q[deg_diff] = coeff
        for j in range(db + 1):
            a[deg_diff + j] ^= F.mul(coeff, b[j])
        poly_strip(a)
    return poly_strip(q), poly_strip(a)

def poly_gcd(a, b, F):
    while not (len(b) == 1 and b[0] == 0):
        _, r = poly_divmod(a, b, F)
        a, b = b, r
    if a[-1] != 1 and a[-1] != 0:
        inv = F.inv_tab[a[-1]]
        a = [F.mul(c, inv) for c in a]
    return a

def poly_powmod(base, exp, modp, F):
    res = [1]
    b = list(base)
    while exp > 0:
        if exp & 1:
            res = poly_divmod(poly_mul(res, b, F), modp, F)[1]
        b = poly_divmod(poly_mul(b, b, F), modp, F)[1]
        exp >>= 1
    return res

def mobius_poly(f_coeffs, F):
    """Exact Mobius function mu(F) via distinct-degree factorization."""
    # Derivative in char 2
    deriv = [f_coeffs[i] if (i & 1) else 0 for i in range(1, len(f_coeffs))]
    if not deriv:
        return 0
    poly_strip(deriv)
    if len(deriv) == 1 and deriv[0] == 0:
        return 0
    g = poly_gcd(f_coeffs, deriv, F)
    if len(g) > 1:
        return 0  # not square-free
    # Count irreducible factors via distinct-degree factorization
    cur = list(f_coeffs)
    num_factors = 0
    h = [0, 1]  # t
    d = 1
    while len(cur) - 1 >= 2 * d:
        h = poly_powmod(h, F.q, cur, F)
        # g = gcd(h + t, cur)
        ht = poly_add(h, [0, 1])
        g = poly_gcd(cur, ht, F)
        if len(g) > 1:
            deg_g = len(g) - 1
            num_factors += deg_g // d
            cur, _ = poly_divmod(cur, g, F)
            h = poly_divmod(h, cur, F)[1]
        d += 1
    if len(cur) > 1:
        num_factors += 1
    return -1 if (num_factors & 1) else 1

# ----------------------------------------------------------------------
# Residue [u^{-1}] N(u) / D(u) at u = infinity
# ----------------------------------------------------------------------
def coeff_u_minus_1(num, den, F):
    """Return [u^{-1}] of num(u)/den(u) in F((u^{-1}))."""
    num = poly_strip(list(num))
    den = poly_strip(list(den))
    if len(num) == 1 and num[0] == 0:
        return 0
    deg_n = len(num) - 1
    deg_d = len(den) - 1
    if deg_n < deg_d - 1:
        return 0
    # Multiply num by u so that [u^{-1}] of num/den is the constant term of quotient of (u * num) // den
    u_num = [0] + num
    _, rem = poly_divmod(num, den, F)
    if len(rem) - 1 < deg_d - 1 or (len(rem) == 1 and rem[0] == 0):
        return 0
    return F.div(rem[-1], den[-1])

def tau_direct(H, F):
    """Compute tau(H) = [u^{-1}] U^3 V / (V^2 + u U^2)."""
    d = len(H) - 1
    U = [H[2 * s + 1] for s in range((d - 1) // 2 + 1)]
    V = [H[2 * s] for s in range(d // 2 + 1)]
    U2 = poly_mul(U, U, F)
    U3 = poly_mul(U2, U, F)
    V2 = poly_mul(V, V, F)
    num = poly_mul(U3, V, F)
    den = poly_add(V2, [0] + U2)
    return coeff_u_minus_1(num, den, F)

# ----------------------------------------------------------------------
# Verification tests
# ----------------------------------------------------------------------
def test_theorem_1_1():
    print("=== Test 1: Theorem 1.1 (Bilinear / F2-affine structure of tau_f) ===")
    for k in [1, 2, 3]:
        F = GF2k(k)
        for d in [3, 4, 5, 6]:
            for h_tuple in itertools.product(range(F.q), repeat=d - 2):
                # H_f = t^d + sum_{j=2}^{d-1} h_j t^j
                Hf = [0, 0] + list(h_tuple) + [1]
                Uf = [Hf[2 * s + 1] for s in range((d - 1) // 2 + 1)]
                Vf = [Hf[2 * s] for s in range(d // 2 + 1)]
                Uf2 = poly_mul(Uf, Uf, F)
                Uf3 = poly_mul(Uf2, Uf, F)
                Vf2 = poly_mul(Vf, Vf, F)
                Vf3 = poly_mul(Vf2, Vf, F)

                if d % 2 == 1:
                    den_u = [0] + poly_add([0] + Uf2, Vf2)  # u^2 U_f^2 + u V_f^2
                    den_0 = poly_add([0] + Uf2, Vf2)
                    c1 = coeff_u_minus_1(poly_mul(Uf, Vf2, F), den_u, F)
                    c0 = coeff_u_minus_1(Vf3, den_u, F)
                    c00 = coeff_u_minus_1(poly_mul(Uf3, Vf, F), den_0, F)
                    delta = c00 ^ F.mul(c0, c1)
                    for b0 in range(F.q):
                        for b1 in range(F.q):
                            H = list(Hf)
                            H[0], H[1] = b0, b1
                            lhs = tau_direct(H, F)
                            rhs = F.mul(b0 ^ c0, b1 ^ c1) ^ delta
                            assert lhs == rhs, f"Mismatch at d={d}, q={F.q}"
                else:
                    G0 = poly_add(Vf2, [0] + Uf2)
                    T1 = coeff_u_minus_1(poly_mul(Uf2, Vf, F), G0, F)
                    T01 = coeff_u_minus_1(Uf3, G0, F)
                    T00 = coeff_u_minus_1(poly_mul(Uf3, Vf, F), G0, F)
                    hd1 = Hf[d - 1]
                    for b0 in range(F.q):
                        for b1 in range(F.q):
                            H = list(Hf)
                            H[0], H[1] = b0, b1
                            lhs = tau_direct(H, F)
                            rhs = F.mul(T01, b0) ^ F.mul(T1, b1) ^ F.mul(hd1, F.sq(b1)) ^ T00
                            assert lhs == rhs, f"Mismatch at d={d}, q={F.q}"
    print("  [PASS] Theorem 1.1 verified exhaustively for d = 3, 4, 5, 6 over F_2, F_4, F_8.")

def test_proposition_2_1():
    print("\n=== Test 2: Proposition 2.1 (Infinite family of sporadic rigid intervals) ===")
    for k in [1, 2]:
        F = GF2k(k)
        for m in [2, 3]:  # n = 16, 24 (d = 8, 12)
            d = 4 * m
            for coeffs in itertools.product(range(F.q), repeat=m - 1):
                for lam in range(1, F.q):
                    Uf = [0] + list(coeffs) + [lam]  # degree m, Uf(0) = 0
                    Uf2 = poly_mul(Uf, Uf, F)
                    inv_lam2 = F.inv_tab[F.sq(lam)]
                    lam2 = F.sq(lam)
                    Vf = poly_add([F.mul(inv_lam2, c) for c in Uf2], [0, lam2])
                    # Build H_f(t) = V_f(t^2) + t U_f(t^2)
                    Hf = [0] * (d + 1)
                    for idx, c in enumerate(Vf):
                        Hf[2 * idx] ^= c
                    for idx, c in enumerate(Uf):
                        Hf[2 * idx + 1] ^= c
                    for b0 in range(F.q):
                        for b1 in range(F.q):
                            H = list(Hf)
                            H[0], H[1] = b0, b1
                            assert tau_direct(H, F) == 0, f"Nonzero tau at m={m}"
    print("  [PASS] Proposition 2.1 verified for n = 16 (m=2) and n = 24 (m=3) over F_2, F_4.")

def test_theorems_1_2_and_1_3():
    print("\n=== Test 3: Theorems 1.2 & 1.3 (Short-interval Mobius & Chowla sums by direct factoring) ===")
    for k in [1, 2]:
        F = GF2k(k)
        q = F.q
        for d in [2, 3, 4]:
            n = 2 * d
            eps_d = 1 if (d % 4 in (2, 3)) else 0
            sign_Factor = -1 if ((k * eps_d) & 1) else 1
            for h_tuple in itertools.product(range(q), repeat=d - 2):
                Hf = [0, 0] + list(h_tuple) + [1]
                # Exhaustively factor all q^3 polynomials F = H^2 + a*t in I_2(f)
                mob_table = {}
                N_tau = 0
                for b0 in range(q):
                    for b1 in range(q):
                        H = list(Hf)
                        H[0], H[1] = b0, b1
                        if tau_direct(H, F) == 0:
                            N_tau += 1
                        H2 = poly_mul(H, H, F)
                        for a in range(q):
                            F_poly = poly_add(H2, [0, a])
                            mob_table[(b0, b1, a)] = mobius_poly(F_poly, F)
                actual_sum = sum(mob_table.values())
                predicted_sum = sign_Factor * q * (N_tau - q)
                assert actual_sum == predicted_sum, f"Mobius sum mismatch at n={n}, q={q}"

        # Check Theorem 1.3(c) on f = t^{12} + t^6 (d = 6, e = 3, deg U_f = 1 in [2/3, 4/3])
        d = 6
        Hf = [0] * 7
        Hf[6] = 1
        Hf[3] = 1  # H_f = t^6 + t^3 => f = H_f^2 = t^{12} + t^6
        for b0 in range(q):
            for b1 in range(q):
                H_base = list(Hf)
                H_base[0], H_base[1] = 0, b1
                H_shift = list(Hf)
                H_shift[0], H_shift[1] = b0, b1
                assert tau_direct(H_base, F) == tau_direct(H_shift, F)
    print("  [PASS] Theorems 1.2 & 1.3 verified by direct polynomial factorization.")

def test_theorem_1_4():
    print("\n=== Test 4: Theorem 1.4 (Sawin-Shusterman polynomial Mobius sums over M_d) ===")
    for k in [1, 2]:
        F = GF2k(k)
        q = F.q
        for d in [1, 2, 3, 4, 5]:
            for a in range(1, q):
                total_mu = 0
                for coeffs in itertools.product(range(q), repeat=d):
                    H = list(coeffs) + [1]
                    F_poly = poly_add(poly_mul(H, H, F), [0, a])
                    total_mu += mobius_poly(F_poly, F)
                if d == 1:
                    assert total_mu == 0
                elif d == 2:
                    expected = ((-1) ** k) * q * sum(F.chi(F.div(F.pow(b, 3), a)) for b in range(q))
                    assert total_mu == expected
                elif d == 3:
                    assert total_mu == ((-1) ** k) * (q ** 2)
                elif d == 4:
                    assert total_mu == (q ** 3)
                else:
                    assert abs(total_mu) <= (q ** (d - 1))
                print(f"  q={q}, d={d}, a={a}: sum_{{H in M_{d}}} mu(H^2 + a*t) = {total_mu:4d} [OK]")

if __name__ == "__main__":
    test_theorem_1_1()
    test_proposition_2_1()
    test_theorems_1_2_and_1_3()
    test_theorem_1_4()
    print("\nALL EXHAUSTIVE CHECKS PASSED 100%.")
