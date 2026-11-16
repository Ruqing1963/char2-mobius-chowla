#!/usr/bin/env python3
"""Supplementary independent checks (direct factorisation) for claims not covered
by verify_mobius_chowla.py: parity law, Thm 1.2 case split, Thm 1.3(a)(b)(c),
Remark 2.2 (n = 16 rigid branch), Thm 1.4(a) cubic sums, Thm 1.4(b)(c),
Thm 1.2(d), Thm 1.3(d), d = 5, 6 closed forms, Prop 2.1 (rigid hypersurface, n = 8m)."""
import itertools, random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify_mobius_chowla as V
from verify_mobius_chowla import (GF2k, poly_mul, poly_add, mobius_poly,
                                  coeff_u_minus_1, tau_direct)
V.IRRED_MOD[6] = 0b1000011  # x^6 + x + 1 (for the cubic-sum sign check at q = 64)

def eps(d): return 1 if d % 4 in (2, 3) else 0
def sgn(k, d): return -1 if (k * eps(d)) & 1 else 1

def UV(Hf, d):
    U = [Hf[2*s+1] for s in range((d-1)//2+1)]
    Vv = [Hf[2*s] for s in range(d//2+1)]
    return U, Vv

def deg(p): return max([i for i, c in enumerate(p) if c], default=-1)

def coeffs_odd(Hf, d, F):
    U, Vv = UV(Hf, d); U2 = poly_mul(U,U,F); U3 = poly_mul(U2,U,F); V2 = poly_mul(Vv,Vv,F)
    den_u = [0] + poly_add([0]+U2, V2); den0 = poly_add([0]+U2, V2)
    c1 = coeff_u_minus_1(poly_mul(U,V2,F), den_u, F)
    c0 = coeff_u_minus_1(poly_mul(V2,Vv,F), den_u, F)
    c00 = coeff_u_minus_1(poly_mul(U3,Vv,F), den0, F)
    return c0, c1, c00 ^ F.mul(c0, c1)

def coeffs_even(Hf, d, F):
    U, Vv = UV(Hf, d); U2 = poly_mul(U,U,F); U3 = poly_mul(U2,U,F); V2 = poly_mul(Vv,Vv,F)
    G0 = poly_add(V2, [0]+U2)
    return (coeff_u_minus_1(poly_mul(U2,Vv,F), G0, F), coeff_u_minus_1(U3, G0, F),
            coeff_u_minus_1(poly_mul(U3,Vv,F), G0, F), Hf[d-1])  # T1, T01, T00, h_{d-1}

def mu_HA(H, a, F):
    return mobius_poly(poly_add(poly_mul(H,H,F), [0, a]), F)

def interval_sum(Hf, F):
    q = F.q; s = 0
    for b0 in range(q):
        for b1 in range(q):
            H = list(Hf); H[0], H[1] = b0, b1
            for a in range(q):
                s += mu_HA(H, a, F)
    return s

fails = []
def check(cond, msg):
    if not cond:
        fails.append(msg); print("  FAIL:", msg)

# ---- 1. parity law
print("[1] parity law mu(H^2+at) = (-1)^{k eps_d} chi(tau(H)/a): d=1..5 q=2,4; d=6 q=2")
for k, dmax in [(1, 6), (2, 5)]:
    F = GF2k(k)
    for d in range(1, dmax+1):
        for c in itertools.product(range(F.q), repeat=d):
            H = list(c) + [1]; t = tau_direct(H, F)
            for a in range(1, F.q):
                check(mu_HA(H, a, F) == sgn(k, d)*F.chi(F.div(t, a)), f"parity k={k} H={H} a={a}")

# ---- 2. Theorem 1.2 case split vs direct sums
print("[2] Thm 1.2 (a)(b)(c) case formulas vs direct factorisation")
for k, ds in [(1, [2,3,4,5,6]), (2, [2,3,4,5])]:
    F = GF2k(k); q = F.q
    for d in ds:
        s = sgn(k, d)
        for h in itertools.product(range(q), repeat=d-2):
            Hf = [0,0] + list(h) + [1]
            act = interval_sum(Hf, F)
            if d == 2:
                pred = 0
            elif d % 2:
                c0, c1, D = coeffs_odd(Hf, d, F)
                if d == 3: check(D == 0, f"Delta!=0 for n=6, h={h}")
                pred = s*q*(q-1) if D == 0 else -s*q
            else:
                T1, T01, T00, hd = coeffs_even(Hf, d, F)
                if T01 or ((T1 != 0) != (hd != 0)): pred = 0
                elif T1 and hd: pred = s*F.chi(F.div(F.mul(hd, T00), F.sq(T1)))*q*q
                elif T00: pred = -s*q*q
                else: pred = s*q*q*(q-1)
            check(act == pred, f"Thm1.2 k={k} d={d} h={h}: actual {act} pred {pred}")

# ---- 3. Theorem 1.3(a)(b): random shift sets, direct correlation sums
print("[3] Thm 1.3(a)(b) Chowla correlations, random shifts; |A|=1 bound q^2 (n>=6), 2q^2 (n=4)")
random.seed(1)
for k, ds in [(1, [2,3,4,5]), (2, [2,3,4])]:
    F = GF2k(k); q = F.q
    for d in ds:
        for trial in range(25 if q == 2 else 12):
            Hf = [0,0] + [random.randrange(q) for _ in range(d-2)] + [1]
            r = random.choice([2,3,4])
            pool = list(itertools.product(range(q), repeat=3))
            if trial % 2 == 0:   # force |A| = 1
                a1 = random.randrange(q); pool = [p for p in pool if p[2] == a1]
            shifts = random.sample(pool, min(r, len(pool)))  # (beta0, beta1, a_i)
            r = len(shifts); A = sorted(set(s_[2] for s_ in shifts))
            C = 0
            for b0 in range(q):
                for b1 in range(q):
                    for a in range(q):
                        p = 1
                        for (be0, be1, ai) in shifts:
                            H = list(Hf); H[0], H[1] = b0 ^ be0, b1 ^ be1
                            p *= mu_HA(H, a ^ ai, F)
                            if p == 0: break
                        C += p
            Theta = {c: {} for c in A}
            for c in A:
                for b0 in range(q):
                    for b1 in range(q):
                        v = 0
                        for (be0, be1, ai) in shifts:
                            if ai == c:
                                H = list(Hf); H[0], H[1] = b0 ^ be0, b1 ^ be1
                                v ^= tau_direct(H, F)
                        Theta[c][(b0, b1)] = v
            allzero = all(v == 0 for c in A for v in Theta[c].values())
            S = (-1)**r if (k*eps(d)) & 1 else 1
            tag = f"k={k} d={d} Hf={Hf} shifts={shifts}"
            if allzero:
                check(C == S*q*q*(q-len(A)), f"1.3(b) {tag}: C={C}")
            elif len(A) == 1:
                N = sum(1 for v in Theta[A[0]].values() if v == 0)
                bound = 2*q*q if d == 2 else q*q
                check(C == S*q*(N-q) and abs(C) <= bound, f"1.3(a) |A|=1 {tag}: C={C}, N={N}")

# ---- 4a. n = 8m+4 family: deg U_f = m = (e-1)/2  =>  T01 = h_{d-1} = 0, T1 = lambda^2
print("[4a] e = 2m+1, deg U_f = m: T01 = h_{d-1} = 0, T1 = lambda^2 != 0 (hence tau_f != 0, sum mu = 0)")
for k, ds in [(1, [6, 10, 14]), (2, [6, 10]), (3, [6])]:
    F = GF2k(k); q = F.q
    for d in ds:
        e = d // 2; m = (e - 1) // 2; cnt = 0
        for h in itertools.product(range(q), repeat=d-2):
            Hf = [0,0] + list(h) + [1]
            U, _ = UV(Hf, d)
            if deg(U) != m: continue
            cnt += 1; lam = U[m]
            T1, T01, T00, hd = coeffs_even(Hf, d, F)
            check(T01 == 0 and hd == 0 and T1 == F.sq(lam), f"4a k={k} d={d} h={h}: T1={T1} T01={T01} hd={hd}")
        print(f"    q={q} n={2*d}: {cnt} f checked")

# ---- 4b. Theorem 1.3(c) (new hypothesis) by direct factorisation
print("[4b] Thm 1.3(c): T01 = 0, tau_f != 0  =>  mu(F+c0) = mu(F), 2-point sum = q^2(q-1)")
for k, ds in [(1, [2, 4, 6, 8]), (2, [2, 4, 6])]:
    F = GF2k(k); q = F.q
    for d in ds:
        cnt = 0
        for h in itertools.product(range(q), repeat=d-2):
            Hf = [0,0] + list(h) + [1]
            if d > 2:
                T1, T01, T00, hd = coeffs_even(Hf, d, F)
                if T01 != 0 or (T1 == 0 and T00 == 0 and hd == 0): continue
            cnt += 1
            if q == 4 and d == 6 and cnt > 40: continue   # cap runtime
            mus = {}
            for b0 in range(q):
                for b1 in range(q):
                    H = list(Hf); H[0], H[1] = b0, b1
                    for a in range(q):
                        mus[(b0, b1, a)] = mu_HA(H, a, F)
            ok = all(mus[(b0, b1, a)] == mus[(0, b1, a)] for (b0, b1, a) in mus)
            check(ok, f"1.3(c) mu(F+c0)!=mu(F) k={k} d={d} h={h}")
            # F + c (c in F_q) = (H + sqrt c)^2 + at ; 2-point sum with c1 != c2
            for c1, c2 in [(0, 1), (1, q-1)] if q > 2 else [(0, 1)]:
                s1, s2 = F.pow(c1, q//2), F.pow(c2, q//2)   # square roots
                tot = sum(mus[(b0 ^ s1, b1, a)] * mus[(b0 ^ s2, b1, a)]
                          for b0 in range(q) for b1 in range(q) for a in range(q))
                check(tot == q*q*(q-1), f"1.3(c) 2-pt k={k} d={d} h={h}: {tot}")
            if d > 2 and (T1 == 0) != (hd == 0):
                check(sum(mus.values()) == 0, f"1.3(c) sum mu != 0 k={k} d={d} h={h}")
            if d > 2 and T1 == 0 and hd == 0:
                check(sum(mus.values()) == -sgn(k, d)*q*q, f"1.3(c) sum mu != -s q^2 k={k} d={d} h={h}")
        print(f"    q={q} n={2*d}: {cnt} intervals with T01=0, tau_f!=0")

# ---- 4c. Remark 2.2: n = 16, deg U_f = 2
print("[4c] Remark 2.2 (n=16, deg U_f=2): T01=0, T1 = lam^2 h6, T00|h6=0 formula, #rigid = (q-1)q^2")
for k in [1, 2, 3]:
    F = GF2k(k); q = F.q; d = 8; rigid = 0
    for h2, h3, h4, h5, h6 in itertools.product(range(q), range(q), range(q), range(1, q), range(q)):
        Hf = [0, 0, h2, h3, h4, h5, h6, 0, 1]   # h7 = 0
        T1, T01, T00, hd = coeffs_even(Hf, d, F)
        lam = h5; l2 = F.sq(lam)
        check(T01 == 0 and T1 == F.mul(l2, h6), f"4c T1/T01 q={q} h={Hf}")
        if h6 == 0:
            pred = F.mul(F.mul(l2, lam), h2 ^ l2) ^ F.mul(h3, F.mul(l2, h4) ^ F.sq(h3))
            check(T00 == pred, f"4c T00 q={q} h={Hf}: {T00} vs {pred}")
        if T1 == T01 == T00 == hd == 0:
            rigid += 1
            # confirm with tau_direct on the whole plane
            for b0 in range(q):
                for b1 in range(q):
                    H = list(Hf); H[0], H[1] = b0, b1
                    check(tau_direct(H, F) == 0, f"4c rigid but tau!=0 q={q} h={Hf}")
    check(rigid == (q-1)*q*q, f"4c rigid count q={q}: {rigid}")
    print(f"    q={q}: {rigid} rigid f with deg U_f = 2 (expected (q-1)q^2 = {(q-1)*q*q})")

# ---- 5. d=2 cubic sums: values {(-1)^{k/2-1} 2 sqrt q, (-1)^{k/2} sqrt q} for k even
print("[5] d=2 cubic character sums sum_b chi(b^3/a)")
for k in [2, 4, 6]:
    F = GF2k(k); q = F.q; r = 1 << (k // 2)
    vals = sorted(set(sum(F.chi(F.div(F.pow(b,3), a)) for b in range(q)) for a in range(1, q)))
    exp = sorted({(-1)**(k//2 - 1) * 2 * r, (-1)**(k//2) * r})
    check(vals == exp, f"cubic sum q={q}: {vals} vs {exp}")
    print(f"    q={q}: values {vals}, predicted {exp}")

# ---- 6. Theorem 1.4(b)(c)
print("[6] Thm 1.4(b) d odd and 1.4(c) d even exact formulas + bounds (direct)")
for k, ds in [(1, [3,4,5,6]), (2, [3,4,5,6])]:
    F = GF2k(k); q = F.q
    for d in ds:
        tables = {}
        for h in itertools.product(range(q), repeat=d-2):
            Hf = [0,0] + list(h) + [1]
            tables[h] = coeffs_odd(Hf, d, F) if d % 2 else coeffs_even(Hf, d, F)
        for a in range(1, q):
            act = 0
            for c in itertools.product(range(q), repeat=d):
                act += mu_HA(list(c)+[1], a, F)
            if d % 2:
                pred = sgn(k,d)*q*sum(F.chi(F.div(D, a)) for (_, _, D) in tables.values())
                bound = q**(d-1)
            else:
                pred = sgn(k,d)*q*q*sum(F.chi(F.div(T00, a)) for (T1, T01, T00, hd) in tables.values()
                                        if T01 == 0 and F.sq(T1) == F.mul(a, hd))
                bound = (d-1)*q**(d-1)
            check(act == pred and abs(act) <= bound, f"Thm1.4 k={k} d={d} a={a}: act {act} pred {pred}")
            print(f"    q={q} d={d} a={a}: sum = {act}, formula = {pred}, bound = {bound}")

# ---- 7. very short intervals I_0(f) = f + F_q and I_1(f) = f + P_{<=1} (Thm 1.2(d))
print("[7] Thm 1.2(d): I_0(f), I_1(f) dichotomy / trichotomy (direct factorisation)")
def tau01(Hf, d, b0, b1, F):
    """(tau_0, tau_1) with tau(H_0 + B_0) = tau_0 + tau_1 B_0, H_0 = H_f + b1 t + b0."""
    if d == 2:
        return F.pow(b1, 3), 0
    if d % 2:
        c0, c1, D = coeffs_odd(Hf, d, F)
        return F.mul(b0 ^ c0, b1 ^ c1) ^ D, b1 ^ c1
    T1, T01, T00, hd = coeffs_even(Hf, d, F)
    return F.mul(T01, b0) ^ F.mul(T1, b1) ^ F.mul(hd, F.sq(b1)) ^ T00, T01
for k, ds in [(1, [2,3,4,5,6]), (2, [2,3,4,5]), (3, [2,3,4])]:
    F = GF2k(k); q = F.q
    for d in ds:
        s = sgn(k, d); n0 = 0
        for h in itertools.product(range(q), repeat=d-2):
            if q == 8 and d == 4 and n0 >= 64: break
            n0 += 1
            Hf = [0,0] + list(h) + [1]
            for b0 in range(q):
                for b1 in range(q):
                    H0 = list(Hf); H0[0], H0[1] = b0, b1
                    t0, t1 = tau01(Hf, d, b0, b1, F)
                    # tau(H_0 + B_0) = t0 + t1 B_0 for all B_0
                    for B0 in range(q):
                        H = list(H0); H[0] ^= B0
                        check(tau_direct(H, F) == t0 ^ F.mul(t1, B0), f"tau affine k={k} d={d} H0={H0}")
                    H02 = poly_mul(H0, H0, F)
                    mu = {(a, c): mobius_poly(poly_add(H02, [c, a]), F) for a in range(q) for c in range(q)}
                    S1 = sum(mu.values())
                    p1 = 0 if t1 else (s*q*(q-1) if t0 == 0 else -s*q)
                    check(S1 == p1, f"I_1 k={k} d={d} H0={H0}: {S1} vs {p1}")
                    for a1 in range(1, q):
                        S0 = sum(mu[(a1, c)] for c in range(q))
                        p0 = 0 if t1 else s*F.chi(F.div(t0, a1))*q
                        check(S0 == p0, f"I_0 k={k} d={d} H0={H0} a1={a1}: {S0} vs {p0}")
        print(f"    q={q} n={2*d}: {n0} H_f checked (all b0, b1, a1)")

# ---- 8. 4-point parallelogram Chowla sums, n = 0 mod 4, n >= 8, arbitrary f
print("[8] 4-point parallelogram shifts (n = 0 mod 4, n >= 8): C = q^2(q-1)")
random.seed(2)
for k, ds, trials in [(1, [4, 6, 8], 30), (2, [4, 6], 25), (3, [4], 6)]:
    F = GF2k(k); q = F.q
    for d in ds:
        nT01 = 0
        for _ in range(trials):
            Hf = [0,0] + [random.randrange(q) for _ in range(d-2)] + [1]
            if coeffs_even(Hf, d, F)[1]: nT01 += 1
            pts = list(itertools.product(range(q), repeat=2))
            while True:
                x, u, v = random.sample(pts, 3)
                u = (u[0] ^ x[0], u[1] ^ x[1]); v = (v[0] ^ x[0], v[1] ^ x[1])
                quad = [x, (x[0]^u[0], x[1]^u[1]), (x[0]^v[0], x[1]^v[1]), (x[0]^u[0]^v[0], x[1]^u[1]^v[1])]
                if len(set(quad)) == 4: break
            c = random.randrange(q)
            mus = {}
            for b0 in range(q):
                for b1 in range(q):
                    H = list(Hf); H[0], H[1] = b0, b1
                    H2 = poly_mul(H, H, F)
                    for a in range(q):
                        mus[(b0, b1, a)] = mobius_poly(poly_add(H2, [0, a]), F)
            C = 0
            for b0 in range(q):
                for b1 in range(q):
                    for a in range(q):
                        p = 1
                        for (be0, be1) in quad:
                            p *= mus[(b0 ^ be0, b1 ^ be1, a ^ c)]
                        C += p
            check(C == q*q*(q-1), f"4-pt k={k} d={d} Hf={Hf} quad={quad} c={c}: {C}")
        print(f"    q={q} n={2*d}: {trials} random (f, parallelogram, c), {nT01} with T01 != 0")

# ---- 9. d = 5, 6 closed forms
print("[9] d=5: Delta = (h2+h3h4)^3, sum = q^3 sum_b chi(b^3/a);  d=6: T01 = h5^2(h3+h5^3), "
      "sum = (-1)^k q^4 + q^3(q-1)")
for k in [1, 2, 3]:
    F = GF2k(k); q = F.q
    for h2, h3, h4 in itertools.product(range(q), repeat=3):
        c0, c1, D = coeffs_odd([0, 0, h2, h3, h4, 1], 5, F)
        check(c1 == F.mul(F.sq(h4), F.sq(h4) ^ h3) and c0 == F.mul(F.sq(h4), F.pow(h4, 3) ^ h2)
              and D == F.pow(h2 ^ F.mul(h3, h4), 3), f"d=5 coeffs q={q} h={(h2,h3,h4)}")
    for h2, h3, h4, h5 in itertools.product(range(q), repeat=4):
        T1, T01, T00, hd = coeffs_even([0, 0, h2, h3, h4, h5, 1], 6, F)
        check(T01 == F.mul(F.sq(h5), h3 ^ F.pow(h5, 3)), f"d=6 T01 q={q} h={(h2,h3,h4,h5)}")
        check(T1 == F.sq(h3) ^ F.mul(h2, F.sq(h5)) ^ F.mul(h4, F.pow(h5, 4))
              ^ F.mul(F.sq(h4), F.sq(h5)) ^ F.pow(h5, 6), f"d=6 T1 q={q} h={(h2,h3,h4,h5)}")
        if h5 and h3 == F.pow(h5, 3):   # locus of the d=6 proof: T00 = c^3 (w^2 + y w + y^2)
            c = h5; w = h2 ^ F.mul(h4, F.sq(c)) ^ F.sq(h4); y = F.pow(c, 4) ^ F.sq(h4)
            check(T00 == F.mul(F.pow(c, 3), F.sq(w) ^ F.mul(y, w) ^ F.sq(y)),
                  f"d=6 T00 on locus q={q} h={(h2,h3,h4,h5)}")
    for a in range(1, q):
        cub = sum(F.chi(F.div(F.pow(b, 3), a)) for b in range(q))
        pred5 = q**3 * cub
        pred6 = (-1)**k * q**4 + q**3 * (q - 1)
        for d, pred in [(5, pred5), (6, pred6)]:
            if q == 8 and not (d == 5 and a == 1):   # parity law (checked in [1]) for the rest
                act = sum(sgn(k, d)*F.chi(F.div(tau_direct(list(c)+[1], F), a))
                          for c in itertools.product(range(q), repeat=d))
                how = "parity"
            else:
                act = sum(mu_HA(list(c)+[1], a, F) for c in itertools.product(range(q), repeat=d))
                how = "direct"
            check(act == pred, f"d={d} closed form q={q} a={a}: {act} vs {pred}")
            print(f"    q={q} d={d} a={a}: sum = {act} ({how}), formula = {pred}")

# ---- 10. Proposition 2.1: V_f = lam^-2 U_f^2 + lam^2 u + W_0, deg U_f = m, deg W_0 <= m
print("[10] Prop 2.1: h_{d-1} = T01 = T1 = 0, T00 = lam^2 [u^-1] U R/G0; = lam^3 w_{m-1} if deg W0 <= m-1")
random.seed(3)
for k, ms, trials in [(2, [2, 3], 120), (3, [2, 3, 4], 80), (4, [2, 3, 4, 5], 60)]:
    F = GF2k(k); q = F.q
    for m in ms:
        d = 4*m; nrig = 0
        for trial in range(trials):
            lam = random.randrange(1, q)
            U = [0] + [random.randrange(q) for _ in range(m-1)] + [lam]
            wdeg = random.choice([m-2, m-1, m]) if m > 2 else random.choice([m-1, m])
            W = [0] + [random.randrange(q) for _ in range(m)]
            for j in range(wdeg+1, m+1): W[j] = 0
            il2 = F.div(1, F.sq(lam))
            U2 = poly_mul(U, U, F)
            Vv = [F.mul(il2, c) for c in U2]
            Vv[1] ^= F.sq(lam)
            for j, c in enumerate(W): Vv[j] ^= c
            Hf = [0]*(d+1)
            for j, c in enumerate(Vv): Hf[2*j] = c
            for j, c in enumerate(U): Hf[2*j+1] = c
            T1, T01, T00, hd = coeffs_even(Hf, d, F)
            R = poly_add(poly_add([F.mul(il2, c) for c in poly_mul(U2, W, F)], poly_mul(W, W, F)),
                         [0, 0, F.pow(lam, 4)])
            G0 = poly_add(poly_mul(Vv, Vv, F), [0] + U2)
            pred = F.mul(F.sq(lam), coeff_u_minus_1(poly_mul(U, R, F), G0, F))
            ok = hd == 0 and T01 == 0 and T1 == 0 and T00 == pred
            if wdeg <= m-1:
                ok = ok and T00 == F.mul(F.pow(lam, 3), W[m-1])
            check(ok, f"Prop2.1 q={q} m={m} U={U} W={W}: T1={T1} T01={T01} T00={T00} pred={pred}")
            if T00 == 0:
                nrig += 1
                if q <= 4 and m == 2:   # confirm rigidity on the whole plane
                    for b0 in range(q):
                        for b1 in range(q):
                            H = list(Hf); H[0], H[1] = b0, b1
                            check(tau_direct(H, F) == 0, f"Prop2.1 rigid but tau!=0 q={q} Hf={Hf}")
        print(f"    q={q} n={8*m}: {trials} random (U_f, W_0), {nrig} rigid")

# ---- 11. Remark 3.2: d even, exact evaluation of Chowla sums by summing over (b0, b1) first
print("[11] Remark 3.2 (d even): C_f = sgn q^2 sum_{a: psi T01 = 0, psi^2 T1^2 = psi h} chi(kappa(a)); bounds")
random.seed(4)
for k, ds, trials in [(2, [4, 6], 150), (3, [4, 6], 120), (4, [4], 40)]:
    F = GF2k(k); q = F.q
    for d in ds:
        ncase = [0, 0, 0]; nzero2 = 0
        for trial in range(trials):
            Hf = [0,0] + [random.randrange(q) for _ in range(d-2)] + [1]
            if trial % 3 == 0:                  # force T01 = 0 fairly often (h_{d-1} = 0, small U_f)
                for j in range(3, d, 2): Hf[j] = 0
                Hf[3] = random.randrange(q)
            tau = [[0]*q for _ in range(q)]
            for b0 in range(q):
                for b1 in range(q):
                    H = list(Hf); H[0], H[1] = b0, b1
                    tau[b0][b1] = tau_direct(H, F)
            T00 = tau[0][0]; T01 = tau[1][0] ^ T00; hd = Hf[d-1]; T1 = tau[0][1] ^ T00 ^ hd
            r = random.choice([2, 2, 3, 4, 5])
            shifts = random.sample(list(itertools.product(range(q), repeat=3)), r)   # (beta0, beta1, a_i)
            A = sorted(set(s_[2] for s_ in shifts))
            rc = {c: sum(1 for s_ in shifts if s_[2] == c) for c in A}
            Aodd = [c for c in A if rc[c] % 2]
            Th0 = {c: 0 for c in A}
            for (be0, be1, ai) in shifts: Th0[ai] ^= tau[be0][be1]
            S = (-1)**r if (k*eps(d)) & 1 else 1
            C = 0
            for a in range(q):
                if a in A: continue
                for b0 in range(q):
                    for b1 in range(q):
                        v = 0
                        for (be0, be1, ai) in shifts:
                            v ^= F.div(tau[b0 ^ be0][b1 ^ be1], a ^ ai)
                        C += F.chi(v)
            C *= S
            pred = 0; ncontrib = 0
            for a in range(q):
                if a in A: continue
                psi = 0; kap = 0
                for c in Aodd: psi ^= F.div(1, a ^ c)
                for c in A: kap ^= F.div(Th0[c], a ^ c)
                if F.mul(psi, T01) == 0 and F.mul(F.sq(psi), F.sq(T1)) == F.mul(psi, hd):
                    pred += F.chi(kap); ncontrib += 1
            pred *= S*q*q
            tag = f"k={k} d={d} Hf={Hf} shifts={shifts}"
            check(C == pred, f"Rem3.2 formula {tag}: C={C} pred={pred}")
            m = len(Aodd)
            if m >= 1:
                if T01 != 0 or ((T1 != 0) != (hd != 0)):
                    ncase[0] += 1
                    check(abs(C) <= ((m-1)//2)*q*q, f"Rem3.2 bound (i) {tag}: C={C}")
                    if r == 2 and len(A) == 2:
                        nzero2 += 1; check(C == 0, f"Rem3.2 2-point {tag}: C={C}")
                elif T1 != 0 and hd != 0:
                    ncase[1] += 1
                    check(abs(C) <= (m + (m-1)//2)*q*q, f"Rem3.2 bound (ii) {tag}: C={C}")
                else:
                    ncase[2] += 1
        print(f"    q={q} n={2*d}: {trials} random (f, shifts); cases (i)/(ii)/const = {ncase}, "
              f"{nzero2} two-point sums with a_1 != a_2 all 0")

print("\nFAILURES:" if fails else "\nALL EXTRA CHECKS PASSED")
for f in fails[:30]: print("  -", f)
