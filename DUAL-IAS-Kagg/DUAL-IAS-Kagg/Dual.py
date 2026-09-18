
import io, math, argparse, contextlib
from math import sqrt, pi, exp, log2, comb, ceil, floor

from MSIS_security import SIS_l2_cost
from MLWE_security import MLWE_optimize_attack, LWE_primal_cost, LWE_dual_cost
from model_BKZ     import svp_classical, svp_quantum, svp_plausible

SQRT_2PI = sqrt(2 * pi)
GAUSS_H  = 0.5 * log2(2 * pi * exp(1))    # bits/coeff to encode a Gaussian 
LOG_INF  = 9999
SVP = {'classical': svp_classical, 'quantum': svp_quantum, 'plausible': svp_plausible}


# ----------------------------------------------------------------------------
# prime:  smallest q > 2^logq with q = 5 (mod 8)  
# ----------------------------------------------------------------------------
def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def get_prime(logq):
    q = 2 ** logq + 1
    while not (_is_prime(q) and q % 8 == 5):
        q += 2
    return q


def _quiet(fn, *a, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)



# ----------------------------------------------------------------------------
# MSIS
# ----------------------------------------------------------------------------
def _msis_bmin(q, w, h, B):
    if B >= q:
        return None
    ok = lambda b: _quiet(SIS_l2_cost, q, w, h, B, b, svp_classical) < LOG_INF
    lo, hi = 50, w
    if not ok(hi):
        return None
    while lo < hi:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo


def msis_bits(ps):
    b = _msis_bmin(ps.q, ps.msis_w * ps.N, ps.msis_k * ps.N, ps.beta)
    if b is None:
        return 0, 0, 0, 0
    return b, int(svp_classical(b)), int(svp_quantum(b)), int(svp_plausible(b))


# ----------------------------------------------------------------------------
# MLWE : 
#   
# ----------------------------------------------------------------------------
def mlwe_bits(N, dim, k, eta, q):
    s = sqrt(eta * (eta + 1) / 3.0)
    nd, mm = N * dim, N * k
    out = {}
    for tag, model in SVP.items():
        _, _, cp = _quiet(MLWE_optimize_attack, q, nd, mm, s, LWE_primal_cost, model, False)
        _, _, cd = _quiet(MLWE_optimize_attack, q, nd, mm, s, LWE_dual_cost,   model, False)
        out[tag] = min(int(floor(cp)), int(floor(cd)))
    return out


def mlwe_primal(N, dim, k, eta, q, model='classical'):
    s = sqrt(eta * (eta + 1) / 3.0)
    _, _, cp = _quiet(MLWE_optimize_attack, q, N * dim, N * k, s,
                      LWE_primal_cost, SVP[model], False)
    return int(floor(cp))


# ----------------------------------------------------------------------------
# Scheme 2 parameter set 
# ----------------------------------------------------------------------------
class ParamSet:
    def __init__(self, N, k, l, lp, logq, eta, n,
                 kappa=60, gamma=1.1, t=13.5, etap=None):
        self.N, self.k, self.l, self.lp, self.logq = N, k, l, lp, logq
        self.eta  = eta
        self.etap = eta if etap is None else etap
        self.n    = n
        self.kappa, self.gamma, self.t = kappa, gamma, t
        self.q = get_prime(logq)

        m, mp = l + k, lp + k
        self.m, self.mp = m, mp
        self.alpha = 8.5 * n                        

        
        self.chal = kappa ** 2

        self.T  = self.chal * self.eta  * sqrt(N * m)    
        self.Tp = self.chal * self.etap * sqrt(N * mp)   

        lq = log2(self.q)
        self.s1  = SQRT_2PI * self.alpha * self.T                       
        self.s2  = 2 ** (1 + log2(N) + lq * (k / m  + 2 / (N * m)))     
        self.sp1 = SQRT_2PI * self.alpha * self.Tp                      
        self.sp2 = 2 ** (1 + log2(N) + lq * (k / mp + 2 / (N * mp)))    
        self.s   = max(self.s1, self.s2)
        self.sp  = max(self.sp1, self.sp2)
        self.log2_s,  self.s_bind  = log2(self.s),  ('rej' if self.s1  >= self.s2  else 'reg')
        self.log2_sp, self.sp_bind = log2(self.sp), ('rej' if self.sp1 >= self.sp2 else 'reg')

        self.M  = exp(t / self.alpha + 1 / (2 * self.alpha ** 2))
        self.Mn = self.M ** n

        self.ew_cur = sqrt(self.s ** 2 + self.sp ** 2)                 
        self.ew_opt = max(sqrt(self.s1 ** 2 + self.sp2 ** 2),          
                          sqrt(self.sp1 ** 2 + self.s2 ** 2))

        self.Bn   = gamma * (self.s      / SQRT_2PI) * sqrt(N * l ) * sqrt(n)
        self.Bpn  = gamma * (self.sp     / SQRT_2PI) * sqrt(N * lp) * sqrt(n)
        self.Bppn = gamma * (self.ew_cur / SQRT_2PI) * sqrt(N * k ) * sqrt(n)
        self.beta = sqrt(self.Bn ** 2 + self.Bpn ** 2 + self.Bppn ** 2)  

        self.log2_C = log2(comb(N, kappa)) + kappa        

        self.msis_w = 1 + l + lp + k       
        self.msis_k = k
        self.mlwe_pk  = (l,      k, self.eta )
        self.mlwe_sim = (lp - 1, k, self.etap)

        # ---- sizes -----------------------------------------------------------
        qbits = ceil(log2(self.q))
        self.sig_z = self.s      * sqrt(n) / SQRT_2PI
        self.sig_r = self.sp     * sqrt(n) / SQRT_2PI
        self.sig_e = self.ew_cur * sqrt(n) / SQRT_2PI

        def _resp(sig_e):
            return (N * l  * (log2(self.sig_z) + GAUSS_H)
                  + N * lp * (log2(self.sig_r) + GAUSS_H)
                  + N * k  * (log2(sig_e)      + GAUSS_H))
        resp       = _resp(self.sig_e)
        resp_appxC = _resp(self.ew_opt * sqrt(n) / SQRT_2PI)

        self.anchor_bits, self.anchor_kind = 32 * 8, 'c'

        self.pk_kB        = (k * N * qbits) / 8 / 1024
        self.resp_kB      = resp / 8 / 1024
        self.sig_kB       = (resp + self.anchor_bits) / 8 / 1024
        self.tot_kB       = self.pk_kB + self.sig_kB
        self.sig_kB_appxC = (resp_appxC + self.anchor_bits) / 8 / 1024
        self.tot_kB_appxC = self.pk_kB + self.sig_kB_appxC

    def valid(self):
        return (self.q % 8 == 5) and (self.log2_C > 256) and (self.beta < self.q)


# ----------------------------------------------------------------------------
def report(ps, name, run_mlwe=False, fast_mlwe=True, appxC=False):
    print('-' * 78)
    print(f'  {name}   [Scheme 2]')
    print(f"  N={ps.N} k={ps.k} l={ps.l} l'={ps.lp} logq={ps.logq} eta={ps.eta} "
          f'n={ps.n} kappa={ps.kappa}')
    print(f'  q={ps.q}  (q mod 8={ps.q % 8})   alpha={ps.alpha:.0f}  '
          f'M={ps.M:.4f}  M^n={ps.Mn:.2f}')
    print(f"  log2(s)={ps.log2_s:.2f} [{ps.s_bind}]   log2(s')={ps.log2_sp:.2f} [{ps.sp_bind}]")
    print(f"  log2(Bn)={log2(ps.Bn):.2f}  log2(B'n)={log2(ps.Bpn):.2f}  "
          f"log2(B''n)={log2(ps.Bppn):.2f}")
    print(f'  log2|C|={ps.log2_C:.1f} (need >256)   '
          f'beta={ps.beta:.3e} (log2={log2(ps.beta):.2f})')

    b, c, q, p = msis_bits(ps)
    print(f'  MSIS : classical={c}  quantum={q}  plausible={p}   BKZ-b={b}')

    if run_mlwe:
        pk  = mlwe_bits(ps.N, *ps.mlwe_pk,  ps.q)
        sim = mlwe_bits(ps.N, *ps.mlwe_sim, ps.q)
        print(f"  MLWE-pk  : classical={pk['classical']}  quantum={pk['quantum']}  "
              f"plausible={pk['plausible']}   [BINDING CONSTRAINT]")
        print(f"  MLWE-sim : classical={sim['classical']}  quantum={sim['quantum']}  "
              f"plausible={sim['plausible']}")
    elif fast_mlwe:
        pk  = mlwe_primal(ps.N, *ps.mlwe_pk,  ps.q, 'classical')
        sim = mlwe_primal(ps.N, *ps.mlwe_sim, ps.q, 'classical')
        print(f'  MLWE-pk  (primal,classical) = {pk}   [BINDING CONSTRAINT]')
        print(f'  MLWE-sim (primal,classical) = {sim}')

    print(f'  PK={ps.pk_kB:.2f} kB   Sig={ps.sig_kB:.2f} kB '
          f'(resp {ps.resp_kB:.2f} + anchor {ps.anchor_kind} {ps.anchor_bits/8/1024:.3f})   '
          f'PK+Sig={ps.tot_kB:.2f} kB')
    if appxC:
        print(f'  [Appendix-C split-error variant]  Sig={ps.sig_kB_appxC:.2f} kB   '
              f'PK+Sig={ps.tot_kB_appxC:.2f} kB   '
              f'(saves {ps.tot_kB - ps.tot_kB_appxC:.2f} kB)')
    print()


# ----------------------------------------------------------------------------
# parameter search : 
# ----------------------------------------------------------------------------
def search(n, target, model='classical',
           logq_range=range(28, 58), k_range=range(4, 11), l_range=range(4, 16),
           dlp_range=(1, 2), eta_range=(1, 2, 3), N=256, kappa=60, gamma=1.1, t=13.5,
           mlwe_gate=True, max_mlwe=60):
    idx = {'classical': 1, 'quantum': 2, 'plausible': 3}[model]

    pool = []
    for eta in eta_range:
        for logq in logq_range:
            for k in k_range:
                for l in l_range:
                    for dlp in dlp_range:
                        ps = ParamSet(N, k, l, l + dlp, logq, eta, n,
                                      kappa=kappa, gamma=gamma, t=t)
                        if ps.valid():
                            pool.append((ps.tot_kB, ps))
    pool.sort(key=lambda x: x[0])

    print('=' * 78)
    print(f'  Scheme 2 parameter search   n={n}   target={target}-bit ({model})')
    print(f'  gates: MSIS>={target} AND MLWE-pk>={target} AND MLWE-sim>={target}'
          f'{"" if mlwe_gate else "   (MLWE GATE OFF -- results may be INSECURE)"}')
    print('=' * 78)

    checked = 0
    for tot, ps in pool:
        mb = msis_bits(ps)[idx]
        if mb < target:
            continue
        if not mlwe_gate:
            print(f'  {tot:8.2f} kB | k={ps.k} l={ps.l} l\'={ps.lp} logq={ps.logq} '
                  f'eta={ps.eta} | MSIS={mb}  (MLWE NOT CHECKED -- INSECURE UNTIL VERIFIED)')
            report(ps, f'Scheme 2 n={n} (MSIS-only -- VERIFY MLWE!)', fast_mlwe=True)
            return ps
        checked += 1
        if checked > max_mlwe:
            print(f'  stopped after {max_mlwe} MLWE checks without a full pass.\n')
            return None
        pk  = mlwe_primal(ps.N, *ps.mlwe_pk,  ps.q, model)
        sim = mlwe_primal(ps.N, *ps.mlwe_sim, ps.q, model)
        ok  = pk >= target and sim >= target
        print(f'  [{checked:2d}] {tot:8.2f} kB | k={ps.k} l={ps.l} l\'={ps.lp} '
              f'logq={ps.logq} eta={ps.eta} | MSIS={mb}  MLWE-pk={pk}  MLWE-sim={sim}  '
              f'{"<== SECURE OPTIMUM" if ok else "(fails MLWE)"}')
        if ok:
            print()
            report(ps, f'Scheme 2 n={n} (secure minimum, primal MLWE)', fast_mlwe=False)
            print(f'  MLWE (primal, {model}): pk={pk}  sim={sim}   [both >= {target}]')
            print('  -> confirm with full min(primal,dual) all models:  '
                  f'python3 Scheme2_estimator.py --mlwe\n')
            return ps
    print('\n  no secure candidate in the given ranges; widen ranges.\n')
    return None


# ----------------------------------------------------------------------------
# fixed reference sets 
# ----------------------------------------------------------------------------
N = 256
FIXED = [
    ('Scheme 2 n=32',   dict(N=N, k=6, l=7, lp=8,  logq=37, eta=1, n=32)),
    ('Scheme 2 n=1024', dict(N=N, k=7, l=9, lp=10, logq=48, eta=1, n=1024)),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mlwe',         action='store_true', help='full min(primal,dual) all models (slow)')
    ap.add_argument('--appxC',        action='store_true', help='also show Appendix-C split-error sizes')
    ap.add_argument('--search',       action='store_true')
    ap.add_argument('--no-mlwe-gate', action='store_true', help='(search) disable MLWE gate -- INSECURE')
    ap.add_argument('--n',      type=int, default=32)
    ap.add_argument('--target', type=int, default=128)
    ap.add_argument('--model',  choices=list(SVP), default='classical')
    args = ap.parse_args()

    if args.search:
        search(args.n, args.target, model=args.model, mlwe_gate=not args.no_mlwe_gate)
        return

    print('=' * 78)
    print('  Scheme 2 security estimator  ')
    print('  MLWE:', 'FULL min(primal,dual) all models' if args.mlwe
          else 'fast primal-classical gate (use --mlwe for full)')
    print('=' * 78, '\n')
    for name, kw in FIXED:
        report(ParamSet(**kw), name, run_mlwe=args.mlwe, fast_mlwe=not args.mlwe,
               appxC=args.appxC)


if __name__ == '__main__':
    main()
