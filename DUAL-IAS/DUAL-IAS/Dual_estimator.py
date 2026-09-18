
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
# prime generation :  
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
        return None                       # insecure: B >= q  (estimator returns 0 bits)
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


def msis_bits(ps, model='classical'):
    b = _msis_bmin(ps.q, ps.msis_w * ps.N, ps.msis_k * ps.N, ps.beta)
    if b is None:
        return 0, 0, 0, 0
    c = int(svp_classical(b)); qb = int(svp_quantum(b)); p = int(svp_plausible(b))
    return b, c, qb, p


# ----------------------------------------------------------------------------
# MLWE :
# ----------------------------------------------------------------------------
def mlwe_bits(N, dim, k, eta, q):
    s = sqrt(eta * (eta + 1) / 3.0)
    nd, mm = N * dim, N * k
    out = {}
    for tag, model in SVP.items():
        _, _, cp = _quiet(MLWE_optimize_attack, q, nd, mm, s, LWE_primal_cost, model, False)
        _, _, cd = _quiet(MLWE_optimize_attack, q, nd, mm, s, LWE_dual_cost,   model, False)
        out[tag] = min(int(floor(cp)), int(floor(cd)))
    return out                             # {'classical':..,'quantum':..,'plausible':..}


def mlwe_primal(N, dim, k, eta, q, model='classical'):
    s = sqrt(eta * (eta + 1) / 3.0)
    _, _, cp = _quiet(MLWE_optimize_attack, q, N * dim, N * k, s,
                      LWE_primal_cost, SVP[model], False)
    return int(floor(cp))


# ----------------------------------------------------------------------------
# AggMS parameter set : 
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

        self.chal = kappa

        self.T  = self.chal * self.eta  * sqrt(N * m)
        self.Tp = self.chal * self.etap * sqrt(N * mp)

        lq = log2(self.q)
        s_rej  = log2(SQRT_2PI * self.alpha * self.T)
        s_reg  = 1 + log2(N) + lq * (k / m  + 2 / (N * m))       
        sp_rej = log2(SQRT_2PI * self.alpha * self.Tp)
        sp_reg = 1 + log2(N) + lq * (k / mp + 2 / (N * mp))
        self.log2_s,  self.s_bind  = (s_rej,  'rej') if s_rej  > s_reg  else (s_reg,  'reg')
        self.log2_sp, self.sp_bind = (sp_rej, 'rej') if sp_rej > sp_reg else (sp_reg, 'reg')
        self.s, self.sp = 2 ** self.log2_s, 2 ** self.log2_sp

        
        self.M  = exp(t / self.alpha + 1 / (2 * self.alpha ** 2))
        self.Mn = self.M ** n

        self.Bn   = gamma * (self.s  / SQRT_2PI) * sqrt(N * l ) * sqrt(n)
        self.Bpn  = gamma * (self.sp / SQRT_2PI) * sqrt(N * lp) * sqrt(n)
        self.Bppn = gamma * (sqrt(self.s ** 2 + self.sp ** 2) / SQRT_2PI) * sqrt(N * k) * sqrt(n)

        core = self.Bn ** 2 + self.Bpn ** 2 + self.Bppn ** 2
        self.beta     = sqrt(core)                       
        self.beta_thm = 2 * sqrt(kappa + core)          

        self.log2_C = log2(comb(N, kappa)) + kappa            

        
        self.msis_w = 1 + l + lp + k     
        self.msis_k = k
        self.mlwe_pk  = (l,      k, self.eta )   # pk pseudorandomness
        self.mlwe_sim = (lp - 1, k, self.etap)   # simulation indistinguishability

        # ---- sizes -----------------------------------------------------------
        qbits = ceil(log2(self.q))          # bits to store one element of Z_q
        sig_z = self.s  * sqrt(n) / SQRT_2PI
        sig_r = self.sp * sqrt(n) / SQRT_2PI
        sig_e = sqrt(n * (self.s ** 2 + self.sp ** 2)) / SQRT_2PI
        resp = (N * l  * (log2(sig_z) + GAUSS_H)
              + N * lp * (log2(sig_r) + GAUSS_H)
              + N * k  * (log2(sig_e) + GAUSS_H))

        self.anchor_bits, self.anchor_kind = k * N * qbits, 'w~'

        self.pk_kB   = (k * N * qbits) / 8 / 1024
        self.resp_kB = resp / 8 / 1024
        self.sig_kB  = (resp + self.anchor_bits) / 8 / 1024
        self.tot_kB  = self.pk_kB + self.sig_kB

    def valid(self):
        return (self.q % 8 == 5) and (self.log2_C > 256) and (self.beta < self.q)


# ----------------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------------
def report(ps, name, run_mlwe=False):
    print('-' * 78)
    print(f'  {name}   [AggMS]')
    print(f"  N={ps.N} k={ps.k} l={ps.l} l'={ps.lp} logq={ps.logq} eta={ps.eta} "
          f'n={ps.n} kappa={ps.kappa}')
    print(f'  q={ps.q}  (q mod 8={ps.q % 8})   alpha={ps.alpha:.0f}  '
          f'M={ps.M:.4f}  M^n={ps.Mn:.2f}')
    print(f"  log2(s)={ps.log2_s:.2f} [{ps.s_bind}]   log2(s')={ps.log2_sp:.2f} [{ps.sp_bind}]")
    print(f"  log2(Bn)={log2(ps.Bn):.2f}  log2(B'n)={log2(ps.Bpn):.2f}  "
          f"log2(B''n)={log2(ps.Bppn):.2f}")
    print(f'  log2|C|={ps.log2_C:.1f} (need >256)')
    print(f"  beta_table  =sqrt(Bn^2+B'n^2+B''n^2)         log2={log2(ps.beta):.2f}  ")
    print(f"  beta_theorem=2*sqrt(kappa+Bn^2+B'n^2+B''n^2) log2={log2(ps.beta_thm):.2f}  ")

    b, c, q, p = msis_bits(ps)
    print(f'  MSIS : classical={c}  quantum={q}  plausible={p}   BKZ-b={b}')

    if run_mlwe:
        pk  = mlwe_bits(ps.N, *ps.mlwe_pk,  ps.q)
        sim = mlwe_bits(ps.N, *ps.mlwe_sim, ps.q)
        print(f"  MLWE-pk  : classical={pk['classical']}  quantum={pk['quantum']}  "
              f"plausible={pk['plausible']}")
        print(f"  MLWE-sim : classical={sim['classical']}  quantum={sim['quantum']}  "
              f"plausible={sim['plausible']}")

    print(f'  PK={ps.pk_kB:.2f} kB   Sig={ps.sig_kB:.2f} kB '
          f'(resp {ps.resp_kB:.2f} + anchor {ps.anchor_kind} {ps.anchor_bits/8/1024:.2f})   '
          f'PK+Sig={ps.tot_kB:.2f} kB')
    print()


# ----------------------------------------------------------------------------
# parameter search : 
# ----------------------------------------------------------------------------
def search(n, target, model='classical',
           logq_range=range(30, 56), k_range=range(4, 11), l_range=range(4, 16),
           dlp=1, eta=1, N=256, kappa=60, gamma=1.1, t=13.5,
           mlwe_verify=False, max_mlwe=8, top=12):
    idx = {'classical': 1, 'quantum': 2, 'plausible': 3}[model]  

    pool = []
    for logq in logq_range:
        for k in k_range:
            for l in l_range:
                for lp in (l, l + dlp):          
                    ps = ParamSet(N, k, l, lp, logq, eta, n,
                                  kappa=kappa, gamma=gamma, t=t)
                    if ps.valid():
                        pool.append((ps.tot_kB, ps))
    pool.sort(key=lambda x: x[0])

    cand, evals = [], 0
    for _, ps in pool:
        evals += 1
        mb = msis_bits(ps)[idx]
        if mb >= target:
            cand.append((ps.tot_kB, mb, ps))
            if len(cand) >= top:
                break
    print(f'  (evaluated MSIS on {evals}/{len(pool)} size-sorted candidates)')

    print('=' * 78)
    print(f'  AggMS parameter search   n={n}   target={target}-bit ({model})')
    print(f'  MSIS-feasible candidates: {len(cand)}   (showing smallest {min(top,len(cand))})')
    print('=' * 78)
    print(f'  {"PK+Sig":>8} {"PK":>6} {"Sig":>7} | {"k":>2} {"l":>2} {"l\'":>2} '
          f'{"logq":>4} | {"MSIS":>4} {"logs":>5} {"logs\'":>5}')
    for tot, mb, ps in cand[:top]:
        print(f'  {tot:8.2f} {ps.pk_kB:6.2f} {ps.sig_kB:7.2f} | {ps.k:2d} {ps.l:2d} '
              f'{ps.lp:2d} {ps.logq:4d} | {mb:4d} {ps.log2_s:5.1f} {ps.log2_sp:5.1f}')
    print()

    if not cand:
        print('  no MSIS-feasible candidate in the given ranges.\n')
        return None

    if not mlwe_verify:
        _, _, best = cand[0]
        print('  smallest MSIS-feasible set (MLWE NOT yet verified -- run --mlwe-verify):')
        report(best, f'AggMS n={n} (MSIS-min)')
        return best

    print('  verifying MLWE on size-sorted candidates (fast primal-only gate)...\n')
    checked = 0
    for tot, mb, ps in cand:
        if checked >= max_mlwe:
            print(f'  stopped after {max_mlwe} MLWE checks without a full pass.')
            break
        checked += 1
        pk  = mlwe_primal(ps.N, *ps.mlwe_pk,  ps.q, model)
        sim = mlwe_primal(ps.N, *ps.mlwe_sim, ps.q, model)
        ok  = pk >= target and sim >= target
        print(f'  [{checked}] k={ps.k} l={ps.l} logq={ps.logq}  PK+Sig={tot:.2f}kB  '
              f'MSIS={mb}  MLWE-pk(primal)={pk}  MLWE-sim(primal)={sim}  '
              f'{"<-- PASS" if ok else "(fails MLWE)"}')
        if ok:
            print(f'\n  WINNER: minimal PK+Sig meeting MSIS>={target} and MLWE(primal)>={target}.')
            report(ps, f'AggMS n={n} (minimal)', run_mlwe=False)
            print(f'  MLWE (primal, fast): pk={pk}  sim={sim}   [both >= {target}]')
            print(f'  -> confirm with the exact primal+dual estimate by pasting these params '
                  f'into AGGMS and running:  python3 {__file__.split("/")[-1]} --mlwe\n')
            return ps
    print('\n  no candidate passed within the window; widen ranges / raise --max_mlwe.\n')
    return None


# ----------------------------------------------------------------------------
# fixed reference sets
# ----------------------------------------------------------------------------
N = 256
AGGMS = [
    ('AggMS n=32',   dict(N=N, k=5, l=7, lp=8, logq=32, eta=1, n=32)),
    ('AggMS n=1024', dict(N=N, k=6, l=8, lp=9, logq=41, eta=1, n=1024)),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mlwe',        action='store_true', help='score MLWE for fixed sets (slow)')
    ap.add_argument('--search',      action='store_true', help='run the parameter search')
    ap.add_argument('--mlwe-verify', action='store_true', help='verify MLWE inside the search (slow)')
    ap.add_argument('--n',      type=int,   default=32)
    ap.add_argument('--target', type=int,   default=128)
    ap.add_argument('--model',  choices=list(SVP), default='classical')
    args = ap.parse_args()

    if args.search:
        search(args.n, args.target, model=args.model, mlwe_verify=args.mlwe_verify)
        return

    print('=' * 78)
    print('  AggMS security estimator   ')
    print('  MLWE:', 'ENABLED' if args.mlwe else 'skipped (use --mlwe)')
    print('=' * 78, '\n')
    for name, kw in AGGMS:
        report(ParamSet(**kw), name, run_mlwe=args.mlwe)


if __name__ == '__main__':
    main()
