
import sys, os
from math import log, log2, sqrt, ceil, floor, pi

HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
for cand in (os.path.join(HERE, "lattice-estimator"), HERE, os.getcwd()):
    if os.path.isdir(os.path.join(cand, "estimator")):
        sys.path.insert(0, cand); break

from estimator import *
from sage.all import oo, log as slog, ZZ, next_prime, is_prime

# ============================ CONFIG =========================================
D_RING   = 512
 
USE_FULL = False                   
MARGIN   = 0.0     
SMOOTH_EPS = 2.0**-40   

FORGERY_BOUND = "paper"

AGG_POW = 3

HINT_MODE = "tight"

   
# =============================================================================


# ------------------------------- helpers -------------------------------------
def eta_eps(dim, eps=SMOOTH_EPS):
    return sqrt(log(2.0 * dim * (1.0 + 1.0 / eps)) / pi)


def ntt_prime(logq, n=D_RING):
    step = 2 * n
    q = (int(ZZ(2)**logq) // step) * step + 1
    if q < int(ZZ(2)**logq):
        q += step
    while not is_prime(q):
        q += step
    return q


def hint_coef(omega, n):
    if HINT_MODE == "s1":
        return 2.0 * omega * omega           
    if HINT_MODE == "hermine_safe":
        return (8.0*omega*omega + 8.0*omega + 4.0) * n  
    return 4.0 * omega * omega               


def sigma_eff(sigma_s, sigma_r, Qsign, omega, n, kappa):
    corr = 1.0 + n * (1.0 / sqrt(Qsign)) * (kappa + 1 + 2 * log(n))
    B = hint_coef(omega, n) * Qsign * corr
    return 1.0 / sqrt(2.0 * (1.0 / sigma_s**2 + B / sigma_r**2))


def scheme_bounds(P, q):
    n, k, l = D_RING, P["k"], P["l"]
    ss, sr = 2.0**P["u_s"], 2.0**P["u_r"]
    omega, rep = P["omega"], P["rep"]
    nu_t, nu_w, nu_wp = P["nu_t"], P["nu_w"], P["nu_wp"]
    ALPHA = max(2.0, 1.0 + sqrt(4.0 * P["kappa"] * log(2.0) / n))

    rt_full = sqrt(n * (k + l))
    rt_k    = sqrt(n * k)
    q_perp  = abs(q - 2.0**nu_w * round(q / 2.0**nu_w))

    
    B2ind = (2.0*omega*ALPHA*ss*sqrt(n * (k + l))
             + ALPHA*sr*sqrt(n*(k+l)*rep)
             + rep*2.0**(nu_wp-1)*rt_k)

    B2 = (P["n_signers"]*B2ind
          + (3.0*2.0**(nu_w-1)+ q_perp
          + omega*2.0**(nu_t-1))*rt_k)

    if FORGERY_BOUND == "handwritten":
        beta = (8.0*omega*B2 + 2.0**(nu_w+2)*omega*rt_k + 2.0**(nu_t+3)*omega**3)

    else: 
        beta =(8.0*omega*B2
                + (2.0**(nu_w+2)*omega + 2.0**(nu_t+2)*omega**AGG_POW)*rt_k
                + 16.0*omega**AGG_POW*ss*rt_full)
 

    return B2ind, B2, beta


def correctness_gate(P, q):
    n, k, l = D_RING, P["k"], P["l"]
    ss, sr = 2.0**P["u_s"], 2.0**P["u_r"]
    logq = log2(q)
    for nm in ("nu_t", "nu_w", "nu_wp"):
        if not (0 < P[nm] < logq):
            return False, f"{nm}={P[nm]} not in (0, logq={logq:.1f})"
    eta_s = eta_eps((l) * n)          
    eta_r = eta_eps((k + l) * n)      
    if ss < eta_s:
        return False, f"sigma_s=2^{P['u_s']} < eta_eps={eta_s:.2f}"
    if sr < eta_r:
        return False, f"sigma_r=2^{P['u_r']} < eta_eps={eta_r:.2f}"
    return True, "ok"


SEED_BYTES = 32   

def sizes(P, q):
    n, k, l = D_RING, P["k"], P["l"]
    logq = ceil(log2(q))
    _, B2, _ = scheme_bounds(P, q)
    vk = SEED_BYTES + k * n * (logq - P["nu_t"]) / 8.0
    z_bits = ceil(log2(2.0 * B2 +1))
    h_bits = ceil(log2(q))-P["nu_w"]
    c_bits = ceil(P["omega"] * (ceil(log2(n)) + 1))
    sig = (l * n * z_bits + k * n * h_bits + c_bits) / 8.0
    return vk, sig


def _bits(cost):
    best = oo
    for _, c in cost.items():
        try:
            r = c["rop"]
        except Exception:
            continue
        if r < best:
            best = r
    return float(slog(best, 2)) if best not in (oo, 0) else float("inf")


def estimate_set(P, verbose=True, final=False):
    n, k, l = D_RING, P["k"], P["l"]
    q = ntt_prime(P["logq"], n) if final else int(next_prime(ZZ(2)**P["logq"]))
    vk_size, sig_size = sizes(P, q)
    ss = 2.0**P["u_s"]
    se = sigma_eff(ss, 2.0**P["u_r"], P["Qsign"], P["omega"], n, P["kappa"])
    _, B2, beta = scheme_bounds(P, q)

    lwe = LWE.Parameters(n=l*n, q=q, Xs=ND.DiscreteGaussian(se),
                         Xe=ND.DiscreteGaussian(se), m=k*n)
    sis = SIS.Parameters(n=k*n, q=q, length_bound=float(beta), m=(k+l)*n, norm=2)
    est = (LWE.estimate, SIS.estimate) if USE_FULL else  (LWE.estimate.rough, SIS.estimate.rough)

    try:    vk = _bits(est[0](lwe))
    except Exception: vk = float("nan")
    if beta >= (q - 1) / 2.0:
        sig = 0.0
    else:
        try:    sig = _bits(est[1](sis))
        except Exception: sig = float("nan")

    if verbose:
        print(f"   (k,l)=({k},{l}) logq={P['logq']} se=2^{log2(se):.2f} "
              f"beta=2^{log2(beta):.2f}  ->  vk={vk:.1f}  sig={sig:.1f}")
        print(f"   VK_size={vk_size:.1f} bytes  SIG_size={sig_size:.1f} bytes")
    return vk, sig, q, beta


def _nu_for(base, logq):
    return dict(nu_t=base["nu_t"], nu_w=base["nu_w"], nu_wp=base["nu_wp"])


BASES = [
    dict(name="IAS-128 Q2^60", kappa=128, logq=64, k=6, l=6, nu_t=28, nu_w=29,
         nu_wp=20, u_s=4, u_r=36.0, omega=19, rep=14, Qsign=2.0**60, n_signers=1024),
    dict(name="IAS-128 Q2^30", kappa=128, logq=56, k=5, l=5, nu_t=34, nu_w=35,
         nu_wp=26, u_s=4, u_r=23, omega=19, rep=14, Qsign=2.0**30, n_signers=1024),
    dict(name="IAS-192 Q2^60", kappa=192, logq=69, k=8, l=8, nu_t=44, nu_w=45,
         nu_wp=36, u_s=4, u_r=39, omega=31, rep=21, Qsign=2.0**60, n_signers=1024),
    dict(name="IAS-192 Q2^30", kappa=192, logq=59, k=7, l=7, nu_t=35, nu_w=36,
         nu_wp=27, u_s=4, u_r=23.5, omega=31, rep=21, Qsign=2.0**30, n_signers=1024),
    dict(name="IAS-256 Q2^60", kappa=256, logq=72, k=11, l=10, nu_t=42, nu_w=48,
         nu_wp=39, u_s=4, u_r=42.0, omega=44, rep=27, Qsign=2.0**60, n_signers=1024),
    dict(name="IAS-256 Q2^30", kappa=256, logq=63, k=9, l=9, nu_t=37, nu_w=38,
         nu_wp=28, u_s=4, u_r=26.0, omega=44, rep=27, Qsign=2.0**30, n_signers=1024),
]

if __name__ == "__main__":
    print("="*80)
    print("IAS Security Estimator")
    print("="*80)

    for P in BASES:
        print(f"\n[{P['name']}]")
        estimate_set(P)

