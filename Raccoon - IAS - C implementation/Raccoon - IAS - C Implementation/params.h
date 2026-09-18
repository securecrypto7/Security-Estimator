#ifndef IAS_PARAMS_H
#define IAS_PARAMS_H

#include <stdint.h>

/* ---- selectable parameter identifiers ---- */
#define IAS128_Q60 1
#define IAS128_Q30 2
#define IAS192_Q60 3
#define IAS192_Q30 4
#define IAS256_Q60 5
#define IAS256_Q30 6

#ifndef IAS_PARAM
#define IAS_PARAM IAS128_Q60
#endif

/* Ring degree  */
#define IAS_N 512

#define IAS_NSIGNERS 1024

#define IAS_SEEDBYTES 32
/* Symmetric security parameter (bytes) for hashes / challenge seeds. */
#define IAS_CRHBYTES 64

#if   IAS_PARAM == IAS128_Q60
  #define IAS_NAME   "IAS-128 Q2^60"
  #define IAS_KAPPA  128
  #define IAS_LOGQ   64
  #define IAS_K      6
  #define IAS_L      6
  #define IAS_NU_T   28
  #define IAS_NU_W   29
  #define IAS_NU_WP  20
  #define IAS_U_S    4        /* sigma_s = 2^u_s */
  #define IAS_U_R    36       /* sigma_r = 2^u_r  */
  #define IAS_OMEGA  19       
  #define IAS_REP    14
#elif IAS_PARAM == IAS128_Q30
  #define IAS_NAME   "IAS-128 Q2^30"
  #define IAS_KAPPA  128
  #define IAS_LOGQ   56
  #define IAS_K      5
  #define IAS_L      5
  #define IAS_NU_T   34
  #define IAS_NU_W   35
  #define IAS_NU_WP  26
  #define IAS_U_S    4
  #define IAS_U_R    23
  #define IAS_OMEGA  19
  #define IAS_REP    14
#elif IAS_PARAM == IAS192_Q60
  #define IAS_NAME   "IAS-192 Q2^60"
  #define IAS_KAPPA  192
  #define IAS_LOGQ   69
  #define IAS_K      8
  #define IAS_L      8
  #define IAS_NU_T   44
  #define IAS_NU_W   45
  #define IAS_NU_WP  36
  #define IAS_U_S    4
  #define IAS_U_R    39
  #define IAS_OMEGA  31
  #define IAS_REP    21
#elif IAS_PARAM == IAS192_Q30
  #define IAS_NAME   "IAS-192 Q2^30"
  #define IAS_KAPPA  192
  #define IAS_LOGQ   59
  #define IAS_K      7
  #define IAS_L      7
  #define IAS_NU_T   35
  #define IAS_NU_W   36
  #define IAS_NU_WP  27
  #define IAS_U_S    4
  #define IAS_U_R    24       
  #define IAS_U_R_HALF 1      
  #define IAS_OMEGA  31
  #define IAS_REP    21
#elif IAS_PARAM == IAS256_Q60
  #define IAS_NAME   "IAS-256 Q2^60"
  #define IAS_KAPPA  256
  #define IAS_LOGQ   72
  #define IAS_K      11
  #define IAS_L      10
  #define IAS_NU_T   42
  #define IAS_NU_W   48
  #define IAS_NU_WP  39
  #define IAS_U_S    4
  #define IAS_U_R    42
  #define IAS_OMEGA  44
  #define IAS_REP    27
#elif IAS_PARAM == IAS256_Q30
  #define IAS_NAME   "IAS-256 Q2^30"
  #define IAS_KAPPA  256
  #define IAS_LOGQ   63
  #define IAS_K      9
  #define IAS_L      9
  #define IAS_NU_T   37
  #define IAS_NU_W   38
  #define IAS_NU_WP  28
  #define IAS_U_S    4
  #define IAS_U_R    26
  #define IAS_OMEGA  44
  #define IAS_REP    27
#else
  #error "Unknown IAS_PARAM"
#endif

#ifndef IAS_U_R_HALF
#define IAS_U_R_HALF 0
#endif

/* sigma_s and sigma_r as doubles */
#include <math.h>
static inline double ias_sigma_s(void) { return ldexp(1.0, IAS_U_S); }
static inline double ias_sigma_r(void) {
    return IAS_U_R_HALF ? ldexp(1.0, IAS_U_R) * 0.7071067811865476  
                        : ldexp(1.0, IAS_U_R);
}

#endif 
