#ifndef IAS_SERIALIZE_H
#define IAS_SERIALIZE_H

#include <stddef.h>
#include "ias.h"

/* sizes (bytes) */
size_t ias_vk_bytes  (const ias_bounds *b);   /* == estimator VK  */
size_t ias_sig_z_bytes(const ias_bounds *b);  /* == estimator SIG */
size_t ias_sig_c_bytes(void);
size_t ias_sig_h_bytes(void);
size_t ias_sig_full_bytes(const ias_bounds *b);

/* pack / unpack */
void ias_pack_vk(uint8_t *out, const ias_pp *pp, const ias_pk *pk,
                 const ias_bounds *b);
void ias_pack_z (uint8_t *out, const ias_sig *sig, const ias_bounds *b);
int  ias_unpack_z(spoly z[IAS_L], const uint8_t *in, const ias_bounds *b);
void ias_pack_c (uint8_t *out, const poly *c);
void ias_pack_h (uint8_t *out, const ias_sig *sig);

#endif
