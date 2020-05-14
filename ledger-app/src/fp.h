#ifndef FP_H
#define FP_H

#include <stdint.h>

//////////////////////////////////////////////////
// Residues mod p
//////////////////////////////////////////////////

typedef struct {
    uint32_t x[12];
} Fp;

// Low-level operations implemented in assembler

int  fp_is_zero(const Fp *);
int  fp_eq(const Fp *, const Fp *);

void fp_cpy (Fp *, const Fp *);
void fp_neg (Fp *, const Fp *);
void fp_sum (Fp *, const Fp *, const Fp *);
void fp_diff(Fp *, const Fp *, const Fp *);

void fp_cset(Fp *, const Fp *, uint8_t);

void fp_to_bytes	(uint8_t *, const Fp *);
int fp_from_bytes	(Fp *, const uint8_t *);

extern void mul378(uint32_t *, const uint32_t *, const uint32_t *);

inline void fp_prod(uint32_t *x, const Fp *y, const Fp *z) { mul378(x, y->x, z->x); }

// Higher-level operations implemented in C

void fp_redc(uint32_t *, uint32_t *);
void fp_inv (Fp *, const Fp *);
void fp_mulred(Fp *, const Fp *, const Fp *);

// Convert between:
// - (external) least non-negative residue in big-endian array of bytes
// - (internal) Montgomery representation in little-endian array of words

void fp_conv_to_bytes	(uint8_t *, const Fp *);
void fp_conv_from_bytes	(Fp *, const uint8_t *);

//////////////////////////////
// Constants
//////////////////////////////

#define FP_ZERO {{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }}
#define FP_ONE	{{ \
    0xffffff68, \
    0x02cdffff, \
    0x7fffffb1, \
    0x51409f83, \
    0x8a7d3ff2, \
    0x9f7db3a9, \
    0x6e7c6305, \
    0x7b4e97b7, \
    0x803c84e8, \
    0x4cf495bf, \
    0xe2fdf49a, \
    0x008d6661  \
}}
#define FP_R	{{ \
    0x9400cd22, \
    0xb786686c, \
    0xb00431b1, \
    0x0329fcaa, \
    0x62d6b46d, \
    0x22a5f111, \
    0x827dc3ac, \
    0xbfdf7d03, \
    0x41790bf9, \
    0x837e92f0, \
    0x1e914b88, \
    0x006dfccb  \
}}
#define FP_P	{{ \
    0x00000001, \
    0x8508C000, \
    0x30000000, \
    0x170B5D44, \
    0xBA094800, \
    0x1EF3622F, \
    0x00F5138F, \
    0x1A22D9F3, \
    0x6CA1493B, \
    0xC63B05C0, \
    0x17C510EA, \
    0x01AE3A46  \
}}
#define FP_PM2	{{ \
    0xFFFFFFFF, \
    0x8508BFFF, \
    0x30000000, \
    0x170B5D44, \
    0xBA094800, \
    0x1EF3622F, \
    0x00F5138F, \
    0x1A22D9F3, \
    0x6CA1493B, \
    0xC63B05C0, \
    0x17C510EA, \
    0x01AE3A46  \
}}

extern const Fp fp_zero;
extern const Fp fp_one;	// R = 2^384 mod p = Montgomery representation of 1
extern const Fp fp_r;	// R^2 = 2^768 mod p = Montgomery representation of R
extern const Fp fp_pm2;	// pm2 = p-2

#endif // FP_H
