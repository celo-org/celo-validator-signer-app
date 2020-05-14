#ifndef FQ_H
#define FQ_H

#include <stdint.h>

typedef struct {
    uint32_t x[8];
} Fq;

////////////////////////////////////////
// Operations on residues mod q
////////////////////////////////////////

// Low-level operations implemented in assembler

int  fq_is_zero(const Fq *);
int  fq_eq(const Fq *, const Fq *);

void fq_cpy (Fq *, const Fq *);
void fq_neg (Fq *, const Fq *);
void fq_sum (Fq *, const Fq *, const Fq *);
void fq_diff(Fq *, const Fq *, const Fq *);

void fq_cset(Fq *, const Fq *, uint8_t);

void fq_to_bytes	(uint8_t *, const Fq *);
void fq_from_bytes	(Fq *, const uint8_t *);

//////////////////////////////
// Constants
//////////////////////////////

#define FQ_ZERO {{ 0, 0, 0, 0, 0, 0, 0, 0 }}
#define FQ_Q	{{ \
    0x00000001, \
    0x0A118000, \
    0xD0000001, \
    0x59AA76FE, \
    0x5C37B001, \
    0x60B44D1E, \
    0x9A2CA556, \
    0x12AB655E  \
}}

#endif // FQ_H
