#ifndef G1_H
#define G1_H

#include "fp.h"
#include "fq.h"

typedef struct {
    Fp x, y;
    unsigned inf : 1;
} G1Affine;

typedef struct {
    Fp x, y, z;
} G1Jacobian;

// Constants

extern const G1Affine	g1a_identity;
extern const G1Jacobian	g1j_identity;

int g1a_eq_a	(const G1Affine *,	const G1Affine *);
int g1a_eq_j	(const G1Affine *,	const G1Jacobian *);
int g1j_eq_a	(const G1Jacobian *,	const G1Affine *);
int g1j_eq_j	(const G1Jacobian *,	const G1Jacobian *);

int g1a_valid	(const G1Affine *);
int g1j_valid	(const G1Jacobian *);

void g1a_sanitize	(G1Jacobian *,	const G1Affine *);

void g1a_cpy	(G1Affine *,	const G1Affine *);
void g1j_cpy	(G1Jacobian *,	const G1Jacobian *);
void g1a_from_j	(G1Affine *,	const G1Jacobian *);
void g1j_from_a	(G1Jacobian *,	const G1Affine *);

void g1a_cset	(G1Affine *,	const G1Affine *,	int);
void g1j_cset	(G1Jacobian *,	const G1Jacobian *,	int);

void g1j_add_a	(G1Jacobian *,	const G1Affine *);
void g1j_add_j	(G1Jacobian *,	const G1Jacobian *);
void g1j_dbl	(G1Jacobian *,	const G1Jacobian *);
void g1a_dbl	(G1Jacobian *,	const G1Affine *);

void g1a_mul	(G1Jacobian *,	const G1Affine *,	const Fq *);
void g1j_mul	(G1Jacobian *,	const G1Jacobian *,	const Fq *);

void g1a_from_bytes	(G1Affine *, const uint8_t *);
void g1a_to_bytes	(uint8_t *, const G1Affine *);
void g1j_to_bytes	(uint8_t *, const G1Jacobian *);

#endif // G1_H
