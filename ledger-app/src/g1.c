#include <assert.h>
#include <stdio.h>

#include "g1.h"

//////////////////////////////
// Constants
//////////////////////////////

const G1Affine g1a_identity =	// (0, 1, inf)
{
    FP_ZERO, FP_ONE, 1
};

const G1Jacobian g1j_identity =	// (0, 1, 0)
{
    FP_ZERO, FP_ONE, FP_ZERO
};

//////////////////////////////
// Conversion
//////////////////////////////

void g1a_from_bytes(G1Affine *a, const uint8_t *p)
{
    fp_conv_from_bytes(&a->x, p);
    fp_conv_from_bytes(&a->y, p+48);

    int x_0 = fp_eq(&a->x, &fp_zero);
    int y_1 = fp_eq(&a->y, &fp_one);

    a->inf = x_0 & y_1;
}

void g1a_to_bytes(uint8_t *p, const G1Affine *a)
{
    fp_conv_to_bytes(p,    &a->x);
    fp_conv_to_bytes(p+48, &a->y);
}

void g1j_to_bytes(uint8_t *p, const G1Jacobian *j)
{
    G1Affine a;

    g1a_from_j(&a, j);

    fp_conv_to_bytes(p,    &a.x);
    fp_conv_to_bytes(p+48, &a.y);
}

void g1a_from_j(G1Affine *a, const G1Jacobian *j)
{
    // x = j.x / j.z^2, y = j.y / j.z^3

    Fp r1, r2, r3;

    const Fp *jx = &j->x;
    const Fp *jy = &j->y;
    const Fp *jz = &j->z;
          Fp *ax = &a->x;
          Fp *ay = &a->y;

    fp_inv(&r1, jz);		// a.y = j.z^(p-2)
    fp_mulred(&r2, &r1, &r1);	// a.x  = a.y * a.y
    fp_mulred(&r3, &r2, &r1);	// a.y *= a.x
    fp_mulred( ax, &r2,  jx);	// a.x *= j.x
    fp_mulred( ay, &r3,  jy);	// a.y *= j.y

    a->inf = fp_is_zero(jz);	// inf = j.z == 0

    fp_cset(ax, &fp_zero, a->inf);	// a.x = inf ? 0 : a.x
    fp_cset(ay, &fp_one, a->inf);	// a.y = inf ? 1 : a.y
}

void g1j_from_a(G1Jacobian *j, const G1Affine *a)
{
    fp_cpy(&j->x, &a->x);
    fp_cpy(&j->y, &a->y);
    fp_cpy(&j->z, &fp_one);

    fp_cset(&j->x, &fp_zero, a->inf);
    fp_cset(&j->y, &fp_one,  a->inf);
    fp_cset(&j->z, &fp_zero, a->inf);
}

//////////////////////////////
// Comparison
//////////////////////////////

int g1a_eq_a(const G1Affine *a, const G1Affine *b)
{
    const Fp *x1 = &a->x;
    const Fp *y1 = &a->y;

    const Fp *x2 = &b->x;
    const Fp *y2 = &b->y;

    int eq_x = fp_eq(x1, x2);
    int eq_y = fp_eq(y1, y2);

    int eq_inf = a->inf ^ b->inf ^ 1;

    return eq_x & eq_y & eq_inf;
}

int g1j_eq_j(const G1Jacobian *j, const G1Jacobian *k)
{
    int eq_x, eq_y;
    Fp t, u, v, w;

    const Fp *x1 = &j->x;
    const Fp *y1 = &j->y;
    const Fp *z1 = &j->z;

    const Fp *x2 = &k->x;
    const Fp *y2 = &k->y;
    const Fp *z2 = &k->z;

    // x1/z1^2 == x2/z2^2 => x1*z2^2 == x2*z1^2

    fp_mulred(&t, z2, z2);	// t = z2^2
    fp_mulred(&u, z1, z1);	// u = z1^2

    fp_mulred(&v, &t, x1);	// v = x1 * z2^2
    fp_mulred(&w, &u, x2);	// w = x2 * z1^2

    eq_x = fp_eq(&v, &w);

    // y1/z1^3 == y2/z2^3 => y1*z2^3 == y2*z1^3

    fp_mulred(&t, &t, z2);	// t = z2^3
    fp_mulred(&u, &u, z1);	// u = z1^3

    fp_mulred(&v, &t, y1);	// v = y1 * z2^3
    fp_mulred(&w, &u, y2);	// w = y2 * z1^3

    eq_y = fp_eq(&v, &w);

    return eq_x & eq_y;
}

int g1j_eq_a(const G1Jacobian *j, const G1Affine *a)
{
    int eq_x, eq_y, eq_inf;
    Fp t, u;

    const Fp *x1 = &j->x;
    const Fp *y1 = &j->y;
    const Fp *z1 = &j->z;

    const Fp *x2 = &a->x;
    const Fp *y2 = &a->y;

    // x1/z1^2 == x2 => x1 == x2*z1^2

    fp_mulred(&t, z1, z1);	// t = z1^2
    fp_mulred(&u, &t, x2);	// u = x2 * z1^2

    eq_x = fp_eq(x1, &u);

    // y1/z1^3 == y2 => y1 == y2*z1^3

    fp_mulred(&t, &t, z1);	// t = z1^3
    fp_mulred(&u, &t, y2);	// u = y2 * z1^3

    eq_y = fp_eq(y1, &u);

    eq_inf = a->inf && fp_is_zero(z1);

    return eq_inf | (eq_x & eq_y);
}

int g1a_eq_j(const G1Affine *a, const G1Jacobian *j)
{
    return g1j_eq_a(j, a);
}

//////////////////////////////
// Basic operations
//////////////////////////////

void g1a_cpy(G1Affine *a, const G1Affine *b)
{
    fp_cpy(&a->x, &b->x);
    fp_cpy(&a->y, &b->y);
    a->inf = b->inf;
}

void g1j_cpy(G1Jacobian *j, const G1Jacobian *k)
{
    fp_cpy(&j->x, &k->x);
    fp_cpy(&j->y, &k->y);
    fp_cpy(&j->z, &k->z);
}

void g1a_cset(G1Affine *a, const G1Affine *b, int c)
{
    fp_cset(&a->x, &b->x, c);
    fp_cset(&a->y, &b->y, c);
    a->inf = (b->inf & c) | (a->inf & ~c);
}

void g1j_cset(G1Jacobian *j, const G1Jacobian *k, int c)
{
    fp_cset(&j->x, &k->x, c);
    fp_cset(&j->y, &k->y, c);
    fp_cset(&j->z, &k->z, c);
}

//////////////////////////////
// Sanitize point
//////////////////////////////

// Multiplies a point by sqrt(3h), the square root of 3 times the cofactor,
// to ensure that the resulting point lies in the correct group.

// sqrt(3h) = 0x8508C00000000000 = 0b1000010100001000110000000000000000000000000000000000000000000000

void g1a_sanitize(G1Jacobian *j, const G1Affine *a)
{
    G1Jacobian t;

    g1a_dbl(j, a);

    for (int i=0; i<4; i++) g1j_dbl(j, j);

    g1j_add_a(j, a); g1j_cpy(&t, j);	// 0b100001 * a

    for (int i=0; i<7; i++) g1j_dbl(j, j);

    g1j_add_j(j, &t);

    for (int i=0; i<4; i++) g1j_dbl(j, j);

    g1j_add_a(j, a);

    g1j_dbl(j, j);

    g1j_add_a(j, a);

    for (int i=0; i<46; i++) g1j_dbl(j, j);
}

//////////////////////////////
// Point arithmetic
//////////////////////////////

void g1a_neg(G1Affine *a, const G1Affine *b)
{
    fp_cpy (&a->x, &b->x);
    fp_neg (&a->y, &b->y);
    a->inf = b->inf;

    fp_cset(&a->y, &b->y, b->inf);
}

void g1j_neg(G1Jacobian *j, const G1Jacobian *k)
{
    fp_cpy(&j->x, &k->x);
    fp_neg(&j->y, &k->y);
    fp_cpy(&j->z, &k->z);

    fp_cset(&j->y, &k->y, fp_is_zero(&k->z));
}

// pub fn double(&self) -> G1Projective
void g1a_dbl(G1Jacobian *j, const G1Affine *self)
{
    // Temporary variables

    Fp t0, t1;
    Fp x3, y3, z3;

    // Placement of intermediate values in temporary variables

    Fp *a = &z3;
    Fp *b = &y3;
    Fp *c = &t1;
    Fp *d = &y3;
    Fp *e = &t0;
    Fp *f = &x3;

    // let a = self.x.square();
    fp_mulred(a, &self->x, &self->x);

    // let b = self.y.square();
    fp_mulred(b, &self->y, &self->y);

    // let c = b.square();
    fp_mulred(c, b, b);

    // let d = self.x + b;
    fp_sum(d, &self->x, b);

    // let d = d.square();
    fp_mulred(d, d, d);

    // let d = d - a - c;
    fp_diff(d, d, a);
    fp_diff(d, d, c);

    // let d = d + d;
    fp_sum(d, d, d);

    // let e = a + a + a;
    fp_sum(e, a, a);
    fp_sum(e, e, a);

    // let f = e.square();
    fp_mulred(f, e, e);

    // let z3 = self.y + self.y;
    fp_sum(&z3, &self->y, &self->y);

    // let x3 = f - (d + d);
    fp_diff(&x3, f, d);
    fp_diff(&x3, &x3, d);

    // let c = c + c;
    // let c = c + c;
    // let c = c + c;
    fp_sum(c, c, c);
    fp_sum(c, c, c);
    fp_sum(c, c, c);

    // let y3 = e * (d - x3) - c;
    fp_diff  (&y3, d, &x3);
    fp_mulred(&y3, e, &y3);
    fp_diff  (&y3, &y3, c);

    fp_cpy(&j->x, &x3);
    fp_cpy(&j->y, &y3);
    fp_cpy(&j->z, &z3);
}

// pub fn double(&self) -> G1Projective
void g1j_dbl(G1Jacobian *k, const G1Jacobian *j)
{
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

    // Temporary variables

    Fp t0, t1;
    Fp x3, y3, z3;

    // Placement of intermediate values in temporary variables

    Fp *a = &z3;
    Fp *b = &y3;
    Fp *c = &t1;
    Fp *d = &y3;
    Fp *e = &t0;
    Fp *f = &x3;

    // let a = self.x.square();
    fp_mulred(a, &j->x, &j->x);

    // let b = self.y.square();
    fp_mulred(b, &j->y, &j->y);

    // let c = b.square();
    fp_mulred(c, b, b);

    // let d = self.x + b;
    fp_sum(d, &j->x, b);

    // let d = d.square();
    fp_mulred(d, d, d);

    // let d = d - a - c;
    fp_diff(d, d, a);
    fp_diff(d, d, c);

    // let d = d + d;
    fp_sum(d, d, d);

    // let e = a + a + a;
    fp_sum(e, a, a);
    fp_sum(e, e, a);

    // let f = e.square();
    fp_mulred(f, e, e);

    // let z3 = self.z * self.y;
    fp_mulred(&z3, &j->z, &j->y);

    // let z3 = z3 + z3;
    fp_sum(&z3, &z3, &z3);

    // let x3 = f - (d + d);
    fp_diff(&x3, f, d);
    fp_diff(&x3, &x3, d);

    // let c = c + c;
    // let c = c + c;
    // let c = c + c;
    fp_sum(c, c, c);
    fp_sum(c, c, c);
    fp_sum(c, c, c);

    // let y3 = e * (d - x3) - c;
    fp_diff(&y3, d, &x3);
    fp_mulred(&y3, e, &y3);
    fp_diff(&y3, &y3, c);

    // let tmp = G1Projective { x: x3, y: y3, z: z3, };

    // G1Projective::conditional_select(&tmp, &G1Projective::identity(), self.is_identity())
    {
	const int f = fp_is_zero(&j->z);
	fp_cset(&x3, &fp_zero, f);
	fp_cset(&y3, &fp_one,  f);
	fp_cset(&z3, &fp_zero, f);
    }

    fp_cpy(&k->x, &x3);
    fp_cpy(&k->y, &y3);
    fp_cpy(&k->z, &z3);
}

// pub fn add_mixed(&self, rhs: &G1Affine) -> G1Projective
void g1j_add_a(G1Jacobian *j, const G1Affine *a)
{
    // This Jacobian point addition technique is based on the implementation in libsecp256k1,
    // which assumes that rhs has z=1. Let's address the case of zero z-coordinates generally.

    // Temporary variables

    Fp t0, t1, t2;
    Fp x3, y3, z3;

    // Placement of intermediate values in temporary variables

    Fp *z	= &t0;
    Fp *s2	= z;
    Fp *m	= s2;
    Fp *T	= m;

    Fp *u2	= &z3;
    Fp *m_alt	= u2;

    Fp *t	= &t1;
    Fp *q	= t;

    Fp *rr	= &y3;
    Fp *n	= rr;

    Fp *tt	= &t2;
    Fp *rr_alt	= tt;

    // If self is the identity, return rhs. Otherwise, return self. The other cases will be
    // predicated on neither self nor rhs being the identity.

    // let f1 = self.is_identity();
    const int f1 = fp_is_zero(&j->z);

    // let res = G1Projective::conditional_select(self, &G1Projective::from(rhs), f1);
    fp_cset(&j->x, &a->x, f1);
    fp_cset(&j->y, &a->y, f1);
    fp_cpy (&t0, &fp_one);
    fp_cset(&t0, &fp_zero, a->inf);
    fp_cset(&j->z, &t0, f1);

    // let f2 = rhs.is_identity();
    const int f2 = a->inf;

    // If neither is the identity but x1 = x2 and y1 != y2, then return the identity

    // let u1 = self.x;
    const Fp *u1 = &j->x;

    // let s1 = self.y;
    const Fp *s1 = &j->y;

    // let z = self.z.square();
    fp_mulred(z, &j->z, &j->z);

    // let u2 = rhs.x * z;
    fp_mulred(u2, &a->x, z);

    // let z = z * self.z;
    fp_mulred(z, z, &j->z);

    // let s2 = rhs.y * z;
    fp_mulred(s2, &a->y, z);

    // let f3 = u1.ct_eq(&u2) & (!s1.ct_eq(&s2));
    const int f3 = fp_eq(u1, u2) & !fp_eq(s1, s2);

    // let res = G1Projective::conditional_select(&res, &G1Projective::identity(), (!f1) & (!f2) & f3);
    {
	int f = (!f1) & (!f2) & f3;

	fp_cset(&j->x, &fp_zero, f);
	fp_cset(&j->y, &fp_one,  f);
	fp_cset(&j->z, &fp_zero, f);
    }

    // let t = u1 + u2;
    fp_sum(t, u1, u2);

    // let m = s1 + s2;
    fp_sum(m, s1, s2);

    // let rr = t.square();
    fp_mulred(rr, t, t);

    // let m_alt = -u2;
    fp_neg(m_alt, u2);

    // let tt = u1 * m_alt;
    fp_mulred(tt, u1, m_alt);

    // let rr = rr + tt;
    fp_sum(rr, rr, tt);

    // Correct for x1 != x2 but y1 = -y2, which can occur because p - 1 is divisible by 3.
    // libsecp256k1 does this by substituting in an alternative (defined) expression for lambda.

    // let degenerate = m.is_zero() & rr.is_zero();
    const int degenerate = fp_is_zero(m) & fp_is_zero(rr);

    // let rr_alt = s1 + s1;
    fp_sum(rr_alt, s1, s1);

    // let m_alt = m_alt + u1;
    fp_sum(m_alt, m_alt, u1);

    // let rr_alt = Fp::conditional_select(&rr_alt, &rr, !degenerate);
    fp_cset(rr_alt, rr, !degenerate);

    // let m_alt = Fp::conditional_select(&m_alt, &m, !degenerate);
    fp_cset(m_alt, m, !degenerate);

    // let n = m_alt.square();
    fp_mulred(n, m_alt, m_alt);

    // let q = n * t;
    fp_mulred(q, n, t);

    // let n = n.square();
    fp_mulred(n, n, n);

    // let n = Fp::conditional_select(&n, &m, degenerate);
    fp_cset(n, m, degenerate);

    // let t = rr_alt.square();
    fp_mulred(T, rr_alt, rr_alt);

    // let z3 = m_alt * self.z;
    fp_mulred(&z3, m_alt, &j->z);

    // let z3 = z3 + z3;
    fp_sum(&z3, &z3, &z3);

    // let q = -q;
    fp_neg(q, q);

    // let t = t + q;
    // let x3 = t;
    // let t = t + t;
    fp_sum(&x3, T, q);
    fp_sum(T, &x3, &x3);

    // let t = t + q;
    fp_sum(T, T, q);

    // let t = t * rr_alt;
    fp_mulred(T, T, rr_alt);

    // let t = t + n;
    fp_sum(T, T, n);

    // let y3 = -t;
    fp_neg(&y3, T);

    // let x3 = x3 + x3;
    // let x3 = x3 + x3;
    fp_sum(&x3, &x3, &x3);
    fp_sum(&x3, &x3, &x3);

    // let y3 = y3 + y3;
    // let y3 = y3 + y3;
    fp_sum(&y3, &y3, &y3);
    fp_sum(&y3, &y3, &y3);

    // let tmp = G1Projective { x: x3, y: y3, z: z3, };

    // G1Projective::conditional_select(&res, &tmp, (!f1) & (!f2) & (!f3))
    {
	int f = (!f1) & (!f2) & (!f3);

	fp_cset(&j->x, &x3, f);
	fp_cset(&j->y, &y3, f);
	fp_cset(&j->z, &z3, f);
    }
}

// pub fn add(&self, rhs: &G1Projective) -> G1Projective
void g1j_add_j(G1Jacobian *self, const G1Jacobian *rhs)
{
    // This Jacobian point addition technique is based on the implementation in libsecp256k1,
    // which assumes that rhs has z=1. Let's address the case of zero z-coordinates generally.

    G1Jacobian *res = self;

    // Temporary variables

    Fp t0, t1, t2, t3;
    Fp x3, y3, z3;

    // Placement of intermediate values in temporary variables

    Fp *z	= &x3;
    Fp *s2	= z;
    Fp *m	= s2;

    Fp *u1	= &y3;
    Fp *T	= u1;

    Fp *u2	= &z3;
    Fp *m_alt	= u2;

    Fp *t	= &t0;
    Fp *q	= t;

    Fp *rr	= &t1;
    Fp *n	= rr;

    Fp *tt	= &t2;

    Fp *rr_alt	= &t3;
    Fp *s1	= rr_alt;

    // If self is the identity, return rhs. Otherwise, return self. The other cases will be
    // predicated on neither self nor rhs being the identity.

    // let f1 = self.is_identity();
    const int f1 = fp_is_zero(&self->z);

    // let res = G1Projective::conditional_select(self, &G1Projective::from(rhs), f1);
    //g1j_cpy (res, self);
    g1j_cset(res, rhs, f1);

    // let f2 = rhs.is_identity();
    const int f2 = fp_is_zero(&rhs->z);

    // If neither is the identity but x1 = x2 and y1 != y2, then return the identity

    // let z = rhs.z.square();
    fp_mulred(z, &rhs->z, &rhs->z);

    // let u1 = self.x * z;
    fp_mulred(u1, &self->x, z);

    // let z = z * rhs.z;
    fp_mulred(z, z, &rhs->z);

    // let s1 = self.y * z;
    fp_mulred(s1, &self->y, z);

    // let z = self.z.square();
    fp_mulred(z, &self->z, &self->z);

    // let u2 = rhs.x * z;
    fp_mulred(u2, &rhs->x, z);

    // let z = z * self.z;
    fp_mulred(z, z, &self->z);

    // let s2 = rhs.y * z;
    fp_mulred(s2, &rhs->y, z);

    // let f3 = u1.ct_eq(&u2) & (!s1.ct_eq(&s2));
    const int f3 = fp_eq(u1, u2) & !fp_eq(s1, s2);

    // let res = G1Projective::conditional_select(&res, &G1Projective::identity(), (!f1) & (!f2) & f3);
    g1j_cset(res, &g1j_identity, (!f1) & (!f2) & f3);

    // let t = u1 + u2;
    fp_sum(t, u1, u2);

    // let m = s1 + s2;
    fp_sum(m, s1, s2);

    // let rr = t.square();
    fp_mulred(rr, t, t);

    // let m_alt = -u2;
    fp_neg(m_alt, u2);

    // let tt = u1 * m_alt;
    fp_mulred(tt, u1, m_alt);

    // let rr = rr + tt;
    fp_sum(rr, rr, tt);

    // Correct for x1 != x2 but y1 = -y2, which can occur because p - 1 is divisible by 3.
    // libsecp256k1 does this by substituting in an alternative (defined) expression for lambda.

    // let degenerate = m.is_zero() & rr.is_zero();
    const int degenerate = fp_is_zero(m) & fp_is_zero(rr);

    // let rr_alt = s1 + s1;
    fp_sum(rr_alt, s1, s1);

    // let m_alt = m_alt + u1;
    fp_sum(m_alt, m_alt, u1);

    // let rr_alt = Fp::conditional_select(&rr_alt, &rr, !degenerate);
    fp_cset(rr_alt, rr, !degenerate);

    // let m_alt = Fp::conditional_select(&m_alt, &m, !degenerate);
    fp_cset(m_alt, m, !degenerate);

    // let n = m_alt.square();
    fp_mulred(n, m_alt, m_alt);

    // let q = n * t;
    fp_mulred(q, n, t);

    // let n = n.square();
    fp_mulred(n, n, n);

    // let n = Fp::conditional_select(&n, &m, degenerate);
    fp_cset(n, m, degenerate);

    // let t = rr_alt.square();
    fp_mulred(T, rr_alt, rr_alt);

    // let z3 = m_alt * self.z * rhs.z; // We allow rhs.z != 1, so we must account for this.
    fp_mulred(&z3, m_alt, &self->z);
    fp_mulred(&z3, &z3, &rhs->z);

    // let z3 = z3 + z3;
    fp_sum(&z3, &z3, &z3);

    // let q = -q;
    fp_neg(q, q);

    // let t = t + q;
    // let x3 = t;
    // let t = t + t;
    fp_sum(&x3, T, q);
    fp_sum(T, &x3, &x3);

    // let t = t + q;
    fp_sum(T, T, q);

    // let t = t * rr_alt;
    fp_mulred(T, T, rr_alt);

    // let t = t + n;
    fp_sum(T, T, n);

    // let y3 = -t;
    fp_neg(&y3, T);

    // let x3 = x3 + x3;
    // let x3 = x3 + x3;
    fp_sum(&x3, &x3, &x3);
    fp_sum(&x3, &x3, &x3);

    // let y3 = y3 + y3;
    // let y3 = y3 + y3;
    fp_sum(&y3, &y3, &y3);
    fp_sum(&y3, &y3, &y3);

    // let tmp = G1Projective { x: x3, y: y3, z: z3, };

    // G1Projective::conditional_select(&res, &tmp, (!f1) & (!f2) & (!f3))
    {
	int f = (!f1) & (!f2) & (!f3);

	fp_cset(&res->x, &x3, f);
	fp_cset(&res->y, &y3, f);
	fp_cset(&res->z, &z3, f);
    }
}

extern void g1a_print(const char *, const G1Affine *);
extern void g1j_print(const char *, const G1Jacobian *);

#if 0
// point multiplication functions with 2-bit window

static G1Affine lut[4];
static G1Jacobian g1j_t0, g1j_t1;

void g1a_mul(G1Jacobian *j, const G1Affine *a, const Fq *by)
{
    {
	Fp r, t, u;
	G1Jacobian *x2 = &g1j_t0;
	G1Jacobian *x3 = &g1j_t1;

	g1a_dbl(x2, a);		// x2 = 2*a
	g1j_cpy(x3, x2);	// x3 = x2
	g1j_add_a(x3, a);	// x3 += a

	fp_mulred(&t, &x2->z, &x3->z);	// t = x2.z * x3.z
	fp_inv(&r, &t);			// r = 1/t

	fp_mulred(&t, &x2->z, &r);	// t = 1 / x3.z
	fp_mulred(&u, &t, &t);		// u = 1 / x3.z^2
	fp_mulred(&t, &u, &t);		// t = 1 / x3.z^3

	fp_mulred(&lut[3].x, &x3->x, &u);
	fp_mulred(&lut[3].y, &x3->y, &t);
	lut[3].inf = fp_is_zero(&r);

	fp_mulred(&t, &x3->z, &r);	// t = 1 / x2.z
	fp_mulred(&u, &t, &t);		// u = 1 / x2.z^2
	fp_mulred(&t, &u, &t);		// t = 1 / x2.z^3

	fp_mulred(&lut[2].x, &x2->x, &u);
	fp_mulred(&lut[2].y, &x2->y, &t);
	lut[2].inf = fp_is_zero(&r);
    }

    g1a_cpy(&lut[1], a);
    g1a_cpy(&lut[0], &g1a_identity);

    uint32_t i = by->x[7];

    g1j_from_a(j, &lut[(i >> 30) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 28) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 26) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 24) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 22) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 20) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 18) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 16) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 14) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 12) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 10) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  8) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  6) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  4) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  2) & 3]); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  0) & 3]);

    for (int w=6; w>=0; w--)
    {
	i = by->x[w];
						g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 30) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 28) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 26) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 24) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 22) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 20) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 18) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 16) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 14) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 12) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 10) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  8) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  6) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  4) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  2) & 3]);	g1j_dbl(j, j);	g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  0) & 3]);
    }
}

void g1j_mul(G1Jacobian *k, const G1Jacobian *j, const Fq *by)
{
    const G1Jacobian *x1 = j;

    {
	Fp r, t, u, v;
	G1Jacobian *x2 = &g1j_t0;
	G1Jacobian *x3 = &g1j_t1;

	g1j_dbl(x2, x1);	// x2 = 2*x1
	g1j_cpy(x3, x2);	// x3 = x2
	g1j_add_j(x3, x1);	// x3 += x1

	fp_mulred(&t, &x1->z, &x2->z);	// t = x1.z * x2.z
	fp_mulred(&t, &x3->z, &t);	// t = x3.z * t
	fp_inv(&r, &t);			// r = 1/t = 1/(x1.z * x2.z * x3.z)

	fp_mulred(&t, &x1->z, &r);	// t = 1 / (x2.z * x3.z)
	fp_mulred(&u, &x2->z, &t);	// u = 1 / x3.z
	fp_mulred(&v, &u, &u);		// v = 1 / x3.z^2
	fp_mulred(&u, &v, &u);		// u = 1 / x3.z^3

	fp_mulred(&lut[3].x, &x3->x, &v);
	fp_mulred(&lut[3].y, &x3->y, &u);
	lut[3].inf = fp_is_zero(&r);

	fp_mulred(&u, &x3->z, &t);	// u = 1 / x2.z
	fp_mulred(&v, &u, &u);		// v = 1 / x2.z^2
	fp_mulred(&u, &v, &u);		// u = 1 / x2.z^3

	fp_mulred(&lut[2].x, &x2->x, &v);
	fp_mulred(&lut[2].y, &x2->y, &u);
	lut[2].inf = fp_is_zero(&r);

	fp_mulred(&t, &x2->z, &r);	// t = 1 / (x1.z * x3.z)
	fp_mulred(&u, &x3->z, &t);	// u = 1 / x1.z
	fp_mulred(&v, &u, &u);		// v = 1 / x1.z^2
	fp_mulred(&u, &v, &u);		// u = 1 / x1.z^3

	fp_mulred(&lut[1].x, &x1->x, &v);
	fp_mulred(&lut[1].y, &x1->y, &u);
	lut[1].inf = fp_is_zero(&r);
    }

    g1a_cpy(&lut[0], &g1a_identity);

    uint32_t i = by->x[7];

    g1j_from_a(k, &lut[(i >> 30) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 28) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 26) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 24) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 22) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 20) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 18) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 16) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 14) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 12) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 10) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  8) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  6) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  4) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  2) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  0) & 3]);

    for (int w=6; w>=0; w--)
    {
	i = by->x[w];
						g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 30) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 28) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 26) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 24) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 22) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 20) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 18) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 16) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 14) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 12) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 10) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  8) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  6) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  4) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  2) & 3]);	g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  0) & 3]);
    }
}
#else
// point multiplication functions with 4-bit window

static G1Affine lut[16];
static G1Jacobian g1j_t0, g1j_t1;

void g1a_mul(G1Jacobian *j, const G1Affine *a, const Fq *by)
{
    g1a_cpy(&lut[1], a);

    // Compute the points in the lookup table twice to remove the need to store
    // them all the same time in both affine and jacobian coordinates. Instead,
    // store all of their z coordinates. Each point must be computed in the
    // exact same way in both passes so that its coordinates stay the same.
    //
    // Calculation sequence used for base point multiples:
    //
    //  2=2*1
    //  4=2*2	store
    //  8=2*4
    //  9=8+1
    //
    //  5=4+1
    // 10=2*5
    // 11=10+1
    //
    //  3=4-1
    //  6=2*3	store (replace 4)
    // 12=2*6
    // 13=12+1
    //
    //  7=6+1
    // 14=2*7
    // 15=14+1

    {
	G1Jacobian  *t = &g1j_t0;
	G1Jacobian  *u = &g1j_t1;
	G1Affine    *b = &lut[0];

	g1a_neg(b, a);

	//////////////////////////////////////////////////
	// First pass - compute and store z coordinates
	//////////////////////////////////////////////////

	g1a_dbl  (t, a);	fp_cpy(&lut[ 2].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 4].x, &t->z);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 8].x, &t->z);
	g1j_add_a(t, a);	fp_cpy(&lut[ 9].x, &t->z);

	g1j_cpy  (t, u);

	g1j_add_a(t, a);	fp_cpy(&lut[ 5].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[10].x, &t->z);
	g1j_add_a(t, a);	fp_cpy(&lut[11].x, &t->z);

	g1j_cpy  (t, u);

	g1j_add_a(t, b);	fp_cpy(&lut[ 3].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 6].x, &t->z);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_cpy(&lut[12].x, &t->z);
	g1j_add_a(t, a);	fp_cpy(&lut[13].x, &t->z);

	g1j_cpy  (t, u);

	g1j_add_a(t, a);	fp_cpy(&lut[ 7].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[14].x, &t->z);
	g1j_add_a(t, a);	fp_cpy(&lut[15].x, &t->z);

	// Compute products of z coordinates of lower point multiples

	fp_cpy   (&lut[ 3].y, &lut[ 2].x);		// z2
	fp_mulred(&lut[ 4].y, &lut[ 3].x, &lut[ 3].y);	// z2 * z3
	fp_mulred(&lut[ 5].y, &lut[ 4].x, &lut[ 4].y);	// z2 * z3 * z4
	fp_mulred(&lut[ 6].y, &lut[ 5].x, &lut[ 5].y);	// z2 * ... * z5
	fp_mulred(&lut[ 7].y, &lut[ 6].x, &lut[ 6].y);	// z2 * ... * z6
	fp_mulred(&lut[ 8].y, &lut[ 7].x, &lut[ 7].y);	// z2 * ... * z7
	fp_mulred(&lut[ 9].y, &lut[ 8].x, &lut[ 8].y);	// z2 * ... * z8
	fp_mulred(&lut[10].y, &lut[ 9].x, &lut[ 9].y);	// z2 * ... * z9
	fp_mulred(&lut[11].y, &lut[10].x, &lut[10].y);	// z2 * ... * z10
	fp_mulred(&lut[12].y, &lut[11].x, &lut[11].y);	// z2 * ... * z11
	fp_mulred(&lut[13].y, &lut[12].x, &lut[12].y);	// z2 * ... * z12
	fp_mulred(&lut[14].y, &lut[13].x, &lut[13].y);	// z2 * ... * z13
	fp_mulred(&lut[15].y, &lut[14].x, &lut[14].y);	// z2 * ... * z14

	// Compute products of z coordinates of same or higher point multiples

	fp_mulred(&lut[14].x, &lut[14].x, &lut[15].x);	// z15 * z14
	fp_mulred(&lut[13].x, &lut[13].x, &lut[14].x);	// z15 * z14 * z13
	fp_mulred(&lut[12].x, &lut[12].x, &lut[13].x);	// z15 * ... * z12
	fp_mulred(&lut[11].x, &lut[11].x, &lut[12].x);	// z15 * ... * z11
	fp_mulred(&lut[10].x, &lut[10].x, &lut[11].x);	// z15 * ... * z10
	fp_mulred(&lut[ 9].x, &lut[ 9].x, &lut[10].x);	// z15 * ... * z9
	fp_mulred(&lut[ 8].x, &lut[ 8].x, &lut[ 9].x);	// z15 * ... * z8
	fp_mulred(&lut[ 7].x, &lut[ 7].x, &lut[ 8].x);	// z15 * ... * z7
	fp_mulred(&lut[ 6].x, &lut[ 6].x, &lut[ 7].x);	// z15 * ... * z6
	fp_mulred(&lut[ 5].x, &lut[ 5].x, &lut[ 6].x);	// z15 * ... * z5
	fp_mulred(&lut[ 4].x, &lut[ 4].x, &lut[ 5].x);	// z15 * ... * z4
	fp_mulred(&lut[ 3].x, &lut[ 3].x, &lut[ 4].x);	// z15 * ... * z3
	fp_mulred(&lut[ 2].x, &lut[ 2].x, &lut[ 3].x);	// z15 * ... * z2

	// Compute the inverse of (z2 * ... * z15)

	fp_inv(&lut[2].y, &lut[2].x);

	// Compute all inverses and store them in the y coordinates

	fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[ 4].x);	// (z2            )   ( z4 * ... * z15)
	fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[ 5].x);	// (z2 * z3       ) * ( z5 * ... * z15)
	fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[ 6].x);	// (z2 * ... *  z4) * ( z6 * ... * z15)
	fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[ 7].x);	// (z2 * ... *  z5) * ( z7 * ... * z15)
	fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[ 8].x);	// (z2 * ... *  z6) * ( z8 * ... * z15)
	fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[ 9].x);	// (z2 * ... *  z7) * ( z9 * ... * z15)
	fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[10].x);	// (z2 * ... *  z8) * (z10 * ... * z15)
	fp_mulred(&lut[10].y, &lut[10].y, &lut[11].x);	// (z2 * ... *  z9) * (z11 * ... * z15)
	fp_mulred(&lut[11].y, &lut[11].y, &lut[12].x);	// (z2 * ... * z10) * (z12 * ... * z15)
	fp_mulred(&lut[12].y, &lut[12].y, &lut[13].x);	// (z2 * ... * z11) * (z13 * ... * z15)
	fp_mulred(&lut[13].y, &lut[13].y, &lut[14].x);	// (z2 * ... * z12) * (      z14 * z15)
	fp_mulred(&lut[14].y, &lut[14].y, &lut[15].x);	// (z2 * ... * z13) * (            z15)

	fp_mulred(&lut[15].y, &lut[15].y, &lut[2].y);	// 1 / z15
	fp_mulred(&lut[14].y, &lut[14].y, &lut[2].y);	// 1 / z14
	fp_mulred(&lut[13].y, &lut[13].y, &lut[2].y);	// 1 / z13
	fp_mulred(&lut[12].y, &lut[12].y, &lut[2].y);	// 1 / z12
	fp_mulred(&lut[11].y, &lut[11].y, &lut[2].y);	// 1 / z11
	fp_mulred(&lut[10].y, &lut[10].y, &lut[2].y);	// 1 / z10
	fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[2].y);	// 1 /  z9
	fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[2].y);	// 1 /  z8
	fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[2].y);	// 1 /  z7
	fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[2].y);	// 1 /  z6
	fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[2].y);	// 1 /  z5
	fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[2].y);	// 1 /  z4
	fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[2].y);	// 1 /  z3
	fp_mulred(&lut[ 2].y, &lut[ 2].y, &lut[3].x);	// 1 /  z2

	// Set infinity flags

	lut[ 1].inf =
	lut[ 2].inf =
	lut[ 3].inf =
	lut[ 4].inf =
	lut[ 5].inf =
	lut[ 6].inf =
	lut[ 7].inf =
	lut[ 8].inf =
	lut[ 9].inf =
	lut[10].inf =
	lut[11].inf =
	lut[12].inf =
	lut[13].inf =
	lut[14].inf =
	lut[15].inf = fp_is_zero(&lut[2].x);

	// Prepare scaling factors (square and cube of each inverse)

	fp_mulred(&lut[ 2].x, &lut[ 2].y, &lut[ 2].y); fp_mulred(&lut[ 2].y, &lut[ 2].y, &lut[ 2].x);
	fp_mulred(&lut[ 3].x, &lut[ 3].y, &lut[ 3].y); fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[ 3].x);
	fp_mulred(&lut[ 4].x, &lut[ 4].y, &lut[ 4].y); fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[ 4].x);
	fp_mulred(&lut[ 5].x, &lut[ 5].y, &lut[ 5].y); fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[ 5].x);
	fp_mulred(&lut[ 6].x, &lut[ 6].y, &lut[ 6].y); fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[ 6].x);
	fp_mulred(&lut[ 7].x, &lut[ 7].y, &lut[ 7].y); fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[ 7].x);
	fp_mulred(&lut[ 8].x, &lut[ 8].y, &lut[ 8].y); fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[ 8].x);
	fp_mulred(&lut[ 9].x, &lut[ 9].y, &lut[ 9].y); fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[ 9].x);
	fp_mulred(&lut[10].x, &lut[10].y, &lut[10].y); fp_mulred(&lut[10].y, &lut[10].y, &lut[10].x);
	fp_mulred(&lut[11].x, &lut[11].y, &lut[11].y); fp_mulred(&lut[11].y, &lut[11].y, &lut[11].x);
	fp_mulred(&lut[12].x, &lut[12].y, &lut[12].y); fp_mulred(&lut[12].y, &lut[12].y, &lut[12].x);
	fp_mulred(&lut[13].x, &lut[13].y, &lut[13].y); fp_mulred(&lut[13].y, &lut[13].y, &lut[13].x);
	fp_mulred(&lut[14].x, &lut[14].y, &lut[14].y); fp_mulred(&lut[14].y, &lut[14].y, &lut[14].x);
	fp_mulred(&lut[15].x, &lut[15].y, &lut[15].y); fp_mulred(&lut[15].y, &lut[15].y, &lut[15].x);

	//////////////////////////////////////////////////
	// Second pass - Affinify and store points
	//////////////////////////////////////////////////

	g1a_dbl  (t, a);	fp_mulred(&lut[ 2].x, &lut[ 2].x, &t->x);	fp_mulred(&lut[ 2].y, &lut[ 2].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 4].x, &lut[ 4].x, &t->x);	fp_mulred(&lut[ 4].y, &lut[ 4].y, &t->y);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 8].x, &lut[ 8].x, &t->x);	fp_mulred(&lut[ 8].y, &lut[ 8].y, &t->y);
	g1j_add_a(t, a);	fp_mulred(&lut[ 9].x, &lut[ 9].x, &t->x);	fp_mulred(&lut[ 9].y, &lut[ 9].y, &t->y);

	g1j_cpy  (t, u);

	g1j_add_a(t, a);	fp_mulred(&lut[ 5].x, &lut[ 5].x, &t->x);	fp_mulred(&lut[ 5].y, &lut[ 5].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[10].x, &lut[10].x, &t->x);	fp_mulred(&lut[10].y, &lut[10].y, &t->y);
	g1j_add_a(t, a);	fp_mulred(&lut[11].x, &lut[11].x, &t->x);	fp_mulred(&lut[11].y, &lut[11].y, &t->y);

	g1j_cpy  (t, u);

	g1j_add_a(t, b);	fp_mulred(&lut[ 3].x, &lut[ 3].x, &t->x);	fp_mulred(&lut[ 3].y, &lut[ 3].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 6].x, &lut[ 6].x, &t->x);	fp_mulred(&lut[ 6].y, &lut[ 6].y, &t->y);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_mulred(&lut[12].x, &lut[12].x, &t->x);	fp_mulred(&lut[12].y, &lut[12].y, &t->y);
	g1j_add_a(t, a);	fp_mulred(&lut[13].x, &lut[13].x, &t->x);	fp_mulred(&lut[13].y, &lut[13].y, &t->y);

	g1j_cpy  (t, u);

	g1j_add_a(t, a);	fp_mulred(&lut[ 7].x, &lut[ 7].x, &t->x);	fp_mulred(&lut[ 7].y, &lut[ 7].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[14].x, &lut[14].x, &t->x);	fp_mulred(&lut[14].y, &lut[14].y, &t->y);
	g1j_add_a(t, a);	fp_mulred(&lut[15].x, &lut[15].x, &t->x);	fp_mulred(&lut[15].y, &lut[15].y, &t->y);
    }

    g1a_cpy(&lut[0], &g1a_identity);

    uint32_t i = by->x[7];

    g1j_from_a(j, &lut[(i >> 28) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 24) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 20) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 16) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >> 12) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  8) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  4) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
    g1j_add_a (j, &lut[(i >>  0) & 15]);

    for (int w=6; w>=0; w--)
    {
	i = by->x[w];
						g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 28) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 24) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 20) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 16) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >> 12) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  8) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  4) & 15]);	g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j); g1j_dbl(j, j);
	g1j_add_a(j, &lut[(i >>  0) & 15]);
    }
}

void g1j_mul(G1Jacobian *k, const G1Jacobian *j, const Fq *by)
{
    // Compute the points in the lookup table twice to remove the need to store
    // them all the same time in both affine and jacobian coordinates. Instead,
    // store all of their z coordinates. Each point must be computed in the
    // exact same way in both passes so that its coordinates stay the same.
    //
    // Calculation sequence used for base point multiples:
    //
    //  2=2*1
    //  4=2*2	store
    //  8=2*4
    //  9=8+1
    //
    //  5=4+1
    // 10=2*5
    // 11=10+1
    //
    //  3=4-1
    //  6=2*3	store (replace 4)
    // 12=2*6
    // 13=12+1
    //
    //  7=6+1
    // 14=2*7
    // 15=14+1

    {
        G1Jacobian  *t = &g1j_t0;
        G1Jacobian  *u = &g1j_t1;

	//////////////////////////////////////////////////
	// First pass - compute and store z coordinates
	//////////////////////////////////////////////////

				fp_cpy(&lut[ 1].x, &j->z);
	g1j_dbl  (t, j);	fp_cpy(&lut[ 2].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 4].x, &t->z);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 8].x, &t->z);
	g1j_add_j(t, j);	fp_cpy(&lut[ 9].x, &t->z);

	g1j_cpy  (t, u);

	g1j_add_j(t, j);	fp_cpy(&lut[ 5].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[10].x, &t->z);
	g1j_add_j(t, j);	fp_cpy(&lut[11].x, &t->z);

	g1j_cpy  (t, u);

	g1j_neg  (u, j);

	g1j_add_j(t, u);	fp_cpy(&lut[ 3].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[ 6].x, &t->z);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_cpy(&lut[12].x, &t->z);
	g1j_add_j(t, j);	fp_cpy(&lut[13].x, &t->z);

	g1j_cpy  (t, u);

	g1j_add_j(t, j);	fp_cpy(&lut[ 7].x, &t->z);
	g1j_dbl  (t, t);	fp_cpy(&lut[14].x, &t->z);
	g1j_add_j(t, j);	fp_cpy(&lut[15].x, &t->z);

	// Compute products of z coordinates of lower point multiples

	fp_cpy   (&lut[ 2].y, &lut[ 1].x);		// z1
	fp_mulred(&lut[ 3].y, &lut[ 2].x, &lut[ 2].y);	// z1 * z2
	fp_mulred(&lut[ 4].y, &lut[ 3].x, &lut[ 3].y);	// z1 * z2 * z3
	fp_mulred(&lut[ 5].y, &lut[ 4].x, &lut[ 4].y);	// z1 * ... * z4
	fp_mulred(&lut[ 6].y, &lut[ 5].x, &lut[ 5].y);	// z1 * ... * z5
	fp_mulred(&lut[ 7].y, &lut[ 6].x, &lut[ 6].y);	// z1 * ... * z6
	fp_mulred(&lut[ 8].y, &lut[ 7].x, &lut[ 7].y);	// z1 * ... * z7
	fp_mulred(&lut[ 9].y, &lut[ 8].x, &lut[ 8].y);	// z1 * ... * z8
	fp_mulred(&lut[10].y, &lut[ 9].x, &lut[ 9].y);	// z1 * ... * z9
	fp_mulred(&lut[11].y, &lut[10].x, &lut[10].y);	// z1 * ... * z10
	fp_mulred(&lut[12].y, &lut[11].x, &lut[11].y);	// z1 * ... * z11
	fp_mulred(&lut[13].y, &lut[12].x, &lut[12].y);	// z1 * ... * z12
	fp_mulred(&lut[14].y, &lut[13].x, &lut[13].y);	// z1 * ... * z13
	fp_mulred(&lut[15].y, &lut[14].x, &lut[14].y);	// z1 * ... * z14

	// Compute products of z coordinates of same or higher point multiples

	fp_mulred(&lut[14].x, &lut[14].x, &lut[15].x);	// z15 * z14
	fp_mulred(&lut[13].x, &lut[13].x, &lut[14].x);	// z15 * z14 * z13
	fp_mulred(&lut[12].x, &lut[12].x, &lut[13].x);	// z15 * ... * z12
	fp_mulred(&lut[11].x, &lut[11].x, &lut[12].x);	// z15 * ... * z11
	fp_mulred(&lut[10].x, &lut[10].x, &lut[11].x);	// z15 * ... * z10
	fp_mulred(&lut[ 9].x, &lut[ 9].x, &lut[10].x);	// z15 * ... * z9
	fp_mulred(&lut[ 8].x, &lut[ 8].x, &lut[ 9].x);	// z15 * ... * z8
	fp_mulred(&lut[ 7].x, &lut[ 7].x, &lut[ 8].x);	// z15 * ... * z7
	fp_mulred(&lut[ 6].x, &lut[ 6].x, &lut[ 7].x);	// z15 * ... * z6
	fp_mulred(&lut[ 5].x, &lut[ 5].x, &lut[ 6].x);	// z15 * ... * z5
	fp_mulred(&lut[ 4].x, &lut[ 4].x, &lut[ 5].x);	// z15 * ... * z4
	fp_mulred(&lut[ 3].x, &lut[ 3].x, &lut[ 4].x);	// z15 * ... * z3
	fp_mulred(&lut[ 2].x, &lut[ 2].x, &lut[ 3].x);	// z15 * ... * z2
	fp_mulred(&lut[ 1].x, &lut[ 1].x, &lut[ 2].x);	// z15 * ... * z1

	// Compute the inverse of (z1 * ... * z15)

	fp_inv(&lut[1].y, &lut[1].x);

	// Compute all inverses and store them in the y coordinates

	fp_mulred(&lut[ 2].y, &lut[ 2].y, &lut[ 3].x);	// (z1            )   ( z3 * ... * z15)
	fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[ 4].x);	// (z1 *  z2      )   ( z4 * ... * z15)
	fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[ 5].x);	// (z1 *  z2 *  z3) * ( z5 * ... * z15)
	fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[ 6].x);	// (z1 * ... *  z4) * ( z6 * ... * z15)
	fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[ 7].x);	// (z1 * ... *  z5) * ( z7 * ... * z15)
	fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[ 8].x);	// (z1 * ... *  z6) * ( z8 * ... * z15)
	fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[ 9].x);	// (z1 * ... *  z7) * ( z9 * ... * z15)
	fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[10].x);	// (z1 * ... *  z8) * (z10 * ... * z15)
	fp_mulred(&lut[10].y, &lut[10].y, &lut[11].x);	// (z1 * ... *  z9) * (z11 * ... * z15)
	fp_mulred(&lut[11].y, &lut[11].y, &lut[12].x);	// (z1 * ... * z10) * (z12 * ... * z15)
	fp_mulred(&lut[12].y, &lut[12].y, &lut[13].x);	// (z1 * ... * z11) * (z13 * ... * z15)
	fp_mulred(&lut[13].y, &lut[13].y, &lut[14].x);	// (z1 * ... * z12) * (      z14 * z15)
	fp_mulred(&lut[14].y, &lut[14].y, &lut[15].x);	// (z1 * ... * z13) * (            z15)

	fp_mulred(&lut[15].y, &lut[15].y, &lut[1].y);	// 1 / z15
	fp_mulred(&lut[14].y, &lut[14].y, &lut[1].y);	// 1 / z14
	fp_mulred(&lut[13].y, &lut[13].y, &lut[1].y);	// 1 / z13
	fp_mulred(&lut[12].y, &lut[12].y, &lut[1].y);	// 1 / z12
	fp_mulred(&lut[11].y, &lut[11].y, &lut[1].y);	// 1 / z11
	fp_mulred(&lut[10].y, &lut[10].y, &lut[1].y);	// 1 / z10
	fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[1].y);	// 1 /  z9
	fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[1].y);	// 1 /  z8
	fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[1].y);	// 1 /  z7
	fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[1].y);	// 1 /  z6
	fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[1].y);	// 1 /  z5
	fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[1].y);	// 1 /  z4
	fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[1].y);	// 1 /  z3
	fp_mulred(&lut[ 2].y, &lut[ 2].y, &lut[1].y);	// 1 /  z2
	fp_mulred(&lut[ 1].y, &lut[ 1].y, &lut[2].x);	// 1 /  z1

	// Set infinity flags

	lut[ 1].inf =
	lut[ 2].inf =
	lut[ 3].inf =
	lut[ 4].inf =
	lut[ 5].inf =
	lut[ 6].inf =
	lut[ 7].inf =
	lut[ 8].inf =
	lut[ 9].inf =
	lut[10].inf =
	lut[11].inf =
	lut[12].inf =
	lut[13].inf =
	lut[14].inf =
	lut[15].inf = fp_is_zero(&lut[1].x);

	// Prepare scaling factors (square and cube of each inverse)

	fp_mulred(&lut[ 1].x, &lut[ 1].y, &lut[ 1].y); fp_mulred(&lut[ 1].y, &lut[ 1].y, &lut[ 1].x);
	fp_mulred(&lut[ 2].x, &lut[ 2].y, &lut[ 2].y); fp_mulred(&lut[ 2].y, &lut[ 2].y, &lut[ 2].x);
	fp_mulred(&lut[ 3].x, &lut[ 3].y, &lut[ 3].y); fp_mulred(&lut[ 3].y, &lut[ 3].y, &lut[ 3].x);
	fp_mulred(&lut[ 4].x, &lut[ 4].y, &lut[ 4].y); fp_mulred(&lut[ 4].y, &lut[ 4].y, &lut[ 4].x);
	fp_mulred(&lut[ 5].x, &lut[ 5].y, &lut[ 5].y); fp_mulred(&lut[ 5].y, &lut[ 5].y, &lut[ 5].x);
	fp_mulred(&lut[ 6].x, &lut[ 6].y, &lut[ 6].y); fp_mulred(&lut[ 6].y, &lut[ 6].y, &lut[ 6].x);
	fp_mulred(&lut[ 7].x, &lut[ 7].y, &lut[ 7].y); fp_mulred(&lut[ 7].y, &lut[ 7].y, &lut[ 7].x);
	fp_mulred(&lut[ 8].x, &lut[ 8].y, &lut[ 8].y); fp_mulred(&lut[ 8].y, &lut[ 8].y, &lut[ 8].x);
	fp_mulred(&lut[ 9].x, &lut[ 9].y, &lut[ 9].y); fp_mulred(&lut[ 9].y, &lut[ 9].y, &lut[ 9].x);
	fp_mulred(&lut[10].x, &lut[10].y, &lut[10].y); fp_mulred(&lut[10].y, &lut[10].y, &lut[10].x);
	fp_mulred(&lut[11].x, &lut[11].y, &lut[11].y); fp_mulred(&lut[11].y, &lut[11].y, &lut[11].x);
	fp_mulred(&lut[12].x, &lut[12].y, &lut[12].y); fp_mulred(&lut[12].y, &lut[12].y, &lut[12].x);
	fp_mulred(&lut[13].x, &lut[13].y, &lut[13].y); fp_mulred(&lut[13].y, &lut[13].y, &lut[13].x);
	fp_mulred(&lut[14].x, &lut[14].y, &lut[14].y); fp_mulred(&lut[14].y, &lut[14].y, &lut[14].x);
	fp_mulred(&lut[15].x, &lut[15].y, &lut[15].y); fp_mulred(&lut[15].y, &lut[15].y, &lut[15].x);

	//////////////////////////////////////////////////
	// Second pass - Affinify and store points
	//////////////////////////////////////////////////

				fp_mulred(&lut[ 1].x, &lut[ 1].x, &j->x);	fp_mulred(&lut[ 1].y, &lut[ 1].y, &j->y);
	g1j_dbl  (t, j);	fp_mulred(&lut[ 2].x, &lut[ 2].x, &t->x);	fp_mulred(&lut[ 2].y, &lut[ 2].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 4].x, &lut[ 4].x, &t->x);	fp_mulred(&lut[ 4].y, &lut[ 4].y, &t->y);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 8].x, &lut[ 8].x, &t->x);	fp_mulred(&lut[ 8].y, &lut[ 8].y, &t->y);
	g1j_add_j(t, j);	fp_mulred(&lut[ 9].x, &lut[ 9].x, &t->x);	fp_mulred(&lut[ 9].y, &lut[ 9].y, &t->y);

	g1j_cpy  (t, u);

	g1j_add_j(t, j);	fp_mulred(&lut[ 5].x, &lut[ 5].x, &t->x);	fp_mulred(&lut[ 5].y, &lut[ 5].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[10].x, &lut[10].x, &t->x);	fp_mulred(&lut[10].y, &lut[10].y, &t->y);
	g1j_add_j(t, j);	fp_mulred(&lut[11].x, &lut[11].x, &t->x);	fp_mulred(&lut[11].y, &lut[11].y, &t->y);

	g1j_cpy  (t, u);

	g1j_neg  (u, j);

	g1j_add_j(t, u);	fp_mulred(&lut[ 3].x, &lut[ 3].x, &t->x);	fp_mulred(&lut[ 3].y, &lut[ 3].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[ 6].x, &lut[ 6].x, &t->x);	fp_mulred(&lut[ 6].y, &lut[ 6].y, &t->y);	g1j_cpy(u, t);
	g1j_dbl  (t, t);	fp_mulred(&lut[12].x, &lut[12].x, &t->x);	fp_mulred(&lut[12].y, &lut[12].y, &t->y);
	g1j_add_j(t, j);	fp_mulred(&lut[13].x, &lut[13].x, &t->x);	fp_mulred(&lut[13].y, &lut[13].y, &t->y);

	g1j_cpy  (t, u);

	g1j_add_j(t, j);	fp_mulred(&lut[ 7].x, &lut[ 7].x, &t->x);	fp_mulred(&lut[ 7].y, &lut[ 7].y, &t->y);
	g1j_dbl  (t, t);	fp_mulred(&lut[14].x, &lut[14].x, &t->x);	fp_mulred(&lut[14].y, &lut[14].y, &t->y);
	g1j_add_j(t, j);	fp_mulred(&lut[15].x, &lut[15].x, &t->x);	fp_mulred(&lut[15].y, &lut[15].y, &t->y);
    }

    g1a_cpy(&lut[0], &g1a_identity);

    uint32_t i = by->x[7];

    g1j_from_a(k, &lut[(i >> 28) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 24) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 20) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 16) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >> 12) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  8) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  4) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
    g1j_add_a (k, &lut[(i >>  0) & 15]);

    for (int w=6; w>=0; w--)
    {
	i = by->x[w];
						g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 28) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 24) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 20) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 16) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >> 12) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  8) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  4) & 15]);	g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k); g1j_dbl(k, k);
	g1j_add_a(k, &lut[(i >>  0) & 15]);
    }
}
#endif
