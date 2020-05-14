#include <stdio.h>
#include "fp.h"
#include "fp_mont.h"

//////////////////////////////
// Constants
//////////////////////////////

const Fp fp_zero	= FP_ZERO;
const Fp fp_one		= FP_ONE;	// R = 2^384 mod p = Montgomery representation of 1
const Fp fp_r		= FP_R;		// R^2 = 2^768 mod p = Montgomery representation of R
const Fp fp_p		= FP_P;		// p (the modulus)
const Fp fp_pm2		= FP_PM2;	// pm2 = p-2

//////////////////////////////
// I/O with conversion
//////////////////////////////

void fp_conv_to_bytes	(uint8_t *p, const Fp *x)
{
    uint32_t tmp[24];
    Fp t;

    for (int i=0; i<12; ++i)
    {
	tmp[i] = x->x[i];
	tmp[i+12] = 0;
    }

    fp_redc(t.x, tmp);
    fp_to_bytes(p, &t);
}

void fp_conv_from_bytes	(Fp *x, const uint8_t *p)
{
    fp_from_bytes(x, p);
    fp_mulred(x, x, &fp_r);
}

//////////////////////////////
// Inversion (x -> x^(p-2))
//////////////////////////////

void fp_inv(Fp *r, const Fp *x)
{
    Fp tmp[4];

    const Fp *x1 = x;
    Fp *x2 = &tmp[2];
    Fp *x3 = &tmp[0];
    Fp *x5 = &tmp[1];
    Fp *x7 = &tmp[2];
    Fp *t  = &tmp[3];	// could use t = r if r != x

    fp_mulred(x2, x1, x1);	// 10
    fp_mulred(x3, x2, x1);	// 11
    fp_mulred(x5, x3, x2);	// 101
    fp_mulred(x7, x5, x2);	// 111
    fp_cpy(t, x7);

// Squaring counts

#define  S1 fp_mulred(t, t, t)
#define  S2 S1;S1
#define  S3 S2;S1
#define  S4 S3;S1
#define  S5 S4;S1
#define  S6 S5;S1
#define  S7 S6;S1
#define  S8 S7;S1
#define  S9 S8;S1
#define S10 S9;S1

// Powers of x to multiply with

#define X1 fp_mulred(t, t, x1)
#define X3 fp_mulred(t, t, x3)
#define X5 fp_mulred(t, t, x5)
#define X7 fp_mulred(t, t, x7)

        X5;		// 1100
    S2; X5;		// 110101
    S2; X3;		// 11010111
    S5; X3;		// 1101011100011
    S3; X5;		// 1101011100011101
    S3; X1;		// 1101011100011101001
    S5; X3;		// ...00011
    S7; X5;		// ...0000101
    S2; X3;		// ...11
    S2; X3;		// ...11
    S6; X5;		// ...000101
    S4; X1;		// ...0001
    S7; X7;		// ...0000111
    S4; X5;		// ...0101
    S3; X3;		// ...011
    S5; X3;		// ...00011
    S6; X7;		// ...000111
    S3; X3;		// ...011
    S8; X5;		// ...00000101
    S2; X3;		// ...11
    S9; X3;		// ...000000011
    S3; X3;		// ...011
    S5; X5;		// ...00101
    S7; X5;		// ...0000101
    S3; X1;		// ...001
    S3; X1;		// ...001
    S5; X7;		// ...00111
    S3; X3;		// ...011
    S5; X3;		// ...00011
    S2; X1;		// ...01
    S4; X1;		// ...0001
    S4; X1;		// ...0001
    S3; X3;		// ...011
    S3; X3;		// ...011
    S4; X3;		// ...0011
    S3; X7;		// ...111
    S4; X3;		// ...0011
    S10; X3;		// ...0000000011
    S2; X3;		// ...11
    S4; X5;		// ...0101
    S4; X1;		// ...0001
    S5; X7;		// ...00111
    S5; X3;		// ...00011
    S2; X3;		// ...11
    S6; X7;		// ...000111
    S3; X5;		// ...101
    S3; X7;		// ...111
    S4; X3;		// ...0011
    S3; X3;		// ...011
    S4; X1;		// ...0001
    S6; X5;		// ...000101
    S2; X3;		// ...11
    S2; X3;		// ...11
    S3; X3;		// ...011
    S3; X5;		// ...101
    S6; X1;		// ...000001
    S5; X5;		// ...00101
    S3; X1;		// ...001
    S10; S7; X5;	// ...00000000000000101
    S2; X3;		// ...11
    S7; X5;		// ...0000101
    S3; X5;		// ...101
    S4; X7;		// ...0111
    S4; X5;		// ...0101
    S4; X1;		// ...0001
    S6; X3;		// ...000011
    S10; S10; S9; X1;	// ...00000000000000000000000000001
    S7; X5;		// ...0000101
    S5; X1;		// ...00001
    S6; X5;		// ...000101
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    S3; X7;		// ...111
    fp_cpy(r, t);
}

#undef S10
#undef S9
#undef S8
#undef S7
#undef S6
#undef S5
#undef S4
#undef S3
#undef S2
#undef S1

#undef X7
#undef X5
#undef X3
#undef X1

void fp_mulred(Fp *output, const Fp *left, const Fp *right)
{
    static uint32_t prod[24];

    fp_prod(prod, left, right);
    fp_redc(output->x, prod);
    fp_diff(output, output, &fp_p);
}
