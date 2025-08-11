/*
	This code is from the ax computer algebra system and is released into the public domain by Roman Pearce, June 2025.
	
	This software and documentation is provided "as is", without warranty of any kind, express or implied,
	including but not limited to the warranties of merchantability, fitness for a particular purpose, and 
	noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, 
	or other liability, whether in an action of contract, tort, or otherwise, arising from, out of or in 
	connection with the software or the use or other dealings in the software.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>


/*	Machine integer routines

	We assume the following non-portable things:
	- two's complement integer arithmetic
	- a signed right shift duplicates the sign bit
	- a flat memory model that is byte addressable
	- malloc aligns memory to the word size
*/

#define INT32 int
#define INT64 long long int
#define UINT32 unsigned int
#define UINT64 unsigned long long int
#define CHAR   unsigned char

/* determine word size */
#if UINTPTR_MAX==0xFFFFFFFF
	#define WORDSIZE 32
	typedef INT32	INT;
	typedef UINT32	UINT;
	typedef UINT32*	DAG;
#elif UINTPTR_MAX==0xFFFFFFFFFFFFFFFF
	#define WORDSIZE 64
	typedef INT64	INT;
	typedef UINT64	UINT;
	typedef UINT64*	DAG;
#else
	#error port WORDSIZE
#endif
#define I(x) ((INT)(x))
#define U(x) ((UINT)(x))
#define D(x) ((DAG)(x))


/* platform specific assembly support */
#if defined(_MSC_VER) && defined(_M_X64)
	#define MSCx64
	#include <intrin.h>
	#pragma intrinsic(_umul128)
	#pragma intrinsic(_addcarry_u64)
	#pragma intrinsic(_subborrow_u64)
#elif defined(__GNUC__) && defined(__x86_64__)
	#define GNUx64
#elif defined(__SIZEOF_INT128__)
	#define GNU128
#endif


/* floating point sign */
#define fsign(x) (((x) < 0) ? -1 : ((x) > 0))

/* print bits */
void uprint(UINT x)
{
	int i;
	i = WORDSIZE-1;
	for (; i >= 0; i--) {
		printf("%s", (x >> i) & 1 ? "1" : "0");
	}
	printf("\n");
}


/* hash x into h */
UINT uhash(UINT h, UINT x)
{
#if WORDSIZE==32
	return (h ^ x)*270566475UL;
#else
	return (h ^ x)*36064050381096011ULL;
#endif
}


/* xorshift from Marsaglia */
UINT urandom()
{
	UINT64 ss;
	static UINT32 rs = 2463534242UL;
	rs ^= (rs << 13);
	rs ^= (rs >> 17);
	rs ^= (rs << 5);
	if (WORDSIZE==32) return rs;
	ss = rs;
	rs ^= (rs << 13);
	rs ^= (rs >> 17);
	rs ^= (rs << 5);
	ss = (ss << 32) | rs;
	return (UINT)ss;
}


/* absolute value */
UINT uabs(INT x)
{
	UINT s;
	s = x >> (WORDSIZE-1);
	return (x + s) ^ s;
}
#define ABS(x) uabs(x)


/* unsigned maximum */
UINT umax(UINT a, UINT b)
{
	return a > b ? a : b;
}


/* unsigned minimum */
UINT umin(UINT a, UINT b)
{
	return a > b ? b : a;
}


/* next power of two */
UINT up2(UINT x)
{
	x--;
	x |= x >> 2;
	x |= x >> 1;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
#if WORDSIZE==64
	x |= x >> 32;
#elif WORDSIZE > 64
	#error "port up2"
#endif
	return x+1;
}


/* count leading zeroes */
UINT ulz(UINT x)
{
	UINT n=0;
	if (!x) return WORDSIZE;
#if WORDSIZE==32
	if (x <= 0x0000FFFF) n += 16, x = x << 16;
	if (x <= 0x00FFFFFF) n += 8, x = x << 8;
	if (x <= 0x0FFFFFFF) n += 4, x = x << 4;
	if (x <= 0x3FFFFFFF) n += 2, x = x << 2;
	if (x <= 0x7FFFFFFF) n++;
#elif WORDSIZE==64
	if (x <= 0x00000000FFFFFFFF) n += 32, x = x << 32;
	if (x <= 0x0000FFFFFFFFFFFF) n += 16, x = x << 16;
	if (x <= 0x00FFFFFFFFFFFFFF) n += 8, x = x << 8;
	if (x <= 0x0FFFFFFFFFFFFFFF) n += 4, x = x << 4;
	if (x <= 0x3FFFFFFFFFFFFFFF) n += 2, x = x << 2;
	if (x <= 0x7FFFFFFFFFFFFFFF) n++;
#else
	#error port ulz
#endif
	return n;
}


/* count trailing zeroes */
UINT utz(UINT x)
{
	UINT n;
	for (n=0; !(x & 1); n++) x = x >> 1;
	return n;
}


/* popcount */
UINT upop(UINT n)
{
	UINT i;
	for (i=0; n; i++) n &= n-1;
	return i;
}


/* floor of log2 for n > 0 */
UINT ulog2(UINT n)
{
	return WORDSIZE-1-ulz(n);
}


/* modified Euclid */
UINT ugcd(UINT a, UINT b)
{
	UINT c;
	while (a && b) {
		if (a < b) c = a, a = b, b = c;
		c = a % b;
		a = b - c;
		b = c;
	}
	return (a | b);
}


/* a^n binary powering */
UINT upow(UINT a, UINT n)
{
	UINT r, s;
	r = 1, s = a;
	if (n & 1) r = s;
	n = n >> 1;
	while (n) {
		s = s * s;
		/* from LSB to MSB */
		if (n & 1) r = r * s;
		n = n >> 1;
	}
	return r;
}

/* 10^n lookup */
UINT upow10(UINT n)
{
#if WORDSIZE==32
	static UINT table[10] = {0};
#elif WORDSIZE==64
	static UINT table[20] = {0};
#else
	#error port upow10
#endif
	if (table[n]) return table[n];
	table[n] = upow(10,n);
	return table[n];
}

/* floor of log(a) base b */
UINT ulog(UINT a, UINT b)
{
	UINT x, y, i;
	if (b <= 1) return 0;
	for (x=1, i=0; i < WORDSIZE; i++) {
		y = x*b;
		if (y > a) break;
		if (y/b != x) break;
		x = y;
	}
	return i;
}


/* floor of a^(1/n) */
UINT uroot(UINT a, UINT n)
{
	UINT i, j, k, x, y, e, w;
	w = WORDSIZE;
	if (n == 0) return 0;
	if (n == 1) return a;
	if (a <= 1) return a;
	if (n >= w) return 1;
	/* binary search */
	i = k = 1;
	k = k << w/n;
	while (i <= k) {
		j = (i + k)/2;
		/* x = j to the nth power */
		for (x=j, e=1; e < n; e++) {
			y = x*j;
			/* overflow check */
			if (y/j != x) break;
			x = y;
		}
		/* x too big or too small */
		if (e < n || x > a) k = j-1;
		else if (x < a) i = j+1;
		else return j;
	}
	return k;
}


/* 1/a mod b from Knuth */
UINT uinvmod(UINT a, UINT b)
{
	UINT p, q, r, x, y, i;
	x = 0, y = 1, p = b;
	/* y < 0 if i odd */
	for (i=0; b; i++) {
		q = a / b;
		r = a % b;
		a = b, b = r;
		r = x;
		x = y - q*x;
		y = r;
	}
	/* not invertible */
	if (a > 1) return 0;
	if (i & 1) y += p;
	return y;
}


/* (a0:a1) += (b0:b1) return carry */
UINT uadd2(UINT *a0, UINT *a1, UINT b0, UINT b1)
{
#if defined(MSCx64)
	unsigned char c;
	c = _addcarry_u64(0,*a0,b0,a0);
	c = _addcarry_u64(c,*a1,b1,a1);
	return c;
#else
	UINT c, d, e;
	*a0 = *a0 + b0;
	  c = *a0 < b0;
	*a1 = *a1 + b1;
	  d = *a1 < b1;
	*a1 = *a1 + c;
	  e = *a1 < c;
	return d | e;
#endif
}


/* (a0:a1) -= (b0:b1) return borrow */
UINT usub2(UINT *a0, UINT *a1, UINT b0, UINT b1)
{
#if defined(MSCx64)
	unsigned char c;
	c = _subborrow_u64(0,*a0,b0,a0);
	c = _subborrow_u64(c,*a1,b1,a1);
	return c;
#else
	UINT t1, c, d, e;
	b0 = *a0 - b0;
	 c = *a0 < b0;
	b1 = *a1 - b1;
	 d = *a1 < b1;
	t1 =  b1 - c;
	 e =  b1 < t1;
	*a0 = b0;
	*a1 = t1;
	return d | e;
#endif
}


/* (lo:hi) = a*b */
UINT umul2(UINT a, UINT b, UINT *hi)
{
#if WORDSIZE==32
	UINT64 t0 = (UINT64)(a)*(b);
	*hi = (t0 >> 32);
	return t0;

#elif defined(MSCx64)
	return _umul128(a,b,hi);

#elif defined(GNUx64)
	__asm__("mulq %1" : "=a"(a), "=d"(b) : "0"(a), "1"(b) : "cc");
	*hi = b; return a;

#elif defined(GNU128)
	unsigned __int128 t0 = (unsigned __int128)(a)*(b);
	*hi = t0 >> 64;
	return t0;

#elif WORDSIZE==64
	UINT64 a0, b0, a1, b1, r0, r1, r2, t0;
	a0 = a & 0xFFFFFFFF;
	b0 = b & 0xFFFFFFFF;
	a1 = (a >> 32);
	b1 = (b >> 32);
	t0 = a0 * b0;
	r0 = t0 & 0xFFFFFFFF;
	t0 = a1 * b0 + (t0 >> 32);
	r1 = t0 & 0xFFFFFFFF;
	r2 = (t0 >> 32);
	t0 = a0 * b1 + r1;
	/* assign upper 64-bits to hi */
	*hi = a1 * b1 + (t0 >> 32) + r2;
	return (t0 << 32) + r0;
#else
	#error port umul2
#endif
}


/* (lo:hi)/v = q,r from Hacker's Delight */
UINT udiv2(UINT lo, UINT hi, UINT v, UINT *r)
{
#if WORDSIZE==32
	UINT64 t0 = hi;
	t0 = (t0 << 32) + lo;
	if (r) *r = t0 % v;
	return (UINT)(t0 / v);

#elif defined(GNUx64)
	__asm__("divq %4" : "=a"(lo), "=d"(hi) : "0"(lo), "1"(hi), "r"(v) : "cc");
	if (r) *r = hi; return lo;

#elif defined(GNU128)
	unsigned __int128 t0 = hi;
	t0 = (t0 << 64) | lo;
	if (r) *r = t0 % v;
	return (UINT)(t0 / v);

#elif WORDSIZE==64
	UINT64 b, un1, un0, vn1, vn0, q1, q0, un32, un21, un10, rhat;
	INT64  s=0;
	b = (UINT64)1 << 32;
	if (hi >= v) {	/* overflow */
		printf("udiv2 overflow\n");
		r = 0;
		*((UINT *)r) = 0;
		if (r) *r = 0xFFFFFFFFFFFFFFFF;
		return 0xFFFFFFFFFFFFFFFF;
	}
	s = ulz(v);
	v = v << s;
	vn1 = v >> 32;
	vn0 = v & 0xFFFFFFFF;
	un32 = (hi << s) | ((lo >> (64 - s)) & ((-s) >> 63));
	un10 = (lo << s);
	un1 = un10 >> 32;
	un0 = un10 & 0xFFFFFFFF;
	q1 = un32/vn1;
	rhat = un32 - q1*vn1;
    q1:	if (q1 >= b || q1*vn0 > b*rhat + un1) {
		q1 = q1 - 1;
		rhat = rhat + vn1;
		if (rhat < b) goto q1;
	}
	un21 = un32*b + un1 - q1*v;
	q0 = un21/vn1;
	rhat = un21 - q0*vn1;
    q0:	if (q0 >= b || q0*vn0 > b*rhat + un0) {
		q0 = q0 - 1;
		rhat = rhat + vn1;
		if (rhat < b) goto q0;
	}
	if (r) *r = (un21*b + un0 - q0*v) >> s;
	return q1*b + q0;
#else
	#error port udiv2
#endif
}


/* reciprocal of divisor */
/* it must be shifted up */
UINT nreciprocal(UINT d)
{
	UINT u0, u1, q;
	static UINT D = 0;
	static UINT Q = 0;
	if (d==D) return Q;
	u0 = -1; u1 = -d-1;
	q = udiv2(u0,u1,d,0);
	D = d, Q = q;
	return q;
}


/* 2/1 division from Moller and Granlund 2011 */
/* d is a normalized divisor with top bit set */
/* v = nreciprocal(d) divide u0 + u1*2^B by d */
UINT ndiv2(UINT u0, UINT u1, UINT d, UINT v, UINT *R)
{
	UINT q0, q1, r;
	q0 = umul2(u1,v,&q1);
	q0 += u0;
	q1 += u1 + (q0 < u0) + 1;
	r = u0 - q1*d;
	if (r > q0) q1--, r += d;
	if (r >= d) q1++, r -= d;
	if (R) *R = r;
	return q1;
}


/* multiply a * b mod c */
/* directly divide by c */
UINT umulmod(UINT a, UINT b, UINT c)
{
	/* (a:b) = a * b */
	/* a = (a:b) % c */
	a = umul2(a,b,&b);
	udiv2(a,b,c,&a);
	return a;
}

/* mulmod a * b mod c where */
/* u = c << w is normalized */
/* and v = nreciprocal of u */
/* and both a < c and b < c */
UINT nmulmod(UINT a, UINT b, UINT u, UINT v, UINT w)
{
	/* (a:b) = (a*b) << w */
	a = umul2(a << w, b, &b);
	ndiv2(a,b,u,v,&a);
	return a >> w;
}


/* a^b mod c, binary method */
/* u = c << w is normalized */
/* and v = nreciprocal of u */
UINT npowmod(UINT a, UINT b, UINT u, UINT v, UINT w)
{
	UINT r=1, s;
	s = a % (u >> w);
	if (b & 1) r = s;
	while (b >>= 1) {
		s = nmulmod(s,s,u,v,w);
		if (b & 1) r = nmulmod(r,s,u,v,w);
	}
	return r;
}


/* a^b mod c, binary powering */
UINT upowmod(UINT a, UINT b, UINT c)
{
	UINT u, v, w;
	w = ulz(c), u = c << w, v = nreciprocal(u);
	return npowmod(a,b,u,v,w);
}


/* a+b mod c, a and b are reduced */
UINT uaddmod(UINT a, UINT b, UINT c)
{
	return (a < c-b ? a+b : a+b-c);
}


/* a-b mod c, a and b are reduced */
UINT usubmod(UINT a, UINT b, UINT c)
{
	return (a < b ? a-b+c : a-b);
}


/* Miller-Rabin composite test */
UINT uiscomposite(UINT n, UINT b)
{
	UINT e, i, j, x;
	e = n - 1;
	/* find (j,e) with n-1 = 2^j*e */
	for (j=0; !(e & 1); j++) e = e/2;
	x = upowmod(b,e,n);
	if (x==1) return 0;
	if (x==n-1) return 0;
	for (i=0; i < j; i++) {
		x = umulmod(x,x,n);
		if (x==n-1) return 0;
	}
	return 1;
}


/* primality test */
UINT uprime(UINT n)
{
	/* the first 12 prime bases suffice for 64-bit n */
	UINT prime[12] = {2,3,5,7,11,13,17,19,23,29,31,37};
	UINT i, j;
	if (n < 2) return 0;
	for (i=0; i < 12; i++) {
		j = prime[i];
		if (n < j*j) return 1;
		if ((n % j)==0) return (n==j);
	}
	/* Miller-Rabin */
	/* OEIS A014233 */
	j = (n < 2047) ? 1
	: (n < 1373653) ? 2
	: (n < 25326001) ? 3
	: (n < 3215031751UL) ? 4
	: 12;	/* 64-bits */
	for (i=0; i < j; i++) {
		if (uiscomposite(n,prime[i])) break;
	}
	return (i==j);
}


/* Pollard Rho method */
UINT ufactor_rho(UINT n)
{
	UINT c, d, g, x, y, z, i, j, u, v, w;
	w = ulz(n), u = n << w, v = nreciprocal(u);
	x = urandom() % n;
	c = urandom() % n;
	y = x, z = 1, g = 1;
	/* bound: O(n^(1/4)) */
	j = U(4) << ulog2(n)/4;
	for (i=0; i < j; i++) {
		/* x = x*x + c mod n */
		/* z = z*(x-y) mod n */
		x = nmulmod(x,x,u,v,w);
		x = uaddmod(x,c,n);
		d = usubmod(x,y,n);
		z = nmulmod(z,d,u,v,w);
		/* gcd at powers of 2 */
		if (i & (i+1)) continue;
		g = ugcd(n,z);
		if (g > 1) break;
		y = x, z = 1;
	}
	return g;
}


/* put prime factorization */
/* into x[0..k-1] return k */
UINT ufactor(UINT n, UINT *x)
{
	/* trial division for all primes up to 5 bits */
	UINT prime[11] = {2,3,5,7,11,13,17,19,23,29,31};
	UINT i, j, k, p;
	if (n < 2) return 0;
	/* trial division */
	for (i=k=0; i < 11; i++) {
		p = prime[i];
		if (n % p) continue;
		x[k++] = p; n = n/p; i--;
	}
	if (n==1) return k;
	/* split all composite factors */
	for (i=k, x[k++]=n; i < k; i++) {
		n = x[i], p = 1;
		if (uprime(n)) continue;
		/* loop: Pollard Rho is probabilistic */
		while (p==1 || p==n) p = ufactor_rho(n);
		x[k++] = n/p;
		x[i--] = p;
	}
	/* selection sort */
	for (i=0; i < k; i++) {
		for (j=i+1; j < k; j++) {
			if (x[i] <= x[j]) continue;
			p = x[i], x[i] = x[j], x[j] = p;
		}
	}
	return k;
}


/* Chinese remainder theorem */
/* compute x = r[i] mod p[i] */
UINT uchinese(UINT *p, UINT *r, UINT n)
{
	UINT i, m, x, y, z;
	m = 1, x = 0;
	for (i=0; i < n; i++) {
		z = uinvmod(m,p[i]);
		y = usubmod(r[i], x % p[i], p[i]);
		x = x + m * umulmod(y,z,p[i]);
		m = m * p[i];
	}
	return x;
}


/* Euler's function */
UINT utotient(UINT n)
{
	UINT x[WORDSIZE];
	UINT p, r, i, k;
	p = 0, r = 1;
	k = ufactor(n,x);
	for (i=0; i < k; i++) {
		if (p==x[i]) r *= p;
		else p=x[i], r *= p-1;
	}
	return r;
}


/* order of a modulo n */
UINT uorder(UINT a, UINT n)
{
	UINT x[WORDSIZE];
	UINT i, j, k, s, t;
	if (ugcd(a,n) > 1) return 0;
	t = utotient(n);
	k = ufactor(t,x);
	/* remove factors of t */
	for (i=j=0; i < k; i++) {
		if (x[i]==j) continue;
		s = t / x[i];
		if (upowmod(a,s,n)==1) t = s;
		else j = x[i];	/* skip j */
	}
	return t;
}


/* Jacobi symbol (a|n) */
INT ujacobi(UINT a, UINT n)
{
	INT r, t=1;
	if (n < 2) return 0;
	while (a %= n) {
		while ((a & 1)==0) {
			a = a / 2;
			r = n & 7;
			/* n = 3 or 5 modulo 8 */
			if (r==3 || r==5) t = -t;
		}
		/* both a and n are 3 modulo 4 */
		if ((a & 3)==3 && (n & 3)==3) t = -t;
		r = a, a = n, n = r;
	}
	return (n==1) ? t : 0;
}


/* solve a = x^b mod c by multiplying */
UINT urootmod(UINT a, UINT b, UINT c)
{
	UINT u, v, w, r, i;
	if (a >= c) a = a % c;
	w = ulz(c), u = c << w;
	v = nreciprocal(u);
	for (i=0; i < c; i++) {
		r = npowmod(i,b,u,v,w);
		if (r==a) return i;
	}
	return -1;
}


/* solve a = b^x mod c by multiplying */
UINT ulogmod_mul(UINT a, UINT b, UINT c)
{
	UINT r, s, u, v, w, i;
	if (a >= c) a = a % c;
	if (b >= c) b = b % c;
	w = ulz(c), u = c << w;
	v = nreciprocal(u);
	r = 1, s = 0;
	for (i=0; r != s; i++) {
		if (r==a) return i;
		/* detect cycle to halt */
		if ((i & (i+1))==0) s = r;
		r = nmulmod(r,b,u,v,w);
	}
	return -1;
}


/* solve a = b^x mod c by Pollard Rho */
UINT ulogmod_rho(UINT a, UINT b, UINT c)
{
	UINT i, g, n, u, v, w, x, y, z, X, Y, Z;
	if (a >= c) a %= c;
	if (b >= c) b %= c;
	if (a==1) return 0;
	if (a==b) return 1;
	n = uorder(b,c);
	/* <b> is cyclic group of order n */
	if (n==0) return ulogmod_mul(a,b,c);
	if (upowmod(a,n,c) != 1) return -1;

	/* choose {x,y,z} with x = a^y*b^z mod c */
	w = ulz(c), u = c << w, v = nreciprocal(u);
	y = urandom() % n; y = y/ugcd(y,n);
	z = urandom() % n; z = z/ugcd(z,n);
	Y = npowmod(a,y,u,v,w);
	Z = npowmod(b,z,u,v,w);
	x = nmulmod(Y,Z,u,v,w);
	X = x, Y = y, Z = z;
	for (i=0; 1; i++) {
		switch (x & 7) {
		case 0: case 1: case 2:
			x = nmulmod(a,x,u,v,w);
			y = uaddmod(y,1,n);
			break;
		case 3: case 4: case 5:
			x = nmulmod(b,x,u,v,w);
			z = uaddmod(z,1,n);
			break;
		case 6: case 7:
			x = nmulmod(x,x,u,v,w);
			y = uaddmod(y,y,n);
			z = uaddmod(z,z,n);
			break;
		}
		if (x==X) break;
		if (i & (i+1)) continue;
		X = x, Y = y, Z = z;
	}

	/* we now have a^y*b^z = a^Y*b^Z mod c */
	/* if there exists x for a = b^x mod c */
	/* then b^(x*y+z) = b^(x*Y+Z) modulo c */
	/* x*(y-Y) = (Z-z) mod n, n = order(b) */
	/* cancel g = gcd(y-Y,Z-z,n) and solve */
	/* for x mod (n/g), then try x+i*(n/g) */

	y = usubmod(y,Y,n);
	z = usubmod(Z,z,n);
	g = ugcd(y,ugcd(z,n));
	y = y/g, z = z/g, n = n/g;
	x = uinvmod(y,n);
	if (!x) return -1;
	x = umulmod(x,z,n);

	/* test every possible x+i*(n/g) mod c */
	/* detect cycles in b^x mod c as we go */
	y = npowmod(b,x,u,v,w);
	z = npowmod(b,n,u,v,w);
	for (Y=-1, i=0; y != Y; i++) {
		if (y==a) return x;
		if ((i & (i+1))==0) Y = y;
		/* x += n, y *= z mod c */
		y = nmulmod(y,z,u,v,w);
		x = uaddmod(x,n,c);
	}
	return -1;
}


/* solve a = b^x mod c by factoring n */		
UINT ulogmod_fac(UINT a, UINT b, UINT c)
{
	UINT x[WORDSIZE];
	UINT d, e, f, g, h, i, j, k, m, n, p, q, r, s, t, z;
	UINT (*ulogmod)(UINT a, UINT b, UINT c);
	if (a >= c) a %= c;
	if (b >= c) b %= c;
	if (a==1) return 0;
	if (a==b) return 1;
	n = uorder(b,c);
	/* <b> is cyclic group of order n */
	if (n==0) return ulogmod_mul(a,b,c);
	if (upowmod(a,n,c) != 1) return -1;

	/* the Pohlig-Hellman algorithm */
	/* factor n = product pi^ei*... */
	/* raise base to the power n/pi */
	/* to get element with order pi */
	/* lift to a solution mod pi^ei */
	/* chinese remainder the result */
	s = 0; m = 1;
	d = uinvmod(b,c);
	k = ufactor(n,x);
	for (i=0; i < k; i=j) {
		p = x[i];
		g = upowmod(b,n/p,c);
		z = a, e = n, q = 1, r = 0;
		/* g has order p : Pollard Rho for large p */
		ulogmod = p >> 9 ? ulogmod_rho : ulogmod_mul;
		for (j=i; j < k; j++) {
			if (p != x[j]) break;
			e = e / p;
			f = upowmod(z,e,c);
			h = ulogmod(f,g,c);
			if (h == -1) return -1;
			h *= q, r += h, q *= p;
			/* lifting step needed */
			if (j < k-1 && p==x[j+1]) {
				z = umulmod(z,upowmod(d,h,c),c);
			}
		}
		/* Chinese remainder theorem */
		/* s = s + m*((r-s)/m mod q) */
		t = usubmod(r, s % q, q);
		t = umulmod(uinvmod(m,q),t,q);
		s += m*t, m *= q;
	}
	return s;
}


/* solve a = b^x mod c */
UINT ulogmod(UINT a, UINT b, UINT c)
{
	UINT i, j;
	if (a >= c) a %= c;
	if (b >= c) b %= c;
	if (a==1) return 0;
	if (a==b) return 1;
	if (c >> 9) {
		/* gcds detect no solution */
		i = ugcd(a,c), j = ugcd(b,c);
		if (i==1 && j==1) return ulogmod_fac(a,b,c);
		if (i==1 || j==1) return -1;
		if ((i % j) != 0) return -1;
	}
	return ulogmod_mul(a,b,c);
}


/* nth prime cache */
UINT *uthprime_array = 0;
UINT  uthprime_bound = 0;

/* nth prime */
UINT uthprime_cache(UINT n)
{
	UINT *P, i, j;
	P = uthprime_array;
	i = uthprime_bound;
	if (n==0) return 0;
	if (n-1 < i) return P[n-1];
	P = realloc(P, n*sizeof(UINT));
	j = (i==0) ? 2 : P[--i];
	for (; i < n; j++) {
		if (uprime(j)) P[i++] = j;
	}
	uthprime_array = P;
	uthprime_bound = i;
	return j-1;
}


/* Legendre's func. */
INT uphi(INT x, INT a)
{
	if (x == 0) return 0;
	if (x <= a) return 1;
	if (a == 1) return (x+1)/2;
	return uphi(x,a-1) - uphi(x/uthprime_cache(a), a-1);
}


/* count primes <= x */
UINT upi(UINT x)
{
	UINT s, b, c, p, i;
	/* Meissel's method */
	if (x < 256) goto count;
	b = upi(uroot(x,2));
	c = upi(uroot(x,3));
	p = uthprime_cache(c+1);
	s = uphi(x,c) + (b+c-2)*(b-c+1)/2;
	for (i=c+1; i <= b; i++) {
		s -= upi(x/p); p+=2;
		while (!uprime(p)) p+=2;
	}
	return s;
count:
	s = (x >= 2);
	for (i=3; i <= x; i+=2) {
		if (uprime(i)) s++;
	}
	return s;
}


/* identify primes <= x */
char * uprimesieve(UINT x)
{
	char *s;
	UINT i, j;
	s = malloc(x+1);
	if (!s) return 0;
	for (i=0; i <= x; i++) s[i] = 1;
	s[0] = s[1] = 0;
	/* sieve of Eratosthenes */
	for (i=2; i*i <= x; i++) {
		if (!s[i]) continue;
		for (j=i*i; j <= x; j+=i) s[j] = 0;
	}
	return s;
}


/* Extended Euclidean algorithm */
/* returns gcd(a,b) = a*x + b*y */
/* note the integers are signed */
INT igcdext(INT a, INT b, INT *x, INT *y)
{
	INT q, r, s, t, u, v;
	s = 1, t = 0;
	u = 0, v = 1;
	while (b != 0) {
		q = a/b;
		r = a - q*b; a = b, b = r;
		r = s - q*t; s = t, t = r;
		r = u - q*v; u = v, v = r;
	}
	if (a < 0) {
		a = -a;
		s = -s;
		u = -u;
	}
	if (x) *x = s;
	if (y) *y = u;
	return a;
}


/* Wang's rational reconstruction */
INT iratrec(INT a, INT b, INT *x, INT *y)
{
	INT q, r, s, t, u, v;
	s = b; u = 0;
	t = a; v = 1;
	while (t*t > b/2) {
		q = s/t;
		r = s - q*t; s = t, t = r;
		r = u - q*v; u = v, v = r;
	}
	/* t/v is the result */
	if (v < 0) t = -t, v = -v;
	if (v*v <= b/2 && ugcd(uabs(t),v)==1) {
		*x = t, *y = v;
		return 1;
	}
	else return 0;
}


/*	Fast Fourier Transform

	FFTs of size 2^n require 2^n to divide p-1.  We use
	2^27 * 3 * 5 + 1	30.906 bits, for 32-bit signed integers
	2^56 * 3 * 29 + 1	62.443 bits, for 64-bit signed integers
*/
#if WORDSIZE==32
#define prime1 2013265921
#define recip1 286331150
#define shift1 1
#elif WORDSIZE==64
#define prime1 6269010681299730433
#define recip1 8693293184161972596
#define shift1 1
#else
#error port fft
#endif

INT fft_prime()
{
	return prime1;
}

/* a * b mod p, a and b are reduced */
INT fft_mulmod(INT a, INT b)
{
	UINT u0, u1, p, r, s;
	p = prime1;
	r = recip1;	
	s = shift1;
	a = a << s;
	u0 = umul2(a, b, &u1);
	ndiv2(u0, u1, p << s, r, &r);
	return r >> s;
}

/* a^n mod p, a is reduced, n > 0 */
INT fft_powmod(INT a, INT n)
{
	UINT r, s;
	if (!a) return 0;
	r = 1;
	s = a;
	if (n & 1) r = a;
	n = n >> 1;
	while (n) {
		s = fft_mulmod(s, s);
		if (n & 1) r = fft_mulmod(r, s);
		n = n >> 1;
	}
	return r;
}

/* 1/a mod p, a is reduced */
INT fft_invmod(INT a)
{
	INT b, s, t, q;
	INT oldb, oldt;
	b = a, a = prime1;
	s = 0, t = 1;
	while (b != 0) {
		oldb = b;
		oldt = t;
		q = a / b;
		b = a % b;
		t = s - q * t;
		a = oldb;
		s = oldt;
	}
	s += (s >> (WORDSIZE-1)) & prime1;
	return s;
}

/* nth root of unity */
INT fft_omega(INT n)
{
	INT p, q, r, i;
	p = prime1;
	q = (p - 1) / n;
	r = (p - 1) % n;
	if (r) return 0;
	/* n divides p-1 */
	for (i=2; i < p; i++) {
		r = fft_powmod(i, (p-1)/2);
		if (r == p-1) return fft_powmod(i, q);
	}
	return 0;
}

/* forward transform */
void fft_forward(INT *a, INT *w, INT n)
{
	INT t0, t1, i;
	n = n / 2;
	for (i=0; i < n; i++) {
		t0 = a[i+0] - a[n+i];
		t1 = a[i+0] + a[n+i] - prime1;
		t0 += (t0 >> (WORDSIZE-1)) & prime1;
		t1 += (t1 >> (WORDSIZE-1)) & prime1;
		t0 = fft_mulmod(t0, w[i]);
		a[i+0] = t1;
		a[i+n] = t0;
	}
	if (n < 2) return;
	fft_forward(a+0, w+n, n);
	fft_forward(a+n, w+n, n);
}

/* reverse transform */
void fft_reverse(INT *a, INT *w, INT n)
{
	INT t0, t1, i;
	if (n < 2) return;
	n = n / 2;
	fft_reverse(a+0, w+n, n);
	fft_reverse(a+n, w+n, n);
	for (i=0; i < n; i++) {
		t1 = fft_mulmod(a[n+i], w[i]);
		t0 = a[i+0] - t1;
		t1 = a[i+0] + t1 - prime1;
		t0 += (t0 >> (WORDSIZE-1)) & prime1;
		t1 += (t1 >> (WORDSIZE-1)) & prime1;
		a[i+n] = t0;
		a[i+0] = t1;
	}
}

/* powers of a root of unity */
/* for fft_forward transform */
void fft_warray(INT *w, INT n)
{
	INT i, r;
	r = fft_omega(n);
	w[0] = 1; n = n/2;
	/* compute initial set of powers needed for fft */
	for (i=1; i < n; i++) w[i] = fft_mulmod(r,w[i-1]);
	while (n > 1) {
		w = w + n;
		n = n / 2;
		/* append powers for smaller transforms */
		/* this reduces cache misses in the fft */
		for (i=0; i < n; i++) w[i] = w[2*i - 2*n];
	}
}

/* invert the roots of unity */
/* for fft_reverse transform */
void fft_invert(INT *w, INT n)
{
	INT i, t;
	while (n > 2) {
		n = n / 2;
		for (i=1; i < n/2; i++) {
			t = w[i];
			w[i+0] = prime1 - w[n-i];
			w[n-i] = prime1 - t;
		}
		w[i] = prime1 - w[i];
		w = w + n;
	}
}

#if 0
static void fft_test()
{
	UINT p, s, r;
	INT n, v, i;
	INT *a, *w;
	#if WORDSIZE==32
		p = 2013265921;	s = ulz(p);
		r = nreciprocal(p << s);
	#elif WORDSIZE==64
        	p = 6269010681299730433; s = ulz(p);
		r = nreciprocal(p << s);
	#endif
	printf("prime1: %lld\n", (long long int)p);
	printf("recip1: %lld\n", (long long int)r);
	printf("shift1: %lld\n", (long long int)s);

	/* check computations */
	if (p != prime1) return;
	if (r != recip1) return;
	if (s != shift1) return;
	n = 16; v = fft_invmod(n);
	a = malloc(2*n*sizeof(INT)); w = a+n;
	for (i=0; i < n; i++) a[i] = i;
	fft_warray(w,n);
	fft_forward(a,w,n);
	fft_invert(w,n);
	fft_reverse(a,w,n);
	for (i=0; i < n; i++) a[i] = fft_mulmod(a[i],v);
	for (i=0; i < n; i++) printf("%lld ",(INT64)a[i]);
	printf("\n");
	free(a);
}
#endif

#undef prime1
#undef recip1
#undef shift1
