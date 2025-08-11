/*	The F4 Groebner Basis Engine from axcas.net
	by Roman Pearce, August 2025

	This code is released into the public domain.

	This software and documentation is provided "as is", without warranty of any kind,
	express or implied, including but not limited to the warranties of merchantability,
	fitness for a particular purpose, and noninfringement.  In no event shall the authors
	or copyright holders be liable for any claim, damages, or other liability, whether in
	an action of contract, tort, or otherwise, arising from, out of or in connection with
	the software or the use or other dealings in the software.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MULTI_THREAD	1
int num_thread = 0;


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
#elif UINTPTR_MAX==0xFFFFFFFFFFFFFFFF
	#define WORDSIZE 64
	typedef INT64	INT;
	typedef UINT64	UINT;
#else
	#error port WORDSIZE
#endif
#define I(x) ((INT)(x))
#define U(x) ((UINT)(x))


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



/*
	Parallel Programming 
*/

/* Windows 32-bit with MSVC */
#if MULTI_THREAD && defined(_MSC_VER) && !defined(_WIN64)
	#include <time.h>
	#include <intrin.h>

	void __stdcall GetSystemInfo(void *sysinfo);
	UINT __stdcall CreateThread(void *ta, UINT ss, void *f, void *x, int dw, UINT *id);
	 int __stdcall WaitForSingleObject(UINT id, int ms);

	#define cas(ptr,old,new)	_InterlockedCompareExchange(ptr,new,old)

	INT numcpu()
	{
		unsigned int x[20] = {0};
		static INT ncpu = 0;
		if (ncpu) goto done;
		GetSystemInfo(&x);
		ncpu = x[5];
		/* override the number of threads */
		if (num_thread > 0) ncpu = num_thread;
	done:	return ncpu;
	}

	INT thread(void *f, void *x)
	{
		return CreateThread(0,0,f,x,0,0);
	}

	INT join(INT id)
	{
		WaitForSingleObject(id,-1);
		return 0;
	}

	double realtime()
	{
		double t = clock();
		return t / CLOCKS_PER_SEC;
	}

/* Windows 64-bit with MSVC */
#elif MULTI_THREAD && defined(_MSC_VER) && defined(_WIN64)
	#include <time.h>
	#include <intrin.h>

	void __stdcall GetSystemInfo(void *sysinfo);
	UINT __stdcall CreateThread(void *ta, UINT ss, void *f, void *x, int dw, UINT *id);
	 int __stdcall WaitForSingleObject(UINT id, int ms);

	#define cas(ptr,old,new)	_InterlockedCompareExchange64(ptr,new,old)

	INT numcpu()
	{
		unsigned int x[20] = {0};
		static INT ncpu = 0;
		if (ncpu) goto done;
		GetSystemInfo(&x);
		ncpu = x[8];
		/* override the number of threads */
		if (num_thread > 0) ncpu = num_thread;
	done:	return ncpu;
	}

	INT thread(void *f, void *x)
	{
		return CreateThread(0,0,f,x,0,0);
	}

	INT join(INT id)
	{
		WaitForSingleObject(id,-1);
		return 0;
	}

	double realtime()
	{
		double t = clock();
		return t / CLOCKS_PER_SEC;
	}

/* Windows 64-bit with GCC */
#elif MULTI_THREAD && defined(__GNUC__) && defined(_WIN64)
	#include <time.h>

	void __stdcall GetSystemInfo(void *sysinfo);
	UINT __stdcall CreateThread(void *ta, UINT ss, void *f, void *x, int dw, UINT *id);
	 int __stdcall WaitForSingleObject(UINT id, int ms);

	#define cas(ptr,old,new)	__sync_val_compare_and_swap(ptr,old,new)

	INT numcpu()
	{
		unsigned int x[20] = {0};
		static INT ncpu = 0;
		if (ncpu) goto done;
		GetSystemInfo(&x);
		ncpu = x[8];
		/* override the number of threads */
		if (num_thread > 0) ncpu = num_thread;
	done:	return ncpu;
	}

	INT thread(void *f, void *x)
	{
		return CreateThread(0,0,f,x,0,0);
	}

	INT join(INT id)
	{
		WaitForSingleObject(id,-1);
		return 0;
	}

	double realtime()
	{
		double t = clock();
		return t / CLOCKS_PER_SEC;
	}

/* Mac OS */
#elif MULTI_THREAD && defined(__APPLE__)
	#include <sys/time.h>
	#include <unistd.h>
	#include <pthread.h>

	#define cas(ptr,old,new)	__sync_val_compare_and_swap(ptr,old,new)

	INT numcpu()
	{
		static INT ncpu = 0;
		if (ncpu) goto done;
		ncpu = sysconf(_SC_NPROCESSORS_CONF);
		/* override the number of threads */
		if (num_thread > 0) ncpu = num_thread;
	done:	return ncpu;
	}

	INT thread(void *f, void *x)
	{
		int err;
		pthread_t pid;
		err = pthread_create(&pid,0,f,x);
		return (INT)pid;
	}

	INT join(INT id)
	{
		INT val;
		int err;
		pthread_t pid = (pthread_t)id;
		err = pthread_join(pid, (void **)(&val));
		return val;
	}

	double realtime()
	{
		struct timeval t;
		gettimeofday(&t,0);
		return 1.0 * t.tv_sec + 0.000001 * t.tv_usec;
	}

/* Linux */
#elif MULTI_THREAD && defined(__linux__)
	#include <sys/time.h>
	#include <unistd.h>
	#include <pthread.h>

	#define cas(ptr,old,new)	__sync_val_compare_and_swap(ptr,old,new)

	INT numcpu()
	{
		static INT ncpu = 0;
		if (ncpu) goto done;
		ncpu = sysconf(_SC_NPROCESSORS_CONF);
		/* override the number of threads */
		if (num_thread > 0) ncpu = num_thread;
	done:	return ncpu;
	}

	INT thread(void *f, void *x)
	{
		int err;
		pthread_t pid;
		err = pthread_create(&pid,0,f,x);
		return (INT)pid;
	}

	INT join(INT id)
	{
		INT val;
		int err;
		pthread_t pid = (pthread_t)id;
		err = pthread_join(pid, (void **)(&val));
		return val;
	}

	double realtime()
	{
		struct timeval t;
		gettimeofday(&t,0);
		return 1.0 * t.tv_sec + 0.000001 * t.tv_usec;
	}

/* sequential */
#else
	#include <time.h>

	INT cas(INT *ptr, INT old, INT new)
	{
		INT tmp = *ptr;
		if (tmp==old) *ptr = new;
		return tmp;
	}

	INT numcpu()
	{
		return 1;
	}

	INT thread(void *f, void *x)
	{
		void (*g)(void *) = f;
		g(x);
		return 0;
	}

	INT join(INT id)
	{
		return 0;
	}

	double realtime()
	{
		double t = clock();
		return t / CLOCKS_PER_SEC;
	}
	
#endif


/*	F4 Algorithm mod p

	Monomials are integers referring to exponent vectors stored in the f4_monom block of memory.
	EXPON(m) returns the start of the exponent vector and COLUM(m) returns the column number.
	EXPON(0) is scratch space.  Monomials are added to a hash table to make them unique.

	Polynomials have an array of monomials and a array of coefficients plus other information.

	Matrix rows have an array of coefficients (not duplicated) and an array of indices encoded as follows:
	it first records an index using an entire word, this is followed by a sequence of unsigned characters
	recording differences from 1 to 255 from the previous index, ending in the zero character.
*/

typedef struct f4row {
	INT   len;	/* number of terms in the matrix row */
	INT   sug;	/* phantom sugar degree for the poly */
	INT   fac;	/* monomial cofactor from the syzygy */
	INT  *cof;	/* an array of coefficients modulo p */
	INT  *mon;	/* array of monomials from the basis */
	CHAR *ind;	/* array of column index differences */
	INT   siz;	/* bytes of encoded sparsity pattern */
} f4row;

typedef struct f4syz {
	INT   lcm;	/* the lead monomial to be cancelled */
	INT   sug;	/* phantom sugar degree for the pair */
	f4row *row0;	/* pointers to the basis polynomials */
	f4row *row1;
} f4syz;

INT	f4_prime;		/* the current prime */
INT	f4_nvars;		/* total # variables */
INT	f4_nelim;		/* # for elimination */

INT    *f4_table, f4_tsize;	/* hash table & size */
INT    *f4_monom, f4_mload;	/* array & num. used */
INT	f4_mutex;		/* hash & array lock */

f4row **f4_basis;		/* basis polynomials */
INT     f4_bload, f4_bsize;	/* num. used / total */

f4syz **f4_pairs;		/* array of syzygies */
INT     f4_pload, f4_psize;	/* num. used / total */

f4row **f4_extra;		/* extra polynomials */
INT	f4_eload, f4_esize;	/* num. used / total */

f4row **f4_array;		/* big sparse matrix */
INT	f4_aload, f4_asize;	/* rows used / total */

INT    *f4_mused;		/* monomials present */
INT	f4_uload, f4_usize;	/* cols used / total */


/* monomial exponent, column index */
#define EXPON(m) (f4_monom + m*(f4_nvars+1))
#define COLUM(m) (f4_monom + m*(f4_nvars+1))[f4_nvars]
#define HASHVAL  (I(1) << (WORDSIZE-1))

/* normalize x with -p <= x < p to 0 <= x < p */
#define NORMAL(x,p) x += (x >> (WORDSIZE-1)) & p


/*
	monomial hash table
*/

/* initialize */
void f4mon_init(INT nvars, INT nelim)
{
	INT m, s, i, j;

	if (nelim > nvars) nelim = nvars;
	if (nelim < 0) nelim = 0;

	f4_nvars = nvars;
	f4_nelim = nelim;

	/* monomial length */
	m = f4_nvars+1;

	/* must be power of two */
	f4_tsize = s = I(1) << 15;
	f4_table = malloc(2*s*sizeof(INT));
	f4_monom = malloc(m*s/2*sizeof(INT));
	f4_mload = 1;	/* EXPON(0) == scratch */
	for (i=0; i < 2*s  ; i++) f4_table[i] = 0;
	for (i=0; i < m*s/2; i++) f4_monom[i] = 0;
	f4_mutex = 0;
}

/* dispose */
void f4mon_free()
{
	free(f4_table); free(f4_monom);
	f4_table = f4_monom = 0;
	f4_tsize = f4_mload = 0;
	f4_nvars = 0;
	f4_mutex = 0;
}

/* reload monoms */
void f4mon_rehash()
{
	INT *e, n, s, h, i, k, m;
	n = f4_nvars; s = f4_tsize;
	for (i=0; i < 2*s; i++) f4_table[i] = 0;
	for (m=1; m < f4_mload; m++) {
		e = EXPON(m);
		/* hash can't be negative */
		for (h=1,i=0; i < n; i++) h = (uhash(h,e[i]) | HASHVAL) >> 1;
		for (k=h,i=0; i < s; i++) {
			k = (k+i) & (s-1);
			if (f4_table[k] == 0) break;
		}
		f4_table[k+0] = h;
		f4_table[k+s] = m;
	}
}

/* enlarge table */
void f4mon_resize()
{
	INT m, s;
	s = f4_tsize;
	m = f4_nvars+1;
	s = f4_tsize = 2*s;
	f4_table = realloc(f4_table, 2*s*sizeof(INT));
	f4_monom = realloc(f4_monom, m*s/2*sizeof(INT));
	f4mon_rehash();
}


/* 
	monomial operations
*/

/* create monomial */
INT f4mon_new(INT *e)
{
	INT *T, *v, n, s, h, i, j, k, m, t;
	T = f4_table; s = f4_tsize; n = f4_nvars;
	for (h=1,i=0; i < n; i++) h = (uhash(h,e[i]) | HASHVAL) >> 1;
	for (k=h,i=0; i < s; i++) {
		k = (k+i) & (s-1);
		if (T[k] == 0) break;
		if (T[k] != h) continue;
		m = T[s+k]; v = EXPON(m);
		/* check that the exponents equal */
		for (j=0; j < n && e[j]==v[j]; j++) ;
		if (j==n) return m;
	}
	m = f4_mload++; 
	T[s+k] = m; T[k] = h; v = EXPON(m);
	for (j=0; j < n; j++) v[j] = e[j]; v[j] = 0;
	if (m+1 == s/2) f4mon_resize();
	return m;
}

/* new monomial in thread */
INT f4mon_new_thread(INT *e)
{
	INT *T, *v, n, s, h, i, j, k, m, z;
	T = f4_table; s = f4_tsize; n = f4_nvars;
	for (h=1,i=0; i < n; i++) h = (uhash(h,e[i]) | HASHVAL) >> 1;
retry:	for (k=h,i=0; i < s; i++) {
		k = (k+i) & (s-1);
		if (T[k] == 0) break;
		if (T[k] != h) continue;
		m = T[s+k]; v = EXPON(m);
		/* check that the exponents equal */
		for (j=0; j < n && e[j]==v[j]; j++) ;
		if (j==n) return m;
	}
	if (cas(&f4_mutex,0,1)) goto retry;
	for (k=h,i=0; i < s; i++) {
		k = (k+i) & (s-1);
		if (T[k] == 0) break;
		if (T[k] != h) continue;
		m = T[s+k]; v = EXPON(m);
		/* check that the exponents equal */
		for (j=0; j < n && e[j]==v[j]; j++) ;
		if (j==n) goto done;
	}
	m = f4_mload++; v = EXPON(m);
	for (j=0; j < n; j++) v[j] = e[j]; v[j] = 0;
	T[s+k] = m; cas(&T[k],0,h); 
done:	f4_mutex = 0;
	return m;
}

/* compare two monomials */
INT f4mon_cmp(INT A, INT B)
{
	INT da, db;
	INT *a, *b, n, e, i, j, k;
	if (A==B) return 0;
	a = EXPON(A);
	b = EXPON(B);
	n = f4_nvars;
	e = f4_nelim;

	if (e >= n) goto plex;
	i = j = 0; k = e;
	if (e) goto comp;
last:	i = j = e; k = n;

comp:	/* grevlex order */
	for (da=db=0; i < k; i++) {
		da += a[i];
		db += b[i];
	}
	if (da != db) return (da-db);
	while (--i > j) {
		if (a[i] != b[i]) return (b[i]-a[i]);
	}
	if (k==n) return 0;
	goto last;

plex:	for (i=0; i < n; i++) {
		if (a[i] != b[i]) return (a[i] - b[i]);
	}
	return 0;
}

/* total degree */
INT f4mon_deg(INT A)
{
	INT *a, n, i, d;
	n = f4_nvars;
	a = EXPON(A);
	for (d=i=0; i < n; i++) d += a[i];
	return d;
}

/* sugar of quotient */
INT f4mon_sug(INT A, INT B)
{
	INT *a, *b, n, i, d;
	n = f4_nvars;
	a = EXPON(A);
	b = EXPON(B);
	for (d=i=0; i < n; i++) d += a[i] - b[i];
	return d;
}

/* division test */
INT f4mon_div(INT A, INT B)
{
	INT *a, *b, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B);
	for (i=0; i < n && a[i] >= b[i]; i++) ;
	return (i==n);
}

/* multiply */
INT f4mon_mul(INT A, INT B)
{
	INT *a, *b, *c, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B); c = EXPON(0);
	for (i=0; i < n; i++) c[i] = a[i] + b[i];
	return f4mon_new(c);
}

/* multiply in thread */
INT f4mon_mul_thread(INT A, INT B, INT *tmp)
{
	INT *a, *b, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B);
	for (i=0; i < n; i++) tmp[i] = a[i] + b[i];
	return f4mon_new_thread(tmp);
}

/* quotient */
INT f4mon_quo(INT A, INT B)
{
	INT *a, *b, *c, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B); c = EXPON(0);
	for (i=0; i < n; i++) c[i] = a[i] - b[i];
	return f4mon_new(c);
}

/* quotient in thread */
INT f4mon_quo_thread(INT A, INT B, INT *tmp)
{
	INT *a, *b, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B);
	for (i=0; i < n; i++) tmp[i] = a[i] - b[i];
	return f4mon_new_thread(tmp);
}

/* greatest common divisor */
INT f4mon_gcd(INT A, INT B)
{
	INT *a, *b, *c, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B); c = EXPON(0);
	for (i=0; i < n; i++) c[i] = a[i] > b[i] ? b[i] : a[i];
	return f4mon_new(c);
}

/* least common multiple */
INT f4mon_lcm(INT A, INT B)
{
	INT *a, *b, *c, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B); c = EXPON(0);
	for (i=0; i < n; i++) c[i] = a[i] > b[i] ? a[i] : b[i];
	return f4mon_new(c);
}

/* test relatively prime */
INT f4mon_prm(INT A, INT B)
{
	INT *a, *b, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B);
	for (i=0; i < n; i++) {
		if (a[i] && b[i]) return 0;
	}
	return 1;
}

/* test A depends on B */
INT f4mon_dep(INT A, INT B)
{
	INT *a, *b, n, i;
	n = f4_nvars;
	a = EXPON(A); b = EXPON(B);
	for (i=0; i < n; i++) {
		if (!a[i] && b[i]) return 0;
	}
	return 1;
}

/* variable x[k] */
INT f4mon_var(INT k)
{
	INT *a, n, i;
	n = f4_nvars;
	a = EXPON(0);
	for (i=0; i < n; i++) a[i] = (i==k);
	return f4mon_new(a);
}

/* constant */
INT f4mon_one()
{
	INT *a, n, i;
	n = f4_nvars;
	a = EXPON(0);
	for (i=0; i < n; i++) a[i] = 0;
	return f4mon_new(a);	
}

/* for debugging */
void f4mon_print(INT A)
{
	INT *e, n, i, d;
	n = f4_nvars;
	e = EXPON(A);
	for (d=i=0; i < n; i++) {
		d += e[i];
		if (e[i]==0) continue;
		printf("*x%lld^%lld",(long long int)i,(long long int)e[i]);
	}
	if (d==0) printf("*1");
}



/* sort array of monomials */
void f4mon_sort(INT *L, INT l)
{
	INT i, j, k, m;
	k = l;
	if (k < 2) return;
sort:	k = 5*(k+1)/13;
	for (i=k-1; i < l; i++) {
		m = L[i];
		for (j=i; j >= k; j-=k) {
			if (f4mon_cmp(m,L[j-k]) > 0) break;
			L[j] = L[j-k];
		}
		L[j] = m;
	}
	if (k > 1) goto sort;
}


/*
	polynomial operations
*/

f4row * f4row_new()
{
	f4row *a;
	a = malloc(sizeof(f4row));
	memset((void *)a,0,sizeof(f4row));
	return a;
}

void f4row_free(f4row *a)
{
	if (!a) return;
	/* Only free cof and mon if they are owned (fac == 0) */
	/* When fac != 0, these arrays are shared with other structures */
	if (!a->fac) {
		free(a->cof);
		free(a->mon);
	}
	free(a->ind);
	free(a);
}

/* sort polynomial */
f4row * f4row_sort(f4row *a)
{
	INT i, j, k, l, m, c;
	l = k = a->len;
	if (k < 2) goto done;
sort:	k = 5*(k+1)/13;
	for (i=k-1; i < l; i++) {
		m = a->mon[i];
		c = a->cof[i];
		for (j=i; j >= k; j-=k) {
			if (f4mon_cmp(a->mon[j-k],m) >= 0) break;
			a->mon[j] = a->mon[j-k];
			a->cof[j] = a->cof[j-k];
		}
		a->mon[j] = m;
		a->cof[j] = c;
	}
	if (k > 1) goto sort;
done:	return a;
}

/* print row data */
void f4row_print(f4row *a)
{
	INT i, j, k, l, z;
	if (a->ind) {
		i = j = k = 0;
		while (j < a->len) {
			z = *(INT *)(a->ind+k); k += sizeof(INT);
		print:	printf("+%lld[%lld]",(long long int)(a->cof[j]),(long long int)z);
			l = z, j++;
			if (a->ind[k]) {
				z = l - a->ind[k]; k++;
				goto print;
			}
			else k++;
		}
		printf("\n");
		return ;
	}
	if (a->fac) {
		f4mon_print(a->fac);
		printf(" : ");
	}
	for (i=0; i < a->len; i++) {
		printf("+%lld",(long long int)(a->cof[i]));
		f4mon_print(a->mon[i]);
	}
	printf("\n");
}


/*
	pair operations
*/


/* free memory */
void f4_memory()
{
	f4row *a, *b;
	INT *e, *f, i, j, k, l, m, s=0;

	/* mark used polynomials */
	for (i=0; i < f4_pload; i++) {
		a = f4_pairs[i]->row0;
		b = f4_pairs[i]->row1;
		a->fac = b->fac = -1;
	}
	/* free unmarked polynomials */
	for (i=j=0; i < f4_eload; i++) {
		a = f4_extra[i];
		if (a->fac == -1) {
			f4_extra[j++] = a;
			a->fac = 0;
		}
		else {
			k = a->len;
			s += 2*k*sizeof(INT) + sizeof(f4row);
			f4row_free(a);
		}
	}
	f4_eload = j;
	for (i=0; i < f4_bload; i++) {
		f4_basis[i]->fac = 0;
	}

	/* garbage collect monomials: important */
	for (i=0; i < f4_mload; i++) COLUM(i) = 0;
	for (i=0; i < f4_bload; i++) {
		a = f4_basis[i];
		for (j=0; j < a->len; j++) {
			m = a->mon[j];
			COLUM(m) = 1;
		}
	}
	for (i=0; i < f4_eload; i++) {
		a = f4_extra[i];
		for (j=0; j < a->len; j++) {
			m = a->mon[j];
			COLUM(m) = 1;
		}
	}
	for (i=0; i < f4_pload; i++) {
		m = f4_pairs[i]->lcm;
		COLUM(m) = 1;
	}

	/* record new positions */
	for (i=j=1; i < f4_mload; i++) {
		if (COLUM(i)) COLUM(i) = j++;
	}
	for (i=0; i < f4_bload; i++) {
		a = f4_basis[i];
		for (j=0; j < a->len; j++) {
			m = COLUM(a->mon[j]);
			a->mon[j] = m;
		}
	}
	for (i=0; i < f4_eload; i++) {
		a = f4_extra[i];
		for (j=0; j < a->len; j++) {
			m = COLUM(a->mon[j]);
			a->mon[j] = m;
		}
	}
	for (i=0; i < f4_pload; i++) {
		m = f4_pairs[i]->lcm;
		f4_pairs[i]->lcm = COLUM(m);
	}
	for (i=j=1; i < f4_mload; i++) {
		if (COLUM(i)) {
			e = EXPON(i);
			f = EXPON(j);
			for (k=0; k < f4_nvars; k++) f[k] = e[k];
			COLUM(j) = 0;
			j++;
		}
	}
	f4_mload = j;
	f4mon_rehash();
}

/* put in array */
void f4_addrow(f4row *a)
{
	if (f4_aload + 1 > f4_asize) {
		f4_asize = 3*f4_asize/2;
		f4_array = realloc(f4_array, f4_asize*sizeof(f4row *));
	}
	f4_array[f4_aload] = a;
	f4_aload++;
}

/* put in basis */
void f4_insert(f4row *a)
{
	if (f4_bload + 1 > f4_bsize) {
		f4_bsize = 3*f4_bsize/2;
		f4_basis = realloc(f4_basis, f4_bsize*sizeof(f4row *));
	}
	f4_basis[f4_bload] = a;
	f4_bload++;	
}

/* add poly to basis */
void f4_update(f4row *a)
{
	f4row *b;
	f4syz *s;
	INT i, j, k, l;
	INT d0, d1, d2;

	/* room for new polys */
	if (f4_bload + 1 > f4_bsize) {
		f4_bsize = 3*f4_bsize/2;
		f4_basis = realloc(f4_basis, f4_bsize*sizeof(f4row *));
	}

	/* room for new pairs */
	if (f4_pload + f4_bload > f4_psize) {
		f4_psize = 2*f4_psize + f4_bload;
		f4_pairs = realloc(f4_pairs, f4_psize*sizeof(f4syz *));
	}

	/* room for junk poly */
	if (f4_eload + f4_bload > f4_esize) {
		f4_esize = 2*f4_esize + f4_bload;
		f4_extra = realloc(f4_extra, f4_esize*sizeof(f4row *));
	}

	/* generate new pairs */
	for (i=l=0; i < f4_bload; i++) {
		b = f4_basis[i];
		/* Buchberger's first criteria */
		if (f4mon_prm(a->mon[0],b->mon[0])) continue;
		s = malloc(sizeof(f4syz));
		s->lcm = f4mon_lcm(b->mon[0],a->mon[0]);
		d0 = f4mon_sug(s->lcm, a->mon[0]) + a->sug;
		d1 = f4mon_sug(s->lcm, b->mon[0]) + b->sug;
		s->sug = umax(d0, d1);
		s->row0 = b; s->row1 = a;
		f4_pairs[f4_pload+l] = s;
		l++;	/* # new pairs */
	}

	/* Gebauer & Moller criteria B */
	for (i=0; i < f4_pload; i++) {
		d0 = f4_pairs[i]->lcm;
		d1 = f4_pairs[i]->row0->mon[0];
		d2 = f4_pairs[i]->row1->mon[0];
		if (f4mon_div(d0,a->mon[0])
		 && f4mon_lcm(d1,a->mon[0]) != d0
		 && f4mon_lcm(d2,a->mon[0]) != d0) {
			f4_pairs[i]->lcm = 0;
		}
	}

	/* Gebauer & Moller criteria M */
	for (i=0; i < l; i++) {
		d0 = f4_pairs[f4_pload+i]->lcm;
		if (d0==0) continue;
		for (j=0; j < l; j++) {
			d1 = f4_pairs[f4_pload+j]->lcm;
			if (i==j || d0==d1 || d1==0) continue;
			if (f4mon_div(d1,d0)) f4_pairs[f4_pload+j]->lcm = 0;
		}
	}

	/* Gebauer & Moller criteria F */
	for (i=0; i < l; i++) {
		d0 = f4_pairs[f4_pload+i]->lcm;
		if (d0==0) continue;
		for (j=i+1; j < l; j++) {
			d1 = f4_pairs[f4_pload+j]->lcm;
			if (d1==0) continue;
			if (d0==d1) f4_pairs[f4_pload+j]->lcm = 0;
		}
	}

	/* compact the set of pairs */
	for (i=j=0; i < f4_pload+l; i++) {
		d0 = f4_pairs[i]->lcm;
		if (d0 == 0) {
			free(f4_pairs[i]);
			continue;
		}
		f4_pairs[j] = f4_pairs[i];
		j++;
	}
	f4_pload = j;


	/* cull basis elements */
	for (i=j=0; i < f4_bload; i++) {
		b = f4_basis[i];
		if (!f4mon_div(b->mon[0],a->mon[0])) {
			f4_basis[j++] = b;
		}
		else if (f4_nelim==0 || a->len <= b->len) {
			f4_extra[f4_eload++] = b;
		}
		else {
			f4_basis[j++] = b;
		}
	}
	f4_bload = j;

	/* add new polynomial */
	f4_basis[f4_bload] = a;
	f4_bload++;
}


/* pair degree */
INT f4_degree()
{
	f4syz *s;
	INT d, e, i;
	d = (INT)1 << (WORDSIZE-2);
	for (i=0; i < f4_pload; i++) {
		s = f4_pairs[i];
		e = f4mon_deg(s->lcm);
		if (e < d) d = e;
	}
	return i ? d : -1;
}

/* sugar degree */
INT f4_sugar()
{
	f4syz *s;
	INT d, e, i;
	d = (INT)1 << (WORDSIZE-2);
	for (i=0; i < f4_pload; i++) {
		s = f4_pairs[i];
		e = s->sug;
		if (e < d) d = e;
	}
	return i ? d : -1;
}

/* Static variables that need resetting */
static INT f4_select_dprev = 0;
static INT f4_select_plimit = 2048;

/* Reset static variables - called during cleanup */
void f4_select_reset_statics(void) {
	f4_select_dprev = 0;
	f4_select_plimit = 2048;
}

/* Initialize static variables - called during init */
void f4_select_init_statics(void) {
	f4_select_dprev = 0;
	f4_select_plimit = 2048;
}

/* select pairs */
void f4_select()
{
	f4syz *s, *t;
	f4row **A, *a, *b;
	INT d, e, i, j, k, l, d0, z;

	/* determine degree, # of pairs */
	d = (INT)1 << (WORDSIZE-2); k = 0;
	for (i=0; i < f4_pload; i++) {
		s = f4_pairs[i];
		e = f4_nelim ? s->sug : f4mon_deg(s->lcm);
		if (e > d) continue;
		else if (e == d) k++;
		else d = e, k = 1;
	}
	l = k;

	if (d == f4_select_dprev) {
		f4_select_plimit *= 2;
	}
	else {
		f4_select_plimit = 2048;
		f4_select_dprev = d;
	}
	if (k > f4_select_plimit) k = f4_select_plimit;

        if (getenv("RUR_VERBOSE_F4"))
            printf("degree=%lld pairs=%lld/%lld with pairs limit %lld\n",
                (long long int)d,
                (long long int)k, 
                (long long int)f4_pload,
                (long long int)f4_select_plimit
            );

	/* allocate matrix */
	if (2*k > f4_asize) {
		f4_array = realloc(f4_array, 2*k*sizeof(f4row *));
		f4_asize = 2*k;
	}
	f4_aload = 0;

	/* seed matrix using pairs */
	for (i=j=0; i < f4_pload; i++) {
		s = f4_pairs[i];
		e = f4_nelim ? s->sug : f4mon_deg(s->lcm);
		if (e > d) continue;
		if (j >= 2*k) continue;

		a = malloc(sizeof(f4row));
		b = malloc(sizeof(f4row));
		a->len = s->row0->len;
		b->len = s->row1->len;
		a->fac = f4mon_quo(s->lcm, s->row0->mon[0]);
		b->fac = f4mon_quo(s->lcm, s->row1->mon[0]);
		a->cof = s->row0->cof;
		b->cof = s->row1->cof;
		a->mon = s->row0->mon;
		b->mon = s->row1->mon;
		a->ind = b->ind = 0;
		f4_array[j++] = a;
		f4_array[j++] = b;

		s->lcm = 0;
	}
	f4_aload = j;

	/* remove selected pairs */
	for (i=j=0; i < f4_pload; i++) {
		s = f4_pairs[i];
		if (s->lcm == 0) free(s);
		else f4_pairs[j++] = s;
	}
	f4_pload = j;
}

#if MULTI_THREAD

void f4_symbol_thread(void *par)
{
	f4row *a, *b, *c;
	INT *param = par;
	INT i, j, k, l, m, q, *tmp;
	INT ss, st, sz;
	INT y, z;

	tmp = 0;
	sz = param[1];
	st = param[0];
start:	ss = cas(param,st,st+1);
	if (ss != st) {
		st = ss;
		goto start;
	}
	if (st >= sz) goto done;

	if (!tmp) tmp = malloc(f4_nvars*sizeof(INT));

	a = f4_array[st];
	for (i=0; i < a->len; i++) {
		m = f4mon_mul_thread(a->fac, a->mon[i], tmp);
		if (COLUM(m)==1) continue;
		if (cas(&COLUM(m),0,1)) continue;

		/* add m to f4_mused atomically */
		z = f4_uload;
	monom:	y = cas(&f4_uload,z,z+1);
		if (y != z) {
			z = y;
			goto monom;
		}
		f4_mused[z] = m;

		/* we acquired m to process */
		for (k=0; k < f4_bload; k++) {
			b = f4_basis[k];
			if (f4mon_div(m,b->mon[0])) break;
		}
		if (k==f4_bload) continue;

		b = f4_basis[k];
		q = f4mon_quo_thread(m,b->mon[0],tmp);

		/* construct matrix row */
		c = malloc(sizeof(f4row));
		c->len = b->len;
		c->fac = q;
		c->cof = b->cof;
		c->mon = b->mon;
		c->ind = 0;
		c->sug = 0;

		/* add c to f4_array atomically */
		z = f4_aload;
	array:	y = cas(&f4_aload,z,z+1);
		if (y != z) {
			z = y;
			goto array;
		}
		f4_array[z] = c;
	}
	goto start;
done:	free(tmp);
}

/* symbolic preprocessing */
void f4_symbol()
{
	INT i, j, k, l, m, n, x, y, z, s, bs, dp, ncpu, big;
	INT param[2] = {0};
	INT tid[256] = {0};

	/* Free any existing f4_mused array */
	if (f4_mused) {
		free(f4_mused);
		f4_mused = NULL;
	}

	f4_uload = 0;
	f4_usize = 2*f4_aload;
	f4_mused = malloc(f4_usize*sizeof(INT));

	for (i=0; i < f4_aload; i=j) {
		/* count remaining monomials */
		for (l=0, j=i; j < f4_aload; j++) {
			l += f4_array[j]->len;
			if (l >= f4_tsize/2) break;
		}

		/* enlarge hash table */
		while (f4_mload + l >= 3*f4_tsize/4) f4mon_resize();

		/* enlarge array of monoms */
		if (f4_uload + l >= f4_usize) {
			f4_usize = 2*f4_usize + l;
			f4_mused = realloc(f4_mused, f4_usize*sizeof(INT));
		}

		/* enlarge the matrix rows */
		if (f4_aload + l >= f4_asize) {
			f4_asize = 2*f4_asize + l;
			f4_array = realloc(f4_array, f4_asize*sizeof(f4row *));
		}

//printf("i=%lld j=%lld l=%lld\n", i, j, l);
		param[0] = i;
		param[1] = j;

		ncpu = umin(numcpu(),256); // here
		for (k=0; k < ncpu; k++) tid[k] = thread(f4_symbol_thread, (void *)param);
		for (k=0; k < ncpu; k++) join(tid[k]);

		while (f4_mload >= f4_tsize/2) f4mon_resize();
	}

		for (l=i=0; i < f4_aload; i++) l += f4_array[i]->len;
		printf("%lld x %lld with %lld non-zero, %.1f per row\n", 
			(long long int)f4_aload,
			(long long int)f4_uload,
			(long long int)l,
			(double)l/f4_aload
		);
}

#else

/* symbolic preprocessing */
void f4_symbol()
{
	f4row *a, *b, *c;
	INT i, j, k, l, m, n, x, y, z;

	/* Free any existing f4_mused array */
	if (f4_mused) {
		free(f4_mused);
		f4_mused = NULL;
	}

	f4_uload = 0;
	f4_usize = 2*f4_aload;
	f4_mused = malloc(f4_usize*sizeof(INT));

	for (l=i=0; i < f4_aload; i++) {
		a = f4_array[i];
		l += a->len;

		/* make room for new cols */
		if (f4_uload + a->len > f4_usize) {
			f4_usize = 2*f4_usize + a->len;
			f4_mused = realloc(f4_mused,f4_usize*sizeof(INT));
		}
		/* make room for new rows */
		if (f4_aload + a->len > f4_asize) {
			f4_asize = 2*f4_asize + a->len;
			f4_array = realloc(f4_array,f4_asize*sizeof(f4row *));
		}
		for (j=0; j < a->len; j++) {
			m = f4mon_mul(a->fac,a->mon[j]);
			if (COLUM(m)==1) continue;
			f4_mused[f4_uload++] = m;
			COLUM(m) = 1;

			/* divide */
			for (k=0; k < f4_bload; k++) {
				b = f4_basis[k];
				if (f4mon_div(m,b->mon[0])) break;
			}
			if (k == f4_bload) continue;
			b = f4_basis[k];

			/* add row to matrix */
			c = malloc(sizeof(f4row));
			c->len = b->len;
			c->fac = f4mon_quo(m,b->mon[0]);
			c->cof = b->cof;
			c->mon = b->mon;
			c->ind = 0;
			c->sug = 0;
			f4_array[f4_aload++] = c;
		}
	}

		printf("%lld x %lld with %lld non-zero, %.1f per row\n", 
			(long long int)f4_aload,
			(long long int)f4_uload,
			(long long int)l,
			(double)l/f4_aload
		);
}

#endif

#if MULTI_THREAD

/* encode matrix rows */
void f4_encode_thread(void *par)
{
	INT *param = par;
	unsigned char *buf;
	INT i, j, k, l, m, n, *tmp;
	INT ss, st, ks, kt;
	f4row *a;

	buf = 0; tmp = 0;
	st = param[0];	/* row to encode */	
row:	ss = cas(param,st,st+1);
	if (ss != st) {
		st = ss;
		goto row;
	}
	if (st >= f4_aload) goto done;

	n = f4_nvars;
	if (!buf) {
		tmp = malloc(n*sizeof(INT));
		buf = malloc(f4_uload*(sizeof(INT)+1));
	}

	/* encode indices for each matrix row, format: */
	/* [ index0, dif0,...,difk, null, index1, ...] */
	/* indices are words and difs / null are bytes */

	a = f4_array[st]; l = -1;
	for (k=j=0; j < a->len; j++) {
		/* its already in monom hash table */
		/* so this won't enlarge the table */
		m = f4mon_mul_thread(a->fac,a->mon[j],tmp);
		m = COLUM(m);
		if (l == -1) {
			/* store first index */
			*(INT *)(buf+k) = m;
			k += sizeof(INT);
		}
		else if (l > m && l-m <= 255) {
			/* store difference from previous */
			/* in sequence terminated by null */
			buf[k++] = (unsigned char)(l-m);
		}
		else {
			/* null byte to stop the sequence */
			buf[k++] = 0;
			/* store new index */
			*(INT *)(buf+k) = m;
			k += sizeof(INT);
		}
		l = m;
	}
	buf[k++] = 0;	/* extra null */
	a->ind = malloc(k*sizeof(char));
	memcpy(a->ind,buf,k);
	a->siz = k;
	goto row;

done:	free(buf);
	free(tmp);
}

/* form matrix */
void f4_encode()
{
	f4row *a;
	INT i, j, k, m;
	unsigned char *buf;
	INT param[1] = {0};
	INT tid[256] = {0};

	/* sort asc, assign columns */
	f4mon_sort(f4_mused, f4_uload);
	for (i=0; i < f4_uload; i++) {
		m = f4_mused[i];
		COLUM(m) = i;
	}

	/* start threads */
	k = umin(numcpu(),256);
	for (i=0; i < k; i++) {
		tid[i] = thread(f4_encode_thread,(void *)(&param));
	}
	for (i=0; i < k; i++) join(tid[i]);

	for (i=j=k=0; i < f4_aload; i++) {
		j += f4_array[i]->len;
		k += f4_array[i]->siz;
	}

		printf("%.3f bytes per non-zero, matrix encoded in %.3f MB\n",
			(double)k/j,
			(double)k/1024/1024
		);

	/* reset the indices */
	for (i=0; i < f4_uload; i++) {
		m = f4_mused[i];
		COLUM(m) = 0;
	}
}

#else

/* make matrix */
void f4_encode()
{
	f4row *a;
	INT i, j, k, l, m, s, t;
	unsigned char *buf;

	/* sort asc, assign columns */
	f4mon_sort(f4_mused, f4_uload);
	for (i=0; i < f4_uload; i++) {
		m = f4_mused[i];
		COLUM(m) = i;
	}

	/* encode indices for each matrix row, format: */
	/* [ index0, dif0,...,difk, null, index1, ...] */
	/* indices are words and difs / null are bytes */
	buf = malloc(f4_uload*(sizeof(INT)+1));
	for (s=t=0, i=0; i < f4_aload; i++) {
		a = f4_array[i]; l = -1;
		for (k=j=0; j < a->len; j++) {
			m = f4mon_mul(a->fac,a->mon[j]);
			m = COLUM(m);
			if (l == -1) {
				/* store first index */
				*(INT *)(buf+k) = m;
				k += sizeof(INT);
			}
			else if (l > m && l-m <= 255) {
				/* store difference from previous */
				/* in sequence terminated by null */
				buf[k++] = (unsigned char)(l-m);
			}
			else {
				/* null byte to stop the sequence */
				buf[k++] = 0;
				/* store new index */
				*(INT *)(buf+k) = m;
				k += sizeof(INT);
			}
			l = m;
		}
		buf[k++] = 0;	/* extra null */
		a->ind = malloc(k*sizeof(char));
		memcpy(a->ind,buf,k);
		a->siz = k; s += k;
		t += a->len;
	}
	free(buf);

	/* reset the cache */
	for (i=0; i < f4_uload; i++) {
		m = f4_mused[i];
		COLUM(m) = 0;
	}

		printf("%.3f bytes per non-zero, matrix encoded in %.3f MB\n",
			(double)s/t,
			(double)s/1024/1024
		);
}

#endif


/* encode dense vector */
f4row * f4_reduce_export(INT n, INT *vec)
{
	f4row *a = 0;
	INT i, j, k, l;
	unsigned char *buf;
	size_t buf_size;

	if (n < 0) return 0;
	
	/* Fix: Allocate for n+1 elements since we process indices from n down to 0 */
	/* Additional fix: Guard against size_t overflow as identified by O3 */
	size_t elem_count = (size_t)n + 1;
	if (elem_count > SIZE_MAX / (sizeof(INT) + 1)) {
		fprintf(stderr, "f4_reduce_export: allocation size would overflow\n");
		return 0;
	}
	buf_size = elem_count * sizeof(INT) + elem_count;
	buf = malloc(buf_size);
	if (!buf) {
		fprintf(stderr, "f4_reduce_export: malloc failed\n");
		return 0;
	}
	j = k = 0; l = -1; 
	for (i=n; i >= 0; i--) {
		if (vec[i]==0) continue;
		if (l == -1) {
			/* Defensive check: ensure we don't overflow */
			if (k + sizeof(INT) > buf_size) {
				fprintf(stderr, "f4_reduce_export: buffer overflow prevented at first element\n");
				free(buf);
				return 0;
			}
			*(INT *)(buf+k) = i;
			k += sizeof(INT);
		}
		else if (l-i <= 255) {
			/* Defensive check for single byte write */
			if (k >= buf_size) {
				fprintf(stderr, "f4_reduce_export: buffer overflow prevented at gap byte\n");
				free(buf);
				return 0;
			}
			buf[k++] = (unsigned char)(l-i);
		}
		else {
			/* Defensive check for null + INT write */
			if (k + 1 + sizeof(INT) > buf_size) {
				fprintf(stderr, "f4_reduce_export: buffer overflow prevented at large gap\n");
				free(buf);
				return 0;
			}
			buf[k++] = 0;
			*(INT *)(buf+k) = i;
			k += sizeof(INT);
		}
		l = i;
		j++;
	}
	if (!j) goto done;

	/* Defensive check for final null byte */
	if (k >= buf_size) {
		fprintf(stderr, "f4_reduce_export: buffer overflow prevented at final null\n");
		free(buf);
		return 0;
	}
	buf[k++] = 0;	/* extra null */
	a = malloc(sizeof(f4row));
	a->len = j;
	a->sug = 0;
	a->fac = 0;
	a->cof = malloc(j*sizeof(INT));
	/* Guard against divide-by-zero as identified by O3 */
	if (f4_prime == 0) {
		fprintf(stderr, "f4_reduce_export: f4_prime is zero - not initialized?\n");
		free(a->cof);
		free(a);
		free(buf);
		return 0;
	}
	for (j=0, i=n; i >= 0; i--) {
		if (vec[i]==0) continue;
		/* Ensure value is reduced modulo p */
		INT val = vec[i] % f4_prime;
		NORMAL(val, f4_prime);
		a->cof[j++] = val;
		vec[i] = 0;
	}
	a->mon = 0;
	a->ind = malloc(k*sizeof(char));
	memcpy(a->ind,buf,k);

done:	free(buf);
	return a;
}

/* add a*c to dense vector */
INT f4_reduce_import(INT *vec, f4row *a, INT c)
{
	UINT u, v, w;
	INT  x, y, z, i, j, k, l, p, p2;
	p = f4_prime;
	if (!a->len) return -1;
	if (0 && WORDSIZE==64 && p <= 2147483647) {
		i = j = k = 0; p2 = p*p;
		while (j < a->len) {
			z = *(INT *)(a->ind+k); k += sizeof(INT);
		mul:	x = (a->cof[j]) * c;
			y = vec[z] + x - p2;
			NORMAL(y,p2);
			/* Ensure final value is in [0, p) */
			if (y >= p) y %= p;
			vec[z] = y;
			l = z, j++;
			if (a->ind[k]) {
				z = l - a->ind[k]; k++;
				goto mul;
			}
			else k++;
		}
	}
	else {	
		w = ulz(p); u = p << w; v = nreciprocal(u);
		i = j = k = 0;
		while (j < a->len) {
			z = *(INT *)(a->ind+k); k += sizeof(INT);
		mulmod:	x = nmulmod(a->cof[j],c,u,v,w);
			y = vec[z] + x - p;
			NORMAL(y,p);
			vec[z] = y;
			l = z, j++;
			if (a->ind[k]) {
				z = l - a->ind[k]; k++;
				goto mulmod;
			}
			else k++;
		}
	}
	return *(INT *)(a->ind);
}

/* reduce vec using pivots */
INT f4_reduce_vector(INT n, INT *vec, f4row **piv)
{
	f4row *a;
	UINT u, v, w;
	INT  x, y, z, i, j, k, l, c, t, p, p2;
	p = f4_prime;
	if (WORDSIZE==64 && p <= 2147483647) {
		p2 = p*p;
		for (i=n; i >= 0; i--) {
			c = vec[i]; if (!c) continue;
			c %= p;
			NORMAL(c, p);
			vec[i] = c;
			a = piv[i]; if (!a) continue;
			j = k = 0;
			while (j < a->len) {
				z = *(INT *)(a->ind+k); k += sizeof(INT);
			mul:	x = (a->cof[j++]) * c;
				y = vec[z] - x;
				t = a->ind[k++];
				NORMAL(y,p2);
				/* Ensure final value is in [0, p) */
				if (y >= p) y %= p;
				vec[z] = y;
				z = z-t;
				if (t) goto mul;
			}
		}
	}
	else {
		w = ulz(p); u = p << w; v = nreciprocal(u);
		for (i=n; i >= 0; i--) {
			c = vec[i]; if (!c) continue;
			a = piv[i]; if (!a) continue;
			j = k = 0;
			while (j < a->len) {
				z = *(INT *)(a->ind+k); k += sizeof(INT);
			mul2:	x = nmulmod(a->cof[j++],c,u,v,w);
				y = vec[z] - x;
				t = a->ind[k++];
				NORMAL(y,p);
				vec[z] = y;
				z = z-t;
				if (t) goto mul2;
			}
		}
	}
	return n;
}

/* make row into a pivot */
f4row * f4_reduce_monic(f4row *a)
{
	UINT u, v, w;
	INT  p, c, i;
	p = f4_prime;
	c = uinvmod(a->cof[0],p);
	w = ulz(p); u = p << w; v = nreciprocal(u);
	for (i=0; i < a->len; i++) {
		a->cof[i] = nmulmod(a->cof[i],c,u,v,w);
	}
	return a;
}

/* normalform */
void f4_normal(INT l)
{
	f4row **piv, *a, *b;
	INT *vec, i, j, k, m, n, o, p, s;

	p = f4_prime;
	m = f4_aload; n = f4_uload;
	vec = malloc(n*sizeof(INT));
	piv = malloc(n*sizeof(f4row *));
	for (i=0; i < n; i++) vec[i] = 0;
	for (i=0; i < n; i++) piv[i] = 0;

	/* insert pivots */
	for (i=l; i < m; i++) {
		a = f4_array[i];
		k = *(INT *)(a->ind);
		if (piv[k]) continue;
		piv[k] = a;
	}
	/* row reduce by piv */
	for (i=0; i < l; i++) {
		a = f4_array[i];
		s = f4_reduce_import(vec,a,1);
		f4_reduce_vector(s,vec,piv);
		b = f4_reduce_export(s,vec);
		if (!b) b = f4row_new();
		f4_array[i] = b;
		f4row_free(a);
	}
	for (i=l; i < m; i++) {
		a = f4_array[i];
		if (!a->fac) free(a->cof);
		if (!a->fac) free(a->mon);
		free(a->ind); free(a);
	}
	f4_aload = l;
	free(vec);
	free(piv);
}

#if MULTI_THREAD

/* multi-threaded code */
void f4_reduce_thread(void *par)
{
	INT *param = par;
	INT ss, st, sz, bs, zb, zr, *vec;
	INT i, j, k, l, m, n, o, p, s, x;
	f4row **piv, *a, *b;

	piv = (f4row **)(param[4]);
	vec = 0;	/* buffer for reduction */
	zb = param[3];	/* zero-reduction bound */
	bs = param[2];	/* block size to reduce */
	sz = param[1];	/* total number of rows */
	st = param[0];	/* first row in a block */

	/* get block to reduce */
block:	ss = cas(param,st,st+bs);
	if (ss != st) {
		st = ss;
		goto block;
	}
	/* blocks all claimed */
	if (st >= sz) goto done;

	p = f4_prime;
	n = f4_uload;
	if (vec == 0) {
		/* alloc buffer if needed */
		vec = malloc(n*sizeof(INT));
		for (i=0; i < n; i++) vec[i] = 0;
	}

	/* reduce combinations of */
	/* rows st..min(st+bs,sz) */
	k = umin(sz-st,bs); zr = 0;
	for (o=j=0; j < k; j++) {
		/* random combination */
		for (l=0; l < k; l++) {
			x = (urandom() % (p-1)) + 1;
			s = f4_reduce_import(vec, f4_array[st+l], x);
			if (s > o) o = s;
		}
		/* row reduce */
	retry:	f4_reduce_vector(o,vec,piv);
		a = f4_reduce_export(o,vec);
		if (!a) zr++;
		if (zr==zb) break;
		if (!a) continue;
		a = f4_reduce_monic(a);
		x = *(INT *)(a->ind);
		/* assign pivot by compare and swap */
		if (cas((INT *)(piv+x),0,(INT)a) == 0) continue;
		o = f4_reduce_import(vec, a, 1);
		goto retry;
	}

	/* free rows in block */
	for (j=0; j < k; j++) {
		a = f4_array[st+j];
		if (!a->fac) free(a->cof);
		if (!a->fac) free(a->mon);
		free(a->ind);
		free(a);
	}

		printf("%.1f/", 100.0*(st+k)/sz);

	goto block;
done:	free(vec);
}

/* row reduce */
void f4_reduce()
{
	double e = 1e-18;
	f4row **piv, *a, *b;
	INT i, j, k, l, m, n, o, p, s;
	INT x, y, z, zb, zr, bs;
	INT tid[256], param[5];

	p = f4_prime;
	m = f4_aload; n = f4_uload;
	piv = malloc(n*sizeof(f4row *));
	for (i=0; i < n; i++) piv[i] = 0;

	/* select sparse pivots */
	for (l=i=0; i < m; i++) {
		a = f4_array[i];
		k = *(INT *)(a->ind);
		b = piv[k];
		if (!b) piv[k] = a;
		else if (a->len < b->len) {
			piv[k] = a;
			f4_array[l++] = b;
		}
		else f4_array[l++] = a;
	}
	k = m = l;

	/* sort non-pivots */
	if (k < 2) goto done;
sort:	k = 5*(k+1)/13;
	for (i=k-1; i < m; i++) {
		a = f4_array[i];
		x = *(INT *)(a->ind);
		for (j=i; j >= k; j-=k) {
			b = f4_array[j-k];
			y = *(INT *)(b->ind);
			if (x >= y) break;
			f4_array[j] = b;
		}
		f4_array[j] = a;
	}
	if (k > 1) goto sort;
done:
	/* blocksize */
	bs = uroot(m,2);
	if (bs > m-1) bs = m-1;
	if (bs > p-1) bs = p-1;

	/* zb = bound zero reductions */
	for (zb=0; e < 1; zb++) e = e*p;
	if (bs <= zb) bs = 1;

		printf("%lld rows to reduce, blocksize %lld, zero reductions %lld\n", 
			(long long int)m,
			(long long int)bs,
			(long long int)zb
		);

	param[0] = 0;		/* start */
	param[1] = m;		/* bound */
	param[2] = bs;		/* block */
	param[3] = zb;		/* zeros */
	param[4] = U(piv);	/* pivot */
	k = umin(numcpu(),256);
	for (i=0; i < k; i++) {
		tid[i] = thread(f4_reduce_thread,(void *)(&param));
	}
	for (i=0; i < k; i++) join(tid[i]);
	printf("\n");

	/* extract pivots */
	for (i=j=0; i < n; i++) {
		a = piv[i];
		if (!a) continue;
		f4_array[j++] = a;
	}
	f4_aload = j;
	free(piv);
}

#else

/* row reduce */
void f4_reduce()
{
	double e = 1e-18;
	f4row **piv, *a, *b;
	INT *vec, i, j, k, l, m, n, o, p, s, x, y, z, zb, zr, bs;

	p = f4_prime;
	m = f4_aload; n = f4_uload;
	vec = malloc(n*sizeof(INT));
	piv = malloc(n*sizeof(f4row *));
	for (i=0; i < n; i++) vec[i] = 0;
	for (i=0; i < n; i++) piv[i] = 0;

	/* select sparse pivots */
	for (l=i=0; i < m; i++) {
		a = f4_array[i];
		k = *(INT *)(a->ind);
		b = piv[k];
		if (!b) piv[k] = a;
		else if (a->len < b->len) {
			piv[k] = a;
			f4_array[l++] = b;
		}
		else f4_array[l++] = a;
	}
	k = m = l;
	if (k < 2) goto done;
sort:	k = 5*(k+1)/13;
	/* sort remaining rows */
	for (i=k-1; i < m; i++) {
		a = f4_array[i];
		x = *(INT *)(a->ind);
		for (j=i; j >= k; j-=k) {
			b = f4_array[j-k];
			y = *(INT *)(b->ind);
			if (x >= y) break;
			f4_array[j] = b;
		}
		f4_array[j] = a;
	}
	if (k > 1) goto sort;
done:
	/* blocksize */
	bs = uroot(m,2);
	if (bs > m-1) bs = m-1;
	if (bs > p-1) bs = p-1;

	/* zb = bound zero reductions */
	for (zb=0; e < 1; zb++) e = e*p;
	if (bs <= zb) bs = 1;

		printf("%lld rows to reduce, blocksize %lld, zero reductions %lld\n", 
			(long long int)m,
			(long long int)bs,
			(long long int)zb
		);

	/* for each block of rows we reduce */
	/* random linear combinations until */
	/* we obtain a zero vector zb times */
	for (i=0; i < m; i+=bs) {
		k = umin(m-i,bs); zr = 0;
		for (o=j=0; j < k; j++) {
			/* random combination */
			for (l=0; l < k; l++) {
				x = (urandom() % (p-1)) + 1;
				s = f4_reduce_import(vec,f4_array[i+l],x);
				if (s > o) o = s;
			}
			/* row reduce */
			f4_reduce_vector(o,vec,piv);
			a = f4_reduce_export(o,vec);
			if (!a) zr++;
			if (zr==zb) break;
			if (!a) continue;
			a = f4_reduce_monic(a);
			x = *(INT *)(a->ind);
			piv[x] = a;
		}
		/* free rows in block */
		for (j=0; j < k; j++) {
			a = f4_array[i+j];
			if (!a->fac) free(a->cof);
			if (!a->fac) free(a->mon);
			free(a->ind);
			free(a);
		}

			printf("%.1f/",100.0*(i+k)/m);
//			fflush(stdout);
	}
		printf("\n");

	/* extract pivots */
	for (i=j=0; i < n; i++) {
		a = piv[i];
		if (!a) continue;
		f4_array[j++] = a;
	}
	f4_aload = j;
	free(vec);
	free(piv);
}

#endif

/* new pivots */
void f4_filter()
{
	f4row *a, *b;
	INT i, j, m, n, z, w;
	/* take new leading monomials */
	for (i=j=0; i < f4_aload; i++) {
		a = f4_array[i];
		/* new pivots are not of the form */
		/* basis[i]*m, so a->fac is unset */
		if (a->fac == 0) {
			f4_array[j++] = a;
		}
		else {
			/* leave a->mon and a->cof alone */
			/* they're in a basis polynomial */
			if (!a->fac) free(a->cof);
			if (!a->fac) free(a->mon);
			free(a->ind);
			free(a);
		}
	}
	f4_aload = j;
}

#if MULTI_THREAD

/* multi-threaded code */
void f4_jordan_thread(void *par)
{
	INT *param = par;
	INT *vec, ss, st, sz;
	INT i, j, k, l, m, n, o, p, s, x;
	f4row **piv, **red, **arr, *a, *b;

	piv = (f4row **)(param[3]);
	red = (f4row **)(param[2]);
	vec = 0;	/* buffer for reduction */
	sz = param[1];	/* total number of rows */
	st = param[0];	/* first row in a block */

	/* get block to reduce */
block:	ss = cas(param,st,st+1);
	if (ss != st) {
		st = ss;
		goto block;
	}
	/* blocks all claimed */
	if (st >= sz) goto done;

	p = f4_prime;
	n = f4_uload;
	if (vec == 0) {
		/* alloc buffer if needed */
		vec = malloc(n*sizeof(INT));
		for (i=0; i < n; i++) vec[i] = 0;
	}

	/* reduce lower order terms using pivots */
	s = f4_reduce_import(vec, f4_array[st], 1);
	f4_reduce_vector(s-1,vec,piv);
	a = f4_reduce_export(s,vec);
	if (!a) {
		/* Handle case where reduction results in zero polynomial */
		goto block;
	}
	x = *(INT *)(a->ind);
	red[x] = a;

	goto block;
done:	free(vec);
}

/* back substitution */
void f4_jordan()
{
	f4row **piv, **red, *a;
	INT *vec, i, j, k, l, m, n, p;
	INT tid[256], param[4];

	p = f4_prime;
	m = f4_aload; n = f4_uload;
	piv = malloc(n*sizeof(f4row *));
	red = malloc(n*sizeof(f4row *));
	for (i=0; i < n; i++) piv[i] = 0;
	for (i=0; i < n; i++) red[i] = 0;

	/* create pivots */
	for (l=i=0; i < m; i++) {
		a = f4_array[i];
		k = *(INT *)(a->ind);
		piv[k] = a;
	}

	param[0] = 0;		/* start */
	param[1] = m;		/* bound */
	param[2] = U(red);	/* result */
	param[3] = U(piv);	/* pivots */
	k = umin(numcpu(),256);
	for (i=0; i < k; i++) {
		tid[i] = thread(f4_jordan_thread,(void *)(&param));
	}
	for (i=0; i < k; i++) join(tid[i]);

	/* extract pivots */
	for (i=j=0; i < n; i++) {
		a = piv[i];
		if (!a) continue;
		if (!a->fac) free(a->cof);
		if (!a->fac) free(a->mon);
		free(a->ind); free(a);
		f4_array[j++] = red[i];
	}
	f4_aload = j;
	free(piv);
	free(red);
}

#else

/* back substitution */
void f4_jordan()
{
	f4row **piv, *a;
	INT *vec, i, j, k, l, m, n, p;
	p = f4_prime;
	m = f4_aload; n = f4_uload;
	vec = malloc(n*sizeof(INT));
	piv = malloc(n*sizeof(f4row *));
	for (i=0; i < n; i++) vec[i] = 0;
	for (i=0; i < n; i++) piv[i] = 0;

	/* create pivots */
	for (l=i=0; i < m; i++) {
		a = f4_array[i];
		k = *(INT *)(a->ind);
		piv[k] = a;
	}
	/* reduce pivots */
	for (i=0; i < n; i++) {
		a = piv[i];
		if (!a) continue;
		f4_reduce_import(vec,a,1);
		if (!a->fac) free(a->cof);
		if (!a->fac) free(a->mon);
		free(a->ind); free(a);
		piv[i] = 0;
		f4_reduce_vector(i,vec,piv);
		a = f4_reduce_export(i,vec);
		if (a) {
			piv[i] = f4_reduce_monic(a);
		} else {
			piv[i] = 0;
		}
	}
	for (j=i=0; i < n; i++) {
		a = piv[i];
		if (!a) continue;
		f4_array[j++] = a;
	}
	f4_aload = j;
	free(vec);
	free(piv);
}

#endif

/* matrix to polynomials */
void f4_decode()
{
	f4row *a;
	INT i, j, k, l, m, n, x, y, z;
	for (i=0; i < f4_aload; i++) {
		a = f4_array[i];
		n = a->len;
		if (n==0) continue;
		a->mon = malloc(n*sizeof(INT));
		j = k = 0;
		while (j < a->len) {
			z = *(INT *)(a->ind+k); k += sizeof(INT);
		monom:	a->mon[j] = f4_mused[z];
			l = z; j++;
			if (a->ind[k]) {
				z = l - a->ind[k]; k++;
				goto monom;
			}
			else k++;
		}
		free(a->ind);
		a->ind = 0;
		a->fac = 0;
		f4_array[i] = a;
	}
}

/* output basis */
void f4_output()
{
	f4row *a, *b;
	INT i, j, k, l, m, x, y;
	m = f4mon_one();
	if (f4_asize < f4_bload) {
		f4_asize = f4_bload;
		f4_array = realloc(f4_array, f4_asize*sizeof(f4row *));
	}
	f4_aload = f4_bload;
	for (i=0; i < f4_bload; i++) {
		b = f4_basis[i];
		a = malloc(sizeof(f4row));
		a->len = b->len;
		a->sug = b->sug;
		a->fac = m;
		a->cof = b->cof;
		a->mon = b->mon;
		a->ind = 0;
		f4_array[i] = a;
	}

	k = l = f4_aload;
	if (k < 2) goto done;
sort:	k = 5*(k+1)/13;
	/* sort remaining rows */
	for (i=k-1; i < l; i++) {
		a = f4_array[i];
		x = a->mon[0];
		for (j=i; j >= k; j-=k) {
			b = f4_array[j-k];
			y = b->mon[0];
			if (f4mon_cmp(x,y) >= 0) break;
			f4_array[j] = b;
		}
		f4_array[j] = a;
	}
	if (k > 1) goto sort;
done:	return ;
}

/* cull basis */
void f4_remove()
{
	f4row *a, *b;
	INT i, j;
	for (i=0; i < f4_aload; i++) {
		b = f4_array[i];
		if (!b) continue;
		for (j=0; j < i; j++) {
			a = f4_array[j];
			if (!a) continue;
			if (f4mon_div(b->mon[0],a->mon[0])) break;
		}
		if (j < i) {
			f4_array[i] = 0;
			if (!b->fac) free(b->cof);
			if (!b->fac) free(b->mon);
			free(b->ind);
			free(b);
		}
	}
	for (i=j=0; i < f4_aload; i++) {
		if (!f4_array[i]) continue;
		f4_array[j++] = f4_array[i];
	}
	f4_aload = j;
}

/* sort basis */
void f4_sort_basis()
{
	f4row *b, *c;
	INT i, j, k;
	/* sort basis */
	k = f4_bload;
	if (k < 2) return;
sort:	k = 5*(k+1)/13;
	for (i=k-1; i < f4_bload; i++) {
		b = f4_basis[i];
		for (j=i; j >= k; j-=k) {
			c = f4_basis[j-k];
			if (b->len > c->len) break;
		//	if (b->sug > c->sug) break;
			f4_basis[j] = f4_basis[j-k];
		}
		f4_basis[j] = b;
	}
	if (k > 1) goto sort;
}

void f4_sort_pairs()
{
	f4syz *b, *c;
	INT i, j, k;
	/* sort pairs */
	k = f4_pload;
	if (k < 2) return;
sort:	k = 5*(k+1)/13;
	for (i=k-1; i < f4_pload; i++) {
		b = f4_pairs[i];
		for (j=i; j >= k; j-=k) {
			c = f4_pairs[j-k];
			if (f4mon_cmp(b->lcm,c->lcm) > 0) break;
			f4_pairs[j] = f4_pairs[j-k];
		}
		f4_pairs[j] = b;
	}
	if (k > 1) goto sort;
}

/* Declare function to reset f4_select statics */
void f4_select_init_statics(void);

void f4mod_init(INT n, INT e, INT p)
{
	f4_prime = p;
	f4mon_init(n,e);
	f4_aload = f4_bload = f4_pload = f4_eload = 0;
	f4_asize = f4_bsize = f4_psize = f4_esize = 30;
	f4_uload = 0;
	f4_usize = 0;
	f4_array = malloc(f4_asize*sizeof(f4row *));
	f4_basis = malloc(f4_bsize*sizeof(f4row *));
	f4_pairs = malloc(f4_psize*sizeof(f4syz *));
	f4_extra = malloc(f4_esize*sizeof(f4row *));
	f4_mused = NULL;  /* Initialize to NULL, will be allocated in f4_symbol */
	f4_select_init_statics();  /* Reset static variables for new computation */
}

/* Forward declaration to reset static variables in f4_select */
void f4_select_reset_statics(void);

void f4mod_free()
{
	INT i;
	for (i=0; i < f4_aload; i++) f4row_free(f4_array[i]);
	for (i=0; i < f4_bload; i++) f4row_free(f4_basis[i]);
	for (i=0; i < f4_eload; i++) f4row_free(f4_extra[i]);
	for (i=0; i < f4_pload; i++) free(f4_pairs[i]);
	free(f4_array);	free(f4_basis); free(f4_pairs); free(f4_extra);
	free(f4_mused);  /* Free the monomials array */
	f4_array = NULL;
	f4_basis = NULL;
	f4_pairs = NULL;
	f4_extra = NULL;
	f4_mused = NULL;
	f4_aload = f4_bload = f4_pload = f4_eload = 0;
	f4_asize = f4_bsize = f4_psize = f4_esize = 0;
	f4_uload = f4_usize = 0;  /* Reset monomial array tracking */
	f4_select_reset_statics();  /* Reset static variables in f4_select */
	f4mon_free();
	f4_prime = 0;
}

/* F4 mod p */
void f4gb_mod()
{
	f4row *b, *c;
	INT d, e, i, j, k, n, p, step;
	double tt, t0, t1;

	n = f4_nvars;
	e = f4_nelim;
	p = f4_prime;
	tt = realtime();

    if (getenv("RUR_VERBOSE_F4")) printf("F4 eliminate %lld/%lld variables mod p=%lld, %lld threads\n", (long long int)e, (long long int)n, (long long int)p, (long long int)numcpu());

	for (i=0; i < f4_aload; i++) {
		b = f4_array[i];
		b = f4_reduce_monic(b);
		f4_update(b);
	}
	step = 0;
	while (f4_pload > 0) {
		t0 = realtime();

		step++;
		printf("\nSTEP %lld\n", (long long int)step);

		f4_sort_pairs();
		d = f4_sugar();

		f4_memory();	/* free up memory for alg */
		f4_select();	/* select pairs to reduce */
		f4_symbol();	/* symbolic preprocessing */
		f4_encode();	/* encode polys in matrix */
		f4_reduce();	/* fast Gauss elimination */
		f4_filter();	/* select new polynomials */
		f4_jordan();	/* fast back substitution */
		f4_decode();	/* convert matrix to poly */

		for (i=0; i < f4_aload; i++) {
			b = f4_array[i];
			b->sug = d;
			f4_update(b);
		}

		t1 = realtime();
		printf("new=%lld basis=%lld extra=%lld, step time=%.3f sec\n", (long long int)f4_aload, (long long int)f4_bload, (long long int)f4_eload, t1-t0);
	}
	printf("\nDONE inter-reduce\n");

	f4_output();
	f4_remove();
	f4_symbol();
	f4_encode();
	f4_jordan();
	f4_decode();
	f4_remove();

	tt = realtime()-tt;

	printf("%lld basis elements, %.3f sec\n", (long long int)f4_aload, tt);
}



/*
	import polynomials
*/

#define MAXVARS 1024
char * vars[MAXVARS] = {0};

/* get integer */
INT getint(char *s, int *l)
{
	int i, j;
	INT c = 0;
	for (i=0; s[i]; i++) {
		if ('0' <= s[i] && s[i] <= '9') {
			c = 10*c + (s[i] - '0');
		}
		else break;
	}
	*l = i; return c;
}

/* get variable */
char * getvar(char *s, int *l)
{
	char *v;
	int i, j;
	for (i=0; s[i]; i++) {
		/* variables begin with a letter */
		/* or underscore but can include */
		/* digits after the first letter */
		if (('_' <= s[i] && s[i] <= 'z')
		 || ('A' <= s[i] && s[i] <= 'Z')
		 || ('0' <= s[i] && s[i] <= '9' && i > 0)
		) continue; else break;
	}
	if (i==0) return 0;
	v = malloc((i+1)*sizeof(char));
	for (j=0; j < i; j++) v[j] = s[j]; v[j] = 0;
	*l = i; return v;
}

/* import one term of polynomial */
INT getmon(char *s, INT *c, int *l)
{
	INT z[MAXVARS] = {0};
	INT b, e, m, n, p, t;
	int i, j, k;
	char *v;

	n = f4_nvars;
	p = f4_prime; 
	*c = 1; i = 0; t = +1;
next:	switch (s[i]) {
	case '+':	if (i > 0) goto done;
			i++; t=+1; break;
	case '-':	if (i > 0) goto done;
			i++; t=-1; break;
	case '*':	i++; break;
	case '/':	goto fail;
	case '\n':	goto done;
	}

	/* coefficient? */
	b = getint(s+i,&j);
	if (j > 0 && t==+1) {
		*c = b % p;
		NORMAL(*c, p);
	}
	if (j > 0 && t==-1) *c = (p-b) % p;
	i += j; if (j) goto next;

	/* var^exponent */
	v = getvar(s+i,&j);
	if (!v) goto fail;
	i += j; e = 1;
	if (s[i] == '^') {
		i++;
		e = getint(s+i,&j);
		if (j==0) goto fail;
		i += j;
	}

	/* put exponent in z */
	for (k=0; k < n; k++) {
		if (strcmp(v,vars[k])) continue;
		z[k] = e; break;
	}
	free(v);
	if (k==n) goto fail;
	goto next;

done:	*l = i;
	return f4mon_new(z);

fail:	printf("error: can't parse term\n");
	return 0;
}


/* import expanded polynomial */
f4row * getpol(char *s, int *l)
{
	f4row *b;
	INT c, m;
	int i, j, k;

	/* count the number of terms */
	for (k=1, i=0; s[i] != '\n'; i++) {
		if (s[i] == '+' || s[i] == '-') k++;
	}

	b = f4row_new();
	b->len = b->fac = 0;
	b->cof = malloc(k*sizeof(INT));
	b->mon = malloc(k*sizeof(INT));
	b->ind = 0;

	for (i=k=0; s[i] != '\n'; i+=j) {
		if (s[i]==0) return 0;
		m = getmon(s+i,&c,&j);
		if (m==0) break;
		if (c==0) continue;
		b->cof[k] = c;
		b->mon[k] = m;
		b->len =  ++k;
	}
	f4row_sort(b);

	*l = i;
	return b;
}


/* write polynomial to file */
void putpol(f4row *b, FILE *out)
{
	INT *e, i, j, m, n, c;
	n = f4_nvars;
	for (i=0; i < b->len; i++) {
		c = b->cof[i];
		m = b->mon[i];
		e = EXPON(m);
		fprintf(out,"%+lld",(long long int)c);
		for (j=0; j < n; j++) {
			if (e[j] == 0) continue;
			if (e[j] == 1) fprintf(out,"*%s",vars[j]);
			if (e[j] >= 2) fprintf(out,"*%s^%lld",vars[j],(long long int)e[j]);
		}
	}
}

#undef NORMAL
#undef EXPON
#undef COLUM
#undef HASHVAL


