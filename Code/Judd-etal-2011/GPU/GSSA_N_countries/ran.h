#include <math.h>

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline REAL doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

struct Normaldev : Ran {
	REAL mu,sig;
	Normaldev(REAL mmu, REAL ssig, Ullong i)
	: Ran(i), mu(mmu), sig(ssig){}
	REAL dev() {
		REAL u,v,x,y,q;
		do {
			u = doub();
			v = 1.7156*(doub()-0.5);
			x = u - 0.449871;
			y = fabs(v) + 0.386595;
			q = pow(x,2) + y*(0.19600*y-0.25472*x);
		} while (q > 0.27597
			 && (q > 0.27846 || pow(v,2) > -4.*log(u)*pow(u,2)));
		return mu + sig*v/u;
	}
};
