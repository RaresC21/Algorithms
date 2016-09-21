/*
Determines the Nth Fibonacci number. Uses exponentiation by squaring.
Utilizes the following form of determining Fibonacci numbers:
The matrix
                        [ F(n+1)    F(n)   ]
                        [ F(n)      F(n-1) ]

is equal to  A^N where A is the following matrix:
                             [ 1   1 ]
                             [ 1   0 ]

Runs in O(LogN) time.
*/

#include "Include.h"

class Fibonacci {

#define MOD 1000000007

private:

	void multiply(lli F[2][2], lli M[2][2])
	{
		lli x = F[0][0] * M[0][0] + F[0][1] * M[1][0];
		lli y = F[0][0] * M[0][1] + F[0][1] * M[1][1];
		lli z = F[1][0] * M[0][0] + F[1][1] * M[1][0];
		lli w = F[1][0] * M[0][1] + F[1][1] * M[1][1];

		F[0][0] = x%MOD;
		F[0][1] = y%MOD;
		F[1][0] = z%MOD;
		F[1][1] = w%MOD;
	}

	void power(lli F[2][2], lli n)
	{
		if (n == 0 || n == 1)
			return;
		lli M[2][2] = { 
			{ 1,1 },
			{ 1,0 } 
		};

		power(F, n / 2);
		multiply(F, F);

		if (n % 2 != 0)
			multiply(F, M);
	}

public :

	lli fib(lli n) {
		lli F[2][2] = { 
			{ 1,1 },
			{ 1,0 } 
		};
		if (n == 0)
			return 0;
		power(F, n - 1);
		return F[0][0];
	}
};
