/*
We can compute nCr using the usual formula n! / (r! * (n - r)!).
To do this, we must store our answer mod P to avoid overflow. 

We can also compute it with Pascal's triangle.
*/

#define MOD 1000000007
#define MAX 2001
lli fact[MAX];

lli exponent_(lli a, lli x) { // returns a^x
	lli ans = 1;
	while (x > 0) {
		if (x % 2 == 1)
			ans = (ans * a) % MOD;
		x /= 2;
		a = (a * a) % MOD;
	}
	return ans;
}

lli modular_inverse_(lli a, lli b) {
	return ((a % MOD) * (exponent_(b, MOD - 2) % MOD)) % MOD;
}

lli comb_(lli n, lli r) {
	calculate_factorial();
	return (modular_inverse_(fact[n], (fact[r] * fact[n - r]) % MOD)) % MOD;
}

int C[MAX];
lli psacal_comb_(lli n, int k)
{
	for (int j = 0; j < MAX; j++)
		C[j] = 0;

	C[0] = 1;  // nC0 is 1

	for (int i = 1; i <= n; i++)
	{
		// Compute next row of pascal triangle using
		// the previous row
		for (lli j = min(i, k); j > 0; j--)
			C[j] = C[j] + C[j - 1];
	}
	return C[k];
}

void calculate_factorial() {
	fact[0] = 1;
	for (int i = 1; i < MAX; i++)
		fact[i] = fact[i - 1] * i;
}
