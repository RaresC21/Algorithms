#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5

// a_0 + a_1*x + a_2 * x^2 + ...
lli hash_string1(const string& s, lli MOD) {
    lli h = 0;
    lli a = 1;
    for (int i =  0; i < s.size(); i++) {
        h = (h + (a * s[i]) % MOD) % MOD;
        a = (a * PRIME) % MOD;
    }
    return h;
}

// a_0 * x^len  +a_1 * x^(len-1) + ... a_(len-1) * x + a_len
lli hash_string2(const string& s, lli MOD) {
    lli h = 0;
    lli a = 1;
    for (int i =  0; i < s.size(); i++) {
        h = ((h * PRIME + s[i]) % MOD) % MOD;
    }
    return h;
}