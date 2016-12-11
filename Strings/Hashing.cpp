#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5

lli hash_string (const string& s, lli MOD) {
    lli h = 0;
    lli a = 1;
    for (int i =  0; i < s.size(); i++) {
        h = (h + (a * s[i]) % MOD) % MOD;
        a = (a * PRIME) % MOD;
    }
    return h;
}