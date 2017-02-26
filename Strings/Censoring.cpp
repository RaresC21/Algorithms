/*
    USACO Contest - 2015 February Contest.
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5
#define MAX 100001

vector<pair<pair<lli, lli>, int> > pref_hash;

// returns a^x mod p
lli exponent_(lli a, lli x, lli p) {
    lli ans = 1;
    while (x > 0) {
        if (x % 2 == 1)
            ans = (ans * a) % p;
        x /= 2;
        a = (a * a) % p;
    }
    return ans;
}

lli modular_inverse_(lli a, lli b, lli p) {
    return ((a % p) * (exponent_(b, p - 2, p) % p)) % p;
}

lli hash_string (const string& s, lli MOD) {
    lli h = 0;
    lli a = 1;
    for (int i =  0; i < s.size(); i++) {
        h = (h + (a * s[i]) % MOD) % MOD;
        a = (a * PRIME) % MOD;
    }
    return h;
}

int main() {

    ifstream cin ("censor.in");
    ofstream cout ("censor.out");

    string S; cin >> S;
    string T; cin >> T;

    lli t_hash = hash_string(T, MOD1);
    lli len = T.length();

    lli XX = 1;
    pref_hash.push_back(make_pair(make_pair(0, 1LL), -1));
    for (int i = 0; i < S.size(); i++) {
        lli XX = pref_hash.back().first.second;
        pref_hash.push_back(make_pair(make_pair(pref_hash.back().first.first + XX * S[i], (XX * PRIME) % MOD1), i));

        lli sz = pref_hash.size();
        if (sz > len) {
            lli cur = pref_hash.back().first.first - pref_hash[sz - 1 - len].first.first;
            cur = modular_inverse_(cur, exponent_(PRIME, sz - len - 1, MOD1), MOD1);
            cur = (cur + MOD1) % MOD1;

            if (cur == t_hash) {
                pref_hash.resize(pref_hash.size() - len);
            }
        }
    }

    for (auto h : pref_hash) {
        if (h.second == -1) continue;
        cout << S[h.second];
    }
    cout << "\n";

    return 0;
}
