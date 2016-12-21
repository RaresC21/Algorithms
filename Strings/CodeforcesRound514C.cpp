#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5
#define MAX 200001

map<lli, bool> hash1, hash2;

lli hash_string(const string& s, lli MOD) {
    lli h = 0;
    lli a = 1;
    for (int i =  0; i < s.size(); i++) {
        h = ((h * PRIME + s[i]) % MOD) % MOD;
    }
    return h;
}

int main() {

    int N, M; cin >> N >> M;

    string s;
    for (int i = 0; i < N; i++) {
        cin >> s;
        hash1[hash_string(s, MOD1)] = hash2[hash_string(s, MOD2)] = true;
    }

    vector<lli> x1(1000000), x2(1000000);
    for (int i = 0; i < M; i++) {
        cin >> s;
        lli h1 = hash_string(s, MOD1);
        lli h2 = hash_string(s, MOD2);

        x1[s.size() - 1] = x2[s.size() - 1] = 1;
        for (int i = s.size() - 2; i >= 0; i--) {
            x1[i] = (x1[i+1] * PRIME) % MOD1;
            x2[i] = (x2[i+1] * PRIME) % MOD2;
        }

        bool found = false;
        for (int k = 0; k < s.size() && !found; k++) {
            for (int c = 'a'; c <= 'c'; c++) {
                if (c == s[k]) continue;

                lli new_1 = (h1 - (x1[k] * (s[k] - c)) % MOD1) % MOD1;
                lli new_2 = (h2 - (x2[k] * (s[k] - c)) % MOD2) % MOD2;

                new_1 = (new_1 + MOD1) % MOD1;
                new_2 = (new_2 + MOD2) % MOD2;

                if (hash1.count(new_1) && hash2.count(new_2)) {
                    found = true;
                    break;
                }
            }
        }

        if (found) {
            cout << "YES\n";
        } else {
            cout << "NO\n";
        }
    }

    return 0;
}