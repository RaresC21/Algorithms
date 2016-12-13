#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5

vector<int> matching_positions(const string& S, const string& T) {

    lli x_power = 1;
    for (int i = 0; i < T.size() - 1; i++)
        x_power = (x_power * PRIME) % MOD1;

    lli s = 0, t = 0; // hash value for S and T respectively
    for (int i = 0; i < T.size(); i++) {
        s = (s * PRIME + S[i]) % MOD1;
        t = (t * PRIME + T[i]) % MOD1;
    }

    vector<int> match;
    int N = S.size(), M = T.size();
    for (int i = M; i < N; i++) {
        if (s == t)
            match.push_back(i - M + 1);

        s = (PRIME * (s - x_power * S[i-M]) + S[i]) % MOD1;
        s = (s + MOD1) % MOD1;
    }
    if (s == t) {
        match.push_back(N - M + 1);
    }


    return match;
}

int main() {

    string S, T; cin >> S >> T; // text S, target T
    for (int p : matching_positions(S, T)) {
        cout << p << "\n";
    }

    return 0;
}