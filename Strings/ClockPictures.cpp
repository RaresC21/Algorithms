/*
    Nordic Collegiate Programming Contest 2014
    - Clock Pictures.
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define mp make_pair
#define pb push_back
#define MAX 200001

int c1[MAX], c2[MAX];

#define MOD1 999999937
#define MOD2 1000000007
#define PRIME 5

bool matching_positions(const vector<int>& S, const vector<int>& T) {

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
        if (s == t) {
            return true;
        }

        s = (PRIME * (s - x_power * S[i - M]) + S[i]) % MOD1;
        s = (s + MOD1) % MOD1;
    }
    if (s == t) {
        match.push_back(N - M + 1);
    }

    return false;
}

int main() {
    ios_base::sync_with_stdio(false);

    int N; cin >> N;
    for (int i = 0; i < N; i++)
        cin >> c1[i];
    for (int i = 0; i < N; i++)
        cin >> c2[i];

    sort(c1, c1 + N);
    sort(c2, c2 + N);

    vector<int> angle1(N), angle2(N * 2);
    for (int i = 1; i < N; i++) {
        angle1[i] = c1[i] - c1[i - 1];
        angle2[i] = c2[i] - c2[i - 1];
    }
    angle1[0] = 360000 - (c1[N - 1] - c1[0]);
    angle2[0] = 360000 - (c2[N - 1] - c2[0]);

    for (int i = N; i < 2 * N; i++) {
        angle2[i] = angle2[i - N];
    }

    if (matching_positions(angle2, angle1)) {
        cout << "possible\n";
    } else {
        cout << "impossible\n";
    }

    return 0;
}