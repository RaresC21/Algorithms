/*
    Problem from hackerrank.
    https://www.hackerrank.com/contests/world-codesprint-8/challenges/prime-digit-sums/copy-from/8166958
*/
#include <bits/stdc++.h>
using namespace std;
#define int long long int

#define MOD 1000000007
#define MAX 400001
int N;
bool prime[200];
vector<int> attach[1000001];
int ans[MAX];
int memo[2][MAX];
vector<int> good, nonzero;

#undef int
int main() {
#define int long long int

    for (int i = 3; i < 100; i++) {
        bool p = true;
        for (int k = 2; k < i; k++) {
            if (i % k == 0) p = false;
        }
        if (p) prime[i] = true;
    }
    prime[2] = true;

    int amnt = 0;
    for (int a = 0; a < 10; a++) {
        if (a != 0) ans[1]++;
        for (int b = 0; b < 10; b++) {
            if (a != 0) ans[2] ++;
            for (int c = 0; c < 10; c++) {
                if (!prime[a+b+c]) continue;
                if (a != 0)
                  ans[3]++;
                for (int d = 0; d < 10; d++) {
                    if (!prime[a+b+c+d]) continue;
                    if (!prime[b+c+d]) continue;
                    if (a != 0)
                        ans[4]++;
                    for (int e = 0; e < 10; e++) {
                        if (!prime[a+b+c+d+e]) continue;
                        if (!prime[b+c+d+e]) continue;
                        if (!prime[c+d+e]) continue;
                        int num = e + 10 * (d + 10 * (c + 10 * (b + 10 * a)));
                        good.push_back(num);
                        if (a != 0) {
                            nonzero.push_back(num);
                            ans[5]++;
                        }
                    }
                }
            }
        }
    }

    for (int a : good) {
        for (int b : good){
            if (a == b) continue;
            // last four digits of a must be same as first four digits of b
            int last_ = a % 10000;
            int first_ = b / 10;
            if (last_ == first_) {
                attach[a].push_back(b);
            }
        }
    }

    for (int n : nonzero) {
        memo[1][n] = 1;
    }
    for (int i = 6; i < MAX; i++) {
        for (int k : good) {
            for (int j : attach[k]) {
                int pp = i % 2;
                memo[pp][j] = (memo[pp][j] + memo[(pp+1)%2][k]) % MOD;
            }
            memo[(i+1)%2][k] = 0;
        }
        for (int g : good) {
            ans[i] = (ans[i] + memo[i%2][g]) % MOD;
        }
        ans[i] = (ans[i] + MOD) % MOD;
    }

    int T; cin >> T;
    while (T) {
        T--;
        cin >> N;
        cout << ans[N] << "\n";
    }

    return 0;
}