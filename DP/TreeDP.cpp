/*
A tree, t, has n vertices numbered from 1 to n and is rooted at vertex 1. Each
vertex i has an integer weight, w_i, associated with it, and t's total weight
is the sum of the weights of its nodes. A single remove operation removes the
subtree rooted at some arbitrary vertex u from tree .
Given t, perform up to k remove operations so that the total weight of the
remaining vertices in t is maximal. Then print t's maximal total weight on a new
line.
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;
#define int long long int

#define MAX 101101
int nxt[MAX], val[MAX];
vector<int> arr;
vector<int> path[MAX];
bool met[MAX];
int dp[MAX][211];

int dfs (int cur) {
    if (met[cur]) return 0;
    met[cur] = true;

    arr.push_back(cur);
    int num = 1;
    for (int adj : path[cur]) {
        num += dfs(adj);
    }

    return nxt[cur] = num;
}

int main() {

    arr.push_back(-1);
    int N, K; cin >> N >> K;
    for (int i = 1; i <= N; i++)
        cin >> val[i];
    for (int i = 1; i < N; i++) {
        int a, b; cin >> a >> b;
        path[a].push_back(b);
        path[b].push_back(a);
    }

    dfs(1);

    for (int i = 0; i <= N; i++)
        for (int k = 0; k <= 200; k++)
            dp[i][k] = -1000000000000000000;

    dp[N][0] = val[arr[N]];
    for (int i = 1; i <= K; i++)
        dp[N][i] = 0;

    int ans = 0;
    for (int i = N - 1; i >= 1; i--) {
        for (int k = 0; k <= K; k++) {
            dp[i][k] = max(dp[i][k], val[arr[i]] + dp[i + 1][k]);
            if (k > 0)
                dp[i][k] = max(dp[i][k], dp[nxt[arr[i]] + i][k - 1]);
            if (i == 1) {
                ans = max (ans, dp[i][k]);
            }
        }
    }

    cout << max(ans, dp[1][0]) << "\n";

    return 0;
}