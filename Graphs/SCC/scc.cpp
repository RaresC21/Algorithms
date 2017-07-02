#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define MAX 100
int energy[MAX];
int dfs_counter = 0, cur = 0;
int dfs_num[MAX], dfs_low[MAX];
bool visited[MAX];
vector<int> S;
vector<int> scc[MAX], path[MAX];

void tarjan(int u) {
    dfs_low[u] = dfs_num[u] = dfs_counter++;
    S.push_back(u);
    visited[u] = true;
    for (int v : path[u]) {
        if (dfs_num[v] == -1)
            tarjan(v);
        if (visited[v])
            dfs_low[u] = min(dfs_low[u], dfs_low[v]);
    }
    if (dfs_low[u] == dfs_num[u]) {
        while (true) {
            int v = S.back();
            S.pop_back();
            visited[v] = false;
            scc[cur].push_back(v);
            if (u == v) break;
        }
        cur++;
    }
}

int main() {

    int N; cin >> N;
    for (int i = 1; i <= N; i++) {
        cin >> energy[i];
        int M; cin >> M;
        for (int k = 0; k < M; k++) {
            int a; cin >> a;
            path[i].push_back(a);
        }
    }

    memset(dfs_num, -1, sizeof(dfs_num));
    for (int i = 1; i <= N; i++) {
        if (dfs_num[i] == -1)
            tarjan(i);
    }

    return 0;
}