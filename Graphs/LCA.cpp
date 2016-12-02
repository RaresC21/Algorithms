#include <bits/stdc++.h>
using namespace std;

#define MAX 100001
int N, M;
int parent[MAX];
int A[MAX][30]; // A[k][i] 2^i-th ancestor of k
int depth[MAX];
vector<int> edge[MAX];

void dfs(int x, int d) {
    if (depth[x] != -1) return;
    depth[x] = d;
    for (int next_ : edge[x]) {
        if (parent[next_] == -1)
            parent[next_] = x;
        dfs(next_, d+1);
    }
}

void preprocessLCA() {
    memset(A, -1, sizeof(A));
    memset(depth, -1, sizeof(depth));
    memset(parent, -1, sizeof(parent));
    dfs(0, 0);

    for (int i = 0; i < N; i++)
        A[i][0] = parent[i];

    for (int k = 0; (1 << k) < N; k++)
        for (int i = 0; i < N; i++)
            if (A[i][k] == -1)
                A[i][k] = A[A[i][k-1]][k-1];
}

int query(int p, int q) {

    if (depth[p] < depth[q])
        swap(p, q);

    int log;
    for (log = 1; (1 << log) <= N; log++);
    log--;

    // we find the ancestor of node p that is situated on the same
    // level with q
    for (int i = log; i >= 0; i--)
        if (depth[p] - (1 << i) >= depth[q])
            p = A[p][i];

    if (p == q)
        return p;

    for (int i = log; i >= 0; i--)
        if (A[p][i] != -1 && A[p][i] != A[q][i])
            p = A[p][i], q = A[q][i];

    return parent[p];
}

int main() {

    cin >> N;
    for (int i = 1; i < N; i++) {
        int a, b; cin >> a >> b;
        edge[a].push_back(b);
        edge[b].push_back(a);
    }

    preprocessLCA();
    // do queries here...

    return 0;
}