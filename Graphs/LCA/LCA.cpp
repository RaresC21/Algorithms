#include <bits/stdc++.h>
using namespace std;

#define MAX 11001
int N, M;
int parent[MAX];
int A[MAX][30]; // A[k][i] 2^i-th ancestor of k
int depth[MAX], dist[MAX];
vector<int> edge[MAX];
map<pair<int, int>, int> weight;

void dfs(int x, int depth_, int dist_) {
    if (depth[x] != -1) return;
    depth[x] = depth_;
    dist[x] = dist_;
    for (int next_ : edge[x]) {
        if (parent[next_] == -1)
            parent[next_] = x;
        dfs(next_, depth_+1, dist_ + weight[make_pair(x, next_)]);
    }
}

void preprocessLCA(int root) {
    for (int i = 0; i <= N; i++) {
        for (int k = 0; k < 30; k++) A[i][k] = -1;
        depth[i] = parent[i] = -1;
    }
    dfs(root, 0, 0);
    depth[0] = 0;
    parent[root] = root;

    for (int i = 1; i <= N; i++)
        A[i][0] = parent[i];

    for (int k = 0; (1 << k) <= N; k++)
        for (int i = 1; i <= N; i++)
            if (A[i][k] == -1)
                A[i][k] = A[A[i][k-1]][k-1];
}

int query(int p, int q) {

    if (depth[p] < depth[q])
        swap(p, q);

    int log;
    for (log = 1; (1 << log) <= N; log++);
    log--;

    for (int i = log; i >= 0; i--)
        if (depth[p] - (1 << i) >= depth[q])
            p = A[p][i];

    if (p == q) return p;

    for (int i = log; i >= 0; i--)
        if (A[p][i] != -1 && A[p][i] != A[q][i])
            p = A[p][i], q = A[q][i];

    return parent[p];
}

int kth_ancestor(int p, int K) {
    if (K == 1) return parent[p];
    else if (K == 0) return p;

    while (K > 0) {
        int i = 0;
        for (; (1 << i) <= K; i++);
        K -= (1 << (i-1));
        return kth_ancestor(A[p][i-1], K);
    }
}

int get_dist(int a, int b) {
    int lca = query(a, b);
    int d = dist[a] + dist[b] - 2 * dist[lca];
    return d;
}

int main() {

    cin >> N;
    for (int i = 1; i < N; i++) {
        int a, b, w; cin >> a >> b >> w; // edges have weight.
        weight[make_pair(a,b)] = weight[make_pair(b,a)] = w;
        edge[a].push_back(b);
        edge[b].push_back(a);
    }

    preprocessLCA(1);
    // do queries here...

    return 0;
}