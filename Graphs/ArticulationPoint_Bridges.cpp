#include <bits/stdc++.h>
using namespace std;

#define MAX 1000
int N, M;
int edge_count = 0, root_children = 0;
int dfs_num[MAX], dfs_low[MAX], dfs_parent[MAX];
vector<int> path[MAX];

struct Edge {
    int a, b;
    Edge(int a, int b) : a(a), b(b) {}
};
vector<Edge> bridge;
vector<int> articulation;

void dfs(int cur) {
    dfs_low[cur] = dfs_num[cur] = edge_count ++;
    for (int v : path[cur]) {
        if (dfs_num[v] == -1) { // unvisited
            dfs_parent[v] = cur;
            if (cur == root) root_children++;

            dfs(v);

            if (dfs_low[v] >= dfs_num[cur]) {
                articulation.push_back(cur);
            }
            if (dfs_low[v] > dfs_num[cur]) {
                bridge.push_back(Edge(v, cur));
            }
            dfs_low[cur] = min(dfs_low[cur], dfs_low[v]);

        } else if (v != dfs_parent[cur]) {
            dfs_low[cur] = min(dfs_low[cur], dfs_num[v]);
        }
    }
}

int main() {

    cin >> N >> M;
    for (int i = 0; i < M; i++) {
        int a, b; cin >> a >> b;
        path[a].push_back(b);
        path[b].push_back(a);
    }

    return 0;
}