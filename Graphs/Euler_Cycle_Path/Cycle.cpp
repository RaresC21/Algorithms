#include <bits/stdc++.h>
using namespace std;

int circuit[200001];
int N, E, pos = 0;

struct Edge {
    int a, b;
    Edge() {}
    Edge(int a, int b) : a(a), b(b) {}
};

Edge edge[200001];
set<int> path[100001];

void find_circuit(int node) {
    if (path[node].empty()) {
        circuit[pos] = node;
        pos++;
    } else {
        while (!path[node].empty()) {
            int indx = *path[node].begin();
            int adj = edge[indx].a;
            if (adj == node) adj = edge[indx].b;

            path[node].erase(indx);
            path[adj].erase(indx);
            find_circuit(adj);
        }
        circuit[pos] = node;
        pos++;
    }
}

bool met[100001];
void connected(int node) {
    if (met[node]) return;
    met[node] = true;
    for (int indx : path[node]) {
        int adj = edge[indx].a;
        if (node == adj) adj = edge[indx].b;
        connected(adj);
    }
}

bool euler_circuit() {
    connected(1);
    int k = 1;
    for (; k <= N; k++)
        if (!met[k]) return false;

    for (int i = 1; i <= N; i++)
        if (path[i].size() % 2 != 0)
            return false;
    return true;
}

int main() {

    cin >> N >> E;
    for (int i = 0; i < E; i++) {
        int a, b; cin >> a >> b;
        path[a].insert(i);
        path[b].insert(i);
        edge[i] = Edge(a, b);
    }

    return 0;
}