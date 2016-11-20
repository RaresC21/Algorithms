#include <bits/stdc++.h>
using namespace std;

// ----------------- DISJOINT SET UNION -------------------------
#define MAX 3001
int parent[MAX];
int rank_[MAX];

int find_set (int a) {
    if (a != parent[a])
        parent[a] = find_set(parent[a]);
    return parent[a];
}

void make_set(int a) {
    parent[a] = a;
    rank_[a] = 0;
}

void merge_set (int a, int b) {
    int root_a = find_set (a);
    int root_b = find_set (b);

    if (root_a == root_b) return;

    if (rank_[root_a] > rank_[root_b])
        parent[root_b] = root_a;
    else
        parent[root_a] = root_b;
    if (rank_[root_a] == rank_[root_b])
        rank_[root_b]++;
}

// -------------- KRUSKAL'S ALGORITHM, MINIMAL SPANNING TREE ------------------

int N, M; // number of nodes and edges respectively
struct Edge {
    int weight;
    int a, b;
    Edge (int a, int b, int weight) : a(a), b(b), weight(weight) {}
    bool operator < (Edge const &other) const {
        return weight < other.weight;
    }
};

vector<Edge> edges, tree;

int kruskal() {
    sort(edges.begin(), edges.end());
    for (int i = 1; i <= N; i++)
        make_set(i);

    int tree_cost = 0;
    for (int i = 0; i < M; i++) {
        if (find_set(edges[i].a) != find_set(edges[i].b)) {
            tree_cost += edges[i].weight;
            tree.push_back(edges[i]);

            merge_set(edges[i].a, edges[i].b);
        }
    }
    return tree_cost;
}

int main() {

    ifstream cin ("input.txt");

    cin >> N >> M;
    for (int i = 0; i < M; i++) {
        int a, b, w; cin >> a >> b >> w;
        edges.push_back(Edge(a, b, w));
    }

    cout << kruskal() << "\n";

    return 0;
}