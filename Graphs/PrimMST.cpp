#include <bits/stdc++.h>
using namespace std;

#define INF 1000000000
#define MAX 1011
int N; // number of nodes.
int weight[MAX][MAX]; // adjacency matrix denoting weight between nodes

bool intree[MAX]; // is node in tree
vector<int> neighbors[MAX]; // vector of nodes that are neighbors of node i

struct Node {
    int n;
    int dist;
    Node (int n, int dist) : n(n), dist(dist) {}
    bool operator < (const Node &other) const {
        return dist > other.dist;
    }
};

void init(int source) {
    for (int i = 1; i <= N; i++) {
        intree[i] = false;
    }
}

int prim(int source) {

    init(source);
    int tree_size = 0;
    int tree_cost = 0; // sum of weights of edges in MST.

    priority_queue<Node> q; // priority queue holding distances from tree to
                                // neighbor nodes
    q.push(Node(source, 0));

    while (tree_size < N) {

        // find node not in tree that is closes to tree.
        Node p = q.top();
        q.pop();

        if (intree[p.n]) continue;
        intree[p.n] = true;

        tree_size ++;
        tree_cost += p.dist;

        for (int i = 0; i < neighbors[p.n].size(); i++) {
            int cur_node = neighbors[p.n][i];
            if (intree[cur_node]) continue;
            q.push(Node(cur_node, weight[p.n][cur_node]));
        }
    }
    return tree_cost;
}

int main() {
    cin >> N;
    for (int i = 1; i <= N; i++) {
        for (int k = 1; k <= N; k++) {
            int x; cin >> x;
            weight[i][k] = x;
            if (x != 0) {
                neighbors[i].push_back(k);
            }
        }
    }

    cout << prim(1) << "\n";

    return 0;
}
