#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define INF 10000000000000000
#define MAX 201010
lli N, M;
lli tree_parent[MAX];
lli A[MAX][30]; // A[k][i] 2^i-th ancestor of k
pair<lli, int> maxim[MAX][30]; // the maximum value on an edge from k to 2^i-th ancestor of k
lli depth[MAX], dist[MAX];
vector<lli> tree_edge[MAX];
map<pair<lli, lli>, lli> weight, to_indx;


// ----------------- DISJOINT SET UNION -------------------------
lli parent[MAX];
lli rank_[MAX];

lli find_set (lli a) {
    if (a != parent[a])
        parent[a] = find_set(parent[a]);
    return parent[a];
}

void make_set(lli a) {
    parent[a] = a;
    rank_[a] = 0;
}

void merge_set (lli a, lli b) {
    lli root_a = find_set (a);
    lli root_b = find_set (b);

    if (root_a == root_b) return;

    if (rank_[root_a] > rank_[root_b])
        parent[root_b] = root_a;
    else
        parent[root_a] = root_b;
    if (rank_[root_a] == rank_[root_b])
        rank_[root_b]++;
}

// -------------- KRUSKAL'S ALGORITHM, MINIMAL SPANNING TREE ------------------

struct Edge {
    int indx;
    lli weight;
    lli cost;
    int a, b;
    Edge (int a, int b, lli weight, lli cost, int indx) :
        a(a),
        b(b),
        weight(weight),
        cost(cost),
        indx(indx) {}

    bool operator < (Edge const &other) const {
        return weight < other.weight;
    }
};

vector<Edge> edges, tree;
bool in_tree[MAX];

lli kruskal() {
    sort(edges.begin(), edges.end());
    for (lli i = 1; i <= N; i++)
        make_set(i);

    lli tree_cost = 0;
    for (lli i = 0; i < M; i++) {
        if (find_set(edges[i].a) != find_set(edges[i].b)) {
            tree_cost += edges[i].weight;

            tree.push_back(edges[i]);
            in_tree[edges[i].indx] = true;

            merge_set(edges[i].a, edges[i].b);
        }
    }
    return tree_cost;
}

// -------------------- LCA ------------------------------------

void dfs(lli x, lli depth_, lli dist_) {
    if (depth[x] != -1) return;
    depth[x] = depth_;
    dist[x] = dist_;
    for (lli next_ : tree_edge[x]) {
        if (tree_parent[next_] == -1)
            tree_parent[next_] = x;
        dfs(next_, depth_+1, dist_ + weight[make_pair(x, next_)]);
    }
}

void preprocessLCA(lli root) {
    for (lli i = 0; i <= N; i++) {
        for (lli k = 0; k < 30; k++) A[i][k] = -1;
        depth[i] = tree_parent[i] = -1;
    }
    dfs(root, 0, 0);
    depth[0] = 0;
    tree_parent[root] = root;

    for (lli i = 1; i <= N; i++) {
        A[i][0] = tree_parent[i];
        maxim[i][0] = make_pair(weight[make_pair(i, tree_parent[i])],
                                to_indx[make_pair(i, tree_parent[i])]);
    }

    for (lli k = 0; (1 << k) <= N; k++)
        for (lli i = 1; i <= N; i++)
            if (A[i][k] == -1) {
                A[i][k] = A[A[i][k-1]][k-1];
                maxim[i][k] = max(maxim[i][k-1], maxim[A[i][k-1]][k-1]);
            }
}

lli query(lli p, lli q) {

    if (depth[p] < depth[q])
        swap(p, q);

    lli log;
    for (log = 1; (1 << log) <= N; log++);
    log--;

    for (lli i = log; i >= 0; i--)
        if (depth[p] - (1 << i) >= depth[q])
            p = A[p][i];

    if (p == q) return p;

    for (lli i = log; i >= 0; i--)
        if (A[p][i] != -1 && A[p][i] != A[q][i])
            p = A[p][i], q = A[q][i];

    return tree_parent[p];
}

pair<lli, int> max_value_kth_ancestor(int p, int K) {
    if (K == 1) return maxim[p][0];
    else if (K == 0) return make_pair(0, 0);

    while (K > 0) {
        lli i = 0;
        for (; (1 << i) <= K; i++); i--;
        K -= (1 << i);
        return max(maxim[p][i], max_value_kth_ancestor(A[p][i], K));
    }
}

lli dissatisfaction[MAX], cost[MAX];

int main() {

    cin >> N >> M;
    for (lli i = 1; i <= M; i++)
        cin >> dissatisfaction[i];
    for (lli i = 1; i <= M; i++)
        cin >> cost[i];

    for (lli i = 1; i <= M; i++) {
        lli a, b, w; cin >> a >> b;

        if (weight[make_pair(a,b)]) {
            if (weight[make_pair(a,b)] > dissatisfaction[i]) {
                weight[make_pair(a,b)] = weight[make_pair(b,a)] = dissatisfaction[i];
                to_indx[make_pair(a,b)] = to_indx[make_pair(b,a)] = i;
            }
        } else {
            weight[make_pair(a,b)] = weight[make_pair(b,a)] = dissatisfaction[i];
            to_indx[make_pair(a,b)] = to_indx[make_pair(b,a)] = i;
        }
        edges.push_back(Edge(a, b, dissatisfaction[i], cost[i], i));
    }
    lli S; cin >> S;

    kruskal();
    lli min_cost = INF, ans = 0, total_diss = 0;
    int change_indx = -1;
    for (auto e : tree) {
        tree_edge[e.a].push_back(e.b);
        tree_edge[e.b].push_back(e.a);

        if (min_cost > e.cost) {
            min_cost = e.cost;
            change_indx = e.indx;
        }

        total_diss += dissatisfaction[e.indx];
    }
    ans = total_diss - S / min_cost;

    preprocessLCA(1);

    int erase_indx = -1;
    for (auto e : edges) {
        if (in_tree[e.indx]) continue;
        if (e.cost >= min_cost) continue;

        int lca = query(e.a, e.b);
        pair<lli, int> max_weight = max(max_value_kth_ancestor(e.a, depth[e.a] - depth[lca]),
                                        max_value_kth_ancestor(e.b, depth[e.b] - depth[lca]));

        lli cur = total_diss - max_weight.first + e.weight - S / e.cost;
        if (ans > cur) {
            ans = cur;
            change_indx = e.indx;
            erase_indx = max_weight.second;
        }
    }

    cout << ans << "\n";
    for (auto e : edges) {
        if (e.indx == erase_indx) continue;
        if (e.indx == change_indx) {
            cout << e.indx <<  " " << (e.weight - S / e.cost);
        } else if (!in_tree[e.indx]) continue;
        else
            cout << e.indx <<  " " << e.weight;
        cout << "\n";
    }

    return 0;
}