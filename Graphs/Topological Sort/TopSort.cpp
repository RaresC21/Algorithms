#include <vector>

#define pb push_back
#define MAX 100

vector<int> topsort;
vector<int> edge[MAX];
bool met[MAX];

void dfs(int node) {
    met[node] = true;
    for (int adj : edge[node]) {
        if (!met[adj])
            dfs(adj);
    }
    topsort.pb(node);
}