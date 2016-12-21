#include <bits/stdc++.h>
using namespace std;

#define MAX 1000
bool met[MAX];
vector<int> path[MAX];

pair<int, int> diameter;
void diam(int cur, int dist) {
    if (met_t[cur]) return;
    met_t[cur] = true;

    diameter = max (diameter, make_pair(dist, cur));
    for (int adj : path[cur]) {
        diam(adj, dist + 1);
    }
}

int main() {

    int N; cin >> N;
    for (int i = 1; i < N; i++) {
        int a, b; cin >> a >> b;
        path[a].push_back(b);
        path[b].push_back(a);
    }

    diam(0, 0);
    memset(met_t, 0, sizeof(met_t));
    diam(diameter.second, 0);

    return 0;
}