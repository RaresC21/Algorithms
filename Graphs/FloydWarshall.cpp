#include <bits/stdc++.h>
using namespace std;

#define INF 1000000000
#define MAX 401
int dist[MAX][MAX];

void init() {
    for (int i = 0; i < MAX; i++)
        for (int k = 0; k < MAX; k++)
            if (i != k)
                dist[i][k] = INF;
}

void floyd_warshall(int N) {
    for(int k = 1; k <= N; k++)
        for(int i = 1; i <= N; i++)
            for(int j = 1; j <= N; j++)
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
}

int main() {

    init();

    int N, M; cin >> N >> M; // number of nodes, edges respectively.
    for (int i = 0; i < M; i++) {
        int a, b, w; cin >> a >> b >> w;
        dist[a][b] = w; // in this case we have directed edges.
    }

    floyd_warshall(N);

    int Q; cin >> Q;
    for (int q = 0; q < Q; q++) {
        int a, b; cin >> a >> b;
        cout << ((dist[a][b] != INF) ? dist[a][b] : -1) << "\n";
    }

    return 0;
}