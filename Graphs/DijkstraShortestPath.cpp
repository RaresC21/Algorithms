/*
Find single source shortest path from source node to all other nodes.
Uses C++ standard library implementation of priority queue.
Runs in O(NLogN) time.
*/

#include <bits/stdc++.h>

class Dijkstra {

#define MAX 1000

private:

	lli parent[MAX];
	vector<vector<int> > path; // stores which edges exist.
	vector<vector<int> > graph; // stores values of edges.

	struct Node {
		int n;
		int parent;
		lli dist;
		Node(int n, int p, lli dist) : n(n), dist(dist), parent(p) {}
		bool operator < (Node const& other) const {
			return dist > other.dist;
		}
	};

public:

	Dijkstra(vector<vector<int> > const& P, vector<vector<int> > const& G) : path(P), graph(G) {
		memset(parent, 0, sizeof(parent));
	}

	lli dijkstra(lli s, int t)
	{
		bool met[MAX];
		for (int i = 0; i < MAX; i++) met[i] = false;

		priority_queue<Node> q;
		q.push(Node(s, s, 0));

		while (q.size()) {
			Node p = q.top();
			q.pop();

			if (!met[p.n]) {
				met[p.n] = true;

				parent[p.n] = p.parent;
				if (p.n == t)
					return p.dist;

				for (int i = 0; i < path[p.n].size(); i++) {
					int cur = path[p.n][i];
					q.push(Node(cur, p.n, p.dist + graph[p.n][cur]));
				}
			}
		}

		return LONG_INF;
	}

};
