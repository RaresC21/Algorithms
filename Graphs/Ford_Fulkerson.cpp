/*
Determines the maximum flow through a directed graph from source to sink.
Implementation of Ford-Fulkerson Algorithm
*/

#include <bits/stdc++.h>

class MaxFlow {

#define INF 1000000000
#define MAX 201

	lli N, M;
	bool visited[MAX];
	lli flow[MAX], prevnode[MAX];
	lli capacity[MAX][MAX];

	lli maximum_flow(lli source, lli sink)
	{
		lli total_flow = 0;
		while (true)
		{
			// find path with highest capacity from source to sink.
			memset(flow, 0, sizeof(flow));
			memset(visited, 0, sizeof(visited));
			memset(prevnode, 0, sizeof(prevnode));
			flow[source] = INF;

			lli max_loc, max_flow;

			while (true)
			{
				max_flow = 0;
				max_loc = 0;

				// find the unvisited node with the highest capacity.
				for (lli i = 1; i <= M; i++)
				{
					if (flow[i] > max_flow && !visited[i])
					{
						max_flow = flow[i];
						max_loc = i;
					}
				}

				if (max_flow == 0 || max_loc == sink)
					break;
				visited[max_loc] = true;

				// update the neighbors.
				for (lli i = 1; i <= M; i++)
				{
					if (capacity[max_loc][i] == 0) continue;
					if (flow[i] < min(max_flow, capacity[max_loc][i])) {
						flow[i] = min(max_flow, capacity[max_loc][i]);
						prevnode[i] = max_loc;
					}
				}
			}

			if (max_loc == 0)
				break;

			lli path_capacity = flow[sink];
			total_flow += path_capacity;

			// add that flow to the network. Update capacity appropriately.
			lli cur_node = sink;
			while (cur_node)
			{
				lli next_node = prevnode[cur_node];
				capacity[next_node][cur_node] -= path_capacity;
				capacity[cur_node][next_node] += path_capacity;
				cur_node = next_node;
			}
		}
		return total_flow;
	}

	int main() {

		cin >> N >> M;
		for (lli i = 0; i < N; i++) {
			lli s, e, c; cin >> s >> e >> c;
			capacity[s][e] += c;
		}

		cout << maximum_flow(1, M) << "\n";

		return 0;
	}
};
