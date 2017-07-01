/*
    Codeforces Round 290. Prob C.
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define pb push_back
#define mp make_pair
#define MOD 100000
#define MAX 201

vector<int> topsort;
vector<int> edge[MAX];
map<int, int> place;
bool met[MAX];

void dfs(int node) {
    met[node] = true;
    for (int adj : edge[node]) {
        if (!met[adj])
            dfs(adj);
    }
    topsort.pb(node);
}

string name[MAX];

int main() {

    int N; cin >> N;
    for (int i = 0; i < N; i++) {
        cin >> name[i];
    }

    for (int i = 0; i < N; i++) {
        for (int k = i + 1; k < N; k++) {
            for (int j = 0; j < min(name[i].size(), name[k].size()); j++) {
                if (name[i][j] == name[k][j]) continue;
                edge[(int)(name[i][j] - 'a')].pb((int)(name[k][j] - 'a'));
                break;
            }
        }
    }

    for (int i = 0; i < 26; i++) {
        if (!met[i]) {
            dfs(i);
        }
    }

    for (int i = 0; i < 26; i++) {
        if (!met[i]) {
            topsort.pb(i);
        }
    }

    reverse(topsort.begin(), topsort.end());
    for (int i = 0; i < topsort.size(); i++) {
        place[topsort[i]] = i;
    }


    // check if correct

    vector<string> sorted;
    for (int i = 0; i < N; i++) {
        string cur = "";
        for (int k = 0; k < name[i].size(); k++) {
            cur = cur + (char)(place[(int)((int)name[i][k] - 'a')] + 'a');
        }
        sorted.pb(cur);
    }

    vector<string> same = sorted;
    sort(sorted.begin(), sorted.end()); // should be no change
    for (int i = 0; i < N; i++) {
        if (sorted[i] != same[i]) {
            cout << "Impossible\n";
            return 0;
        }
    }

    for (int x : topsort) {
        cout << (char)(x + 'a');
    } cout << "\n";

    return 0;
}