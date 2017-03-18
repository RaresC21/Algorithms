/*

Given an array of size N. All elements of array <= N. You need to answer
M queries. Each query is of the form L, R. You need to answer the count of
distinct values in range [L, R]
*/

#include <bits/stdc++.h>
using namespace std;

#define MAX 30010
#define BLOCK 555
int arr[MAX];
vector<pair<int, int> > sqrt_decomp[MAX];
int amnt[1000011];

// handle square root decomposition in sorting. ----------------------------
inline bool cmp (const pair<int, int>& a, const pair<int, int>& b) {
    if (a.first / BLOCK != b.first / BLOCK)
        return a.first/BLOCK < b.first/BLOCK;
    return a.second < b.second;
}

void add(int k, int& ans) {
    if (amnt[arr[k]] == 0)
        ans++;
    amnt[arr[k]]++;
}

void remove(int k, int& ans) {
    amnt[arr[k]]--;
    if (amnt[arr[k]] == 0)
        ans--;
}

int main() {

    std::ios_base::sync_with_stdio(false);

    int N; cin >> N;
    for (int i = 1; i <= N; i++) {
        cin >> arr[i];
    }

    int Q; cin >> Q;
    vector<pair<int, int> > v, orig;
    for (int q = 0; q < Q; q++) {
        int a, b; cin >> a >> b;
        v.push_back(make_pair(a, b));
    }

    orig = v;
    sort(v.begin(), v.end(), cmp);

    int l = 1, r = 0;
    int ans = 0;
    map<pair<int, int>, int> sol;
    for (int i = 0; i < Q; i++) {
        int a = v[i].first, b = v[i].second;
        while (r < b) {
            r++;
            add(r, ans);
        }

        while (l > a) {
            l--;
            add(l, ans);
        }

        while (r > b) {
            remove(r, ans);
            r--;
        }

        while (l < a) {
            remove(l, ans);
            l++;
        }

        l = a;
        r = b;

        sol[make_pair(a, b)] = ans;
    }

    for (int i = 0; i < Q; i++) {
        cout << sol[orig[i]] << "\n";
    }

    return 0;
}