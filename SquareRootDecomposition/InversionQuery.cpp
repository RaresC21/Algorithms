/*
hackerearth.com
Sherlock and Inversions.

Watson gives to Sherlock an array of N integers denoted by A1, A2 ... AN.
Now he gives him Q queries of form Li, Ri. For each such query Sherlock has to
report the number of inversions in subarray denoted by [Li, Ri].
*/

#include <bits/stdc++.h>
using namespace std;
#define int long long

#define MAX 100005

class BitTree {
public:

    int N;
    int BIT[MAX];

    BitTree() {}

    BitTree(int N) : N(N) {
        memset(BIT, 0, sizeof(BIT));
    }

    void update(int x, int val) {
        for(; x <= N; x += x & -x)
            BIT[x] += val;
    }

    int query(int x) {
        int sum = 0;
        for(; x > 0; x -= x & -x)
            sum += BIT[x];
        return sum;
    }
};

#define BLOCK 200
int N, M;
int arr[MAX];
int amnt[MAX];

BitTree *bit;

// handle square root decomposition in sorting. ----------------------------
inline bool cmp (const pair<int, int>& a, const pair<int, int>& b) {
    if (a.first / BLOCK != b.first / BLOCK)
        return a.first/BLOCK < b.first/BLOCK;
    return a.second < b.second;
}

void add(int k, int& ans, int sgn) {
    bit->update(arr[k], 1);
    ans -= amnt[k];
    if (sgn == -1)
        amnt[k] = bit->query(N) - bit->query(arr[k]);
    else
        amnt[k] = bit->query(arr[k] - 1);
    ans += amnt[k];
}

void remove_(int k, int& ans, int sgn) {
    bit->update(arr[k], -1);
    if (sgn == 1)
        ans -= bit->query(arr[k] - 1);
    else
        ans -= bit->query(N) - bit->query(arr[k]);
    amnt[k] = 0;
}

void compress(set<int> &nums) {
    int cur = 1;
    map<int, int> compression;
    for (auto n : nums) {
        compression[n] = cur;
        cur++;
    }

    for (int i = 1; i <= N; i++)
        arr[i] = compression[arr[i]];
}

#undef int
int main() {
#define int long long

    ios_base::sync_with_stdio(false);

    cin >> N;
    set<int> nums;
    for (int i = 1; i <= N; i++) {
        cin >> arr[i];
        nums.insert(arr[i]);
    }
    compress(nums);

    bit = new BitTree(N);

    vector<pair<int, int> > v, orig;
    cin >> M;
    for (int q = 0; q < M; q++) {
        int a, b; cin >> a >> b;
        v.push_back(make_pair(a, b));
    }

    orig = v;
    sort(v.begin(), v.end(), cmp);

    int l = 1, r = 0;
    int ans = 0;
    map<pair<int, int>, int> sol;

    for (int i = 0; i < M; i++) {
        int a = v[i].first, b = v[i].second;

        while (r < b) {
            r++;
            add(r, ans, -1);
        }

        while (l > a) {
            l--;
            add(l, ans, 1);
        }

        while (r > b) {
            remove_(r, ans, -1);
            r--;
        }

        while (l < a) {
            remove_(l, ans, 1);
            l++;
        }

        l = a;
        r = b;

        sol[make_pair(a, b)] = ans;
    }

    for (int i = 0; i < M; i++) {
        cout << sol[orig[i]] << "\n";
    }

    delete bit;
    return 0;
}