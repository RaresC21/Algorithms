/*
  USACO Gold 2014 December.
  Problem: Guard Mark.
  Bitmask DP using recurision, memoization.
*/

#include <bits/stdc++.h>
using namespace std;

typedef long long int lli;

#define INF 1000000000000000
#define MAX (1 << 20) + 1
struct Cow {
    lli h, w, s;
    Cow () {}
    Cow (lli  h, lli  w, lli  s) : h(h), w(w), s(s) {}
};

lli N, H, ans = -1;
vector<Cow> cows;
lli state_height[MAX];
lli dp[MAX], heigth[MAX]; // for a state, max amount we can add on.

// start at state 11...11
lli  DP (lli state) {

    if (state_height[state] == 0) {
        lli height = 0;
        for (int i = 0; i < N; i++) {
            if (state & (1 << i)) {
                height += cows[i].h;
            }
        }
        state_height[state] = height;
    }

    if (state == 0) {
        return dp[state] = INF;
    }

    if (dp[state] != -1)
        return dp[state];

    lli best = -INF;
    for (int i = 0; i < N; i++) {
        lli  cur = state & (1 << i);
        if (cur == 0) continue;

        lli new_state = state ^ (1 << i);
        lli sub_strength = DP(new_state);

        lli cur_strength = min (sub_strength - cows[i].w, cows[i].s);
        best = max(best, cur_strength);
    }

    if (state_height[state] >= H && best >= 0) {
        ans = max (ans, best);
    }

    return dp[state] = best;
}

int main() {

    ifstream cin ("guard.in");
    ofstream cout ("guard.out");

    cin >> N >> H;
    for (int i = 0; i < N; i++) {
        lli h, w, s; cin >> h >> w >> s;
        cows.push_back(Cow(h, w, s));
    }

    memset (dp, -1, sizeof(dp));
    lli start_ = (1 << N) - 1;
    DP(start_);

    if (ans < 0) {
        cout << "Mark is too tall\n";
    } else {
        cout << ans << "\n";
    }

    return 0;
}
