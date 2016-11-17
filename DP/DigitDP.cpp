/*
    2014 Southeast Regionals ACM ICPC : Hill Number.

    A Hill Number is a positive integer, the digits of which possibly rise and
    then possibly fall, but never fall and then rise.

    Given a positive integer, if it is a hill number, print the number of
    positive hill numbers less than or equal to it. If it is not a hill number,
    print -1.
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli N;
int digs[20]; // digits
lli dp[20][10][2][2];

lli DP(int pos, int cur, int rf, bool smaller, int len) {
    lli d = dp[pos][cur][rf][smaller];
    if (d != -1)
        return d;

    d = 0;

    if (pos == len) {
        if (smaller)
            return d = 1;
    }
    else {
        for (int i = 0; i <= (smaller ? 9 : digs[pos]); i++) {
            bool rise_fall = rf;

            // it was falling, but now rising - no good.
            if (rf == 0 && i > cur)
                break;

            if (i < cur)
                rise_fall = 0;

            bool s = smaller || (i < digs[pos]);
            d += DP(pos + 1, i, rise_fall, s, len);
        }
    }
    return dp[pos][cur][rf][smaller] = d;
}

bool hill_(int len) {
    bool rf = 1;
    for (int i = 1; i < len; i++) {
        if (rf == 0 && digs[i] > digs[i - 1]) return false;
        if (digs[i] < digs[i - 1]) rf = 0;
    }
    return true;
}

void get_digits(lli N) {
    int i = log10(N);
    int len = i + 1;
    while (N > 0) {
        digs[i] = N % 10;
        N /= 10;
        i--;
    }
}

int main() {
    int T; cin >> T;
    for (int t = 0; t < T; t++) {

        cin >> N;
        memset(dp, -1, sizeof(dp));

        get_digits(N);
        if (hill_(log10(N) + 1)) {
            get_digits(N + 1);
            cout << DP(0, 0, 1, false, log10(N+1) + 1) - 1 << "\n";
        }
        else {
            cout << "-1\n";
        }
    }

    return 0;
}