// infoarena.com - Cel mai lung subsir comun

#include <bits/stdc++.h>
using namespace std;

#define MAX 1100
int A[MAX], B[MAX];
int dp[MAX][MAX];

int main() {

    ifstream cin ("cmlsc.in");
    ofstream cout ("cmlsc.out");

    int N, M; cin >> N >> M;
    for (int i = 1; i <= N; i++) {
        cin >> A[i];
    }

    for (int i = 1; i <= M; i++) {
        cin >> B[i];
    }


    int ans = 0;
    pair<int, int> aa;
    for (int i = 1; i <= N; i++) {
        for (int k = 1; k <= M; k++) {
            dp[i][k] = max (dp[i - 1][k], dp[i][k - 1]);
            if (dp[i][k] < dp[i - 1][k]) {
                dp[i][k] = dp[i - 1][k];
            }
            if (dp[i][k] < dp[i][k - 1]) {
                dp[i][k] = dp[i][k - 1];
            }

            if (A[i] == B[k]) {
                if (dp[i][k] < dp[i - 1][k - 1] + 1) {
                    dp[i][k] = dp[i - 1][k - 1] + 1;
                }
            }

            if (dp[i][k] > ans) {
                ans = dp[i][k];
                aa = {i, k};
            }
        }
    }

    cout << ans << "\n";
    int xx = aa.first;
    int yy = aa.second;

    vector<int> sol;
    while (xx >= 0 && yy >= 0 && dp[xx][yy] != 0) {
        if (A[xx] == B[yy]) {
            sol.push_back(A[xx]);
            --xx; --yy;
        }
        else if (dp[xx - 1][yy] < dp[xx][yy - 1]) {
            --yy;
        } else {
            --xx;
        }
    }

    for (int i = sol.size() - 1; i >= 0; i--) {
        cout << sol[i] << " ";
    }
    cout << "\n";

    return 0;
}
