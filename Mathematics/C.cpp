#include <bits/stdc++.h>
using namespace std;

double sq(double x) {
    return x*x;
}

int main() {

    cout << fixed << setprecision(10);

    int N; cin >> N;
    string s; cin >> s;

    int left_ = 0;
    for (int i = 0; i < N; i++)
        if (s[i] == 'T') left_++;
        else break;

    if (left_ == N) {
        cout << 2 * N + 1 << "\n";
        return 0;
    }

    int right_ = 0;
    for (int i = N - 1; i >= 0; i--) {
        if (s[i] == 'T') right_++;
        else break;
    }

    double ans = N; // bottom
    if (left != 0) {
        if (s[left_] == 'S') {
            ans += sqrt(sq(1 - sin(PI/3)) + sq(left*1.0 - 0.5));
        } else {
            double l = sqrt(sq(sin(PI/3) - 0.5) + sq(left_*1.0 - 0.5));
            double straight = sqrt(sq(l) - sq(0.5));
            double x = sqrt(sq(sin(PI/3)) + sq(left(1.0 - 0.5)));

            double theta = acos(0.5 / l);
            double beta = acos((sq(x) - sq(l) - sq(0.5)) / (-l));
            double alpha = PI - theta - beta;

            ans += straight + 0.5 * alpha;
        }
    }

    cout << ans << "\n";

    return 0;
}