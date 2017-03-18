/*
    code golf, SWERC 2014.
    https://open.kattis.com/problems/golfbot
*/
#include <iostream>
#include <algorithm>
#include <queue>
#include <set>
#include <vector>
#include <complex>

#define PA pair<int,int>
#define PA2 pair<PA, int>
#define MAX_N 524288
#define VCD vector<complex<double> >

using namespace std;

int n, m, k[MAX_N], d[MAX_N], ans;
bool is[MAX_N];

const double PI = acos(-1);
VCD FFT(VCD a, int t) {
    if (a.size() == 1) return a;
    VCD x, y;
    for (int i = 0 ; i < a.size() ; i ++) {
        if (i&1) y.push_back(a[i]); else x.push_back(a[i]);
    }
    x = FFT(x,t);
    y = FFT(y,t);
    int n = a.size();
    double arc = 2 * PI * t / n;
    complex<double> wn(cos(arc), sin(arc)), w = 1;
    for (int i = 0 ; i < x.size() ; ++i) {
        a[i] = x[i] + w * y[i];
        a[i + n/2] = x[i] - w * y[i];
        w = w * wn;
    }
    return a;
}

int main() {
    cin >> n;
    for (int i = 0 ; i < n ; i ++) {
        cin >> k[i];
        is[k[i]] = 1;
    }
    cin >> m;
    for (int i = 0 ; i < m ; i ++) cin >> d[i];
    VCD X, Y;
    for (int i = 0 ; i < MAX_N ; i ++) {
        X.push_back(is[i]);
        Y.push_back(is[i]);
    }
    X = FFT(X, 1);
    Y = FFT(Y, 1);
    VCD C;
    for (int i = 0 ; i < MAX_N ; i ++) C.push_back(X[i] * Y[i]);
    C = FFT(C, -1);
    ans = 0;

    for (int i = 0 ; i < m ; i ++) {
        if (C[d[i]].real() > 0.5 || is[d[i]]) ans ++;
    }
    cout << ans << endl;
}
