#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

#define pb push_back
#define mp make_pair
#define MOD 100000
#define MAX 201


int main() {

    int S; cin >> S;
    int a, b, c; cin >> a >> b >> c;
    int sum = a + b + c;

    double x = (S * a) * 1.0 / sum * 1.0;
    double y = (S * b) * 1.0 / sum * 1.0;
    double z = (S * c) * 1.0 / sum * 1.0;

    if (sum == 0) {
        cout << "0 0 0\n";
        return 0;
    }

    cout << fixed << setprecision(15) << x <<  " " << y << " " << z << "\n";

    return 0;
}