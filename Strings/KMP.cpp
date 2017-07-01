#include <bits/stdc++.h>
using namespace std;

#define MAX 100000
int lps[MAX];

void longest_prefix_suffix(const string& pattern) {
    lps[0] = lps[1] = 0;
    int m = pattern.size();

    for (int i = 2; i <= m; i++) {
        int j = lps[i - 1];
        while (true) {
            if (pattern[j] == pattern[i - 1]) {
                lps[i] = j + 1;
                break;
            }
            if (j == 0) {
                lps[i] = 0;
                break;
            }
            j = lps[j];
        }
    }
}

vector<int> kmp(const string& pattern, const string& search) {
    vector<int> ans;
    longest_prefix_suffix(pattern);

    int i = 0, j = 0;
    while (i < text.size()) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }
        if (j == pattern.size()) {
            ans.push_back(i - j);
            j = lps[j - 1];
        } else if (i < text.size() && pattern[j] != text[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i++;
        }
    }
    return ans;
}