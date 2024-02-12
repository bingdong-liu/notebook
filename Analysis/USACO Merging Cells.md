Note: let $\texttt{sum}(l, r)$ denote $s_l + s_{l+1} + \dots + s_{r-1} + s_r$.

**Subtask 1 ($N\leq 8$)**

It suffices to brute force all $(N-1)!$ permutations of merges.

**Subtask 2 ($N\leq 100$)**

Let's try to calculate the answer for a given cell $c$. At any step, $c$ must be
the label of a cell that represents the sum of a range $[l, r]$ in the original
array. Let $\texttt{dp}[l][r]$ represent the probability that $c$ becomes the
final label if we only consider the cells in the range $[l, r]$. For
transitions, we can consider enumerating the last merge while making sure that
the $c$ remains as the label. Also, note that each of the $r-l$ possible merges
has an equal probability of being the last merge.

The base case is $\texttt{dp}[c][c] = 1$. Otherwise, we have the following:

$$
\texttt{dp}[l][r] =
\sum_{\substack{m\in (l, c] \\ \texttt{sum}(l, m-1)\leq \texttt{sum}(m, r)}}{\frac{1}{r-l}\cdot\texttt{dp}[m][r]} +
\sum_{\substack{m\in [c, r] \\ \texttt{sum}(l, m) > \texttt{sum}(m+1, r)}}{\frac{1}{r-l}\cdot\texttt{dp}[l][m]}
$$

This runs in $\mathcal O(N^3)$ for a single cell $c$, resulting in a
$\mathcal O(N^4)$ solution.

```cpp
#include <bits/stdc++.h>
using namespace std;

const int MOD = 1e9 + 7;

int main() {
    int N; cin >> N;
    vector S(N + 1, 0);
    for (int i = 0; i < N; i++) {
        int x; cin >> x;
        S[i + 1] = S[i] + x;
    }
    auto sum = [&](int l, int r) {
        return S.at(r + 1) - S.at(l);
    };
    vector inv(N, 0);
    inv[1] = 1;
    for (int i = 2; i <= N; i++)
        inv[i] = (long long) inv[MOD % i] * (MOD - MOD / i) % MOD;
    for (int c = 0; c < N; c++) {
        vector dp(N, vector(N, 0));
        dp[c][c] = 1;
        for (int l = c; l >= 0; l--) for (int r = c; r < N; r++) {
            if (l == c && r == c) continue;
            long long x = 0;
            for (int m = l + 1; m <= c; m++) // last merge was leftwards
                if (sum(m, r) >= sum(l, m - 1)) x += dp[m][r];
            for (int m = c; m < r; m++) // last merge was rightwards
                if (sum(l, m) > sum(m + 1, r)) x += dp[l][m]; 
            dp[l][r] = x % MOD * inv[r - l] % MOD; // probability of fixed last merge
        }
        cout << dp[0][N - 1] << '\n';
    }
}
```

**Subtask 3 ($N\leq 500$)**

It is possible to optimize the solution from subtask 2 to $\mathcal O(N^3)$
using prefix sums and two pointers (described more in detail in subtask 4's
notes).

However, I will discuss another approach that is closer to the full solution.
Notice that recalculating the DP for each $c$ can be suboptimal in many ways;
for example, regardless of which label a cell $(l, r)$ contains, it has the same
contribution to $\texttt{dp}[1][N]$. Therefore, let's consider calculating the
DP backwards.

Let $\texttt{dp}[l][r]$ represent the probability that the label of the cell
representing the range $[l, r]$ becomes the label of the final cell after all
merges. The base case is $\texttt{dp}[1][N] = 1$. Otherwise, transitions are
similar to the ones described in the subtask 2 solution:

$$
\texttt{dp}[l][r] =
\sum_{\substack{p\in [1, l) \\ \texttt{sum}(l, r)\geq \texttt{sum}(p, l-1)}}{\frac{1}{r-p}\cdot \texttt{dp}[p][r]} +
\sum_{\substack{p\in (r, N] \\ \texttt{sum}(l, r) > \texttt{sum}(r+1, p)}}{\frac{1}{p-l}\cdot \texttt{dp}[l][p]}
$$

The answer for a cell $c$ is $\texttt{dp}[c][c]$. This runs in
$\mathcal O(N^3)$.

```cpp
#include <bits/stdc++.h>
using namespace std;

const int MOD = 1e9 + 7;

int main() {
    int N; cin >> N;
    vector S(N + 1, 0);
    for (int i = 0; i < N; i++) {
        int x; cin >> x;
        S[i + 1] = S[i] + x;
    }
    auto sum = [&](int l, int r) {
        return S.at(r + 1) - S.at(l);
    };
    vector inv(N, 0);
    inv[1] = 1;
    for (int i = 2; i <= N; i++)
        inv[i] = (long long) inv[MOD % i] * (MOD - MOD / i) % MOD;
    vector dp(N, vector(N, 0));
    dp[0][N - 1] = inv[N - 1];
    for (int l = 0; l < N; l++) for (int r = N - 1; r >= l; r--) {
        if (l == 0 && r == N - 1) continue;
        long long x = 0;
        for (int p = 0; p < l; p++) // last merge was leftwards
            if (sum(l, r) >= sum(p, l - 1)) x += dp[p][r];
        for (int p = r + 1; p < N; p++) // last merge was rightwards
            if (sum(l, r) > sum(r + 1, p)) x += dp[l][p];
        dp[l][r] = x % MOD;
        if (l < r) // probability of fixed last merge
            dp[l][r] = (long long) dp[l][r] * inv[r - l] % MOD;
    }
    for (int c = 0; c < N; c++) cout << dp[c][c] << '\n';
}
```

**Subtask 4 $(N\leq 5000$)**

Let's try to accelerate the solution from subtask 3. For example, say the last
merge is leftwards. We query the sum of $\frac{\texttt{dp}[p][r]}{r-p}$ over a
fixed $r$ and a range of $p$. The range of $p$ is determined by the leftmost
position such that $\texttt{sum}(l, r)\geq \texttt{sum}(p, l-1)$ holds. Let this
bound be denoted as $\texttt{pL}[l][r]$. It can be proven that
$\texttt{pL}[l][r]\leq \texttt{pL}[l+1][r]$, which allows us to compute all
$\texttt{pL}[l][r]$ in $\mathcal O(N^2)$ time using two pointers. Next, to
calculate the sum, we need a data-structure that can support adding
$\frac{\texttt{dp}[l][r]}{r-l}$ and querying the sum over
$p\in [\texttt{pL}[l][r], l)$. Since we enumerate ranges in decreasing order of
size, for a fixed $r$, $\frac{\texttt{dp}[l][r]}{r-l}$ is added in increasing
order of $l$. This allows us to maintain prefix sums and answer range sum
queries in constant time.

This runs in $\mathcal O(N^2)$.

```cpp
#include <bits/stdc++.h>
using namespace std;

#define sz(v) int(std::size(v))
const int MOD = 1e9 + 7;

int add(int a, int b) { return a + b >= MOD ? a + b - MOD : a + b; }
int sub(int a, int b) { return add(a, MOD - b); }
int mul(int a, int b) { return (long long) a * b % MOD; }
void add_self(int &a, int b) { a = add(a, b); }

struct sumL { // fixed l => updates (r) in decreasing order
    vector<int> t; sumL(int N) : t(N) {}
    void set(int i, int x) { t.at(i) = add(x, i + 1 < sz(t) ? t.at(i + 1) : 0); }
    int sum(int l, int r) { return l > r ? 0 : sub(t.at(l), r + 1 < sz(t) ? t.at(r + 1) : 0); }
};
struct sumR { // fixed r => updates (l) in increasing order
    vector<int> t; sumR(int N) : t(N) {}
    void set(int i, int x) { t.at(i) = add(x, i ? t.at(i - 1) : 0); }
    int sum(int l, int r) { return l > r ? 0 : sub(t.at(r), l ? t.at(l - 1) : 0); }
};

int main() {
    int N; cin >> N;
    vector<int> S(N + 1);
    for (int i = 0; i < N; i++) {
        int x; cin >> x;
        S[i + 1] = S[i] + x;
    }
    auto sum = [&](int l, int r) {
        return S.at(r + 1) - S.at(l);
    };
    vector<int> inv(N);
    inv[1] = 1;
    for (int i = 2; i < N; i++) inv[i] = mul(inv[MOD % i], MOD - MOD / i);
    vector dp(N, vector<int>(N));
    vector sL(N, sumL(N)); vector sR(N, sumR(N));
    // pL[l][r] = farthest we can merge leftwards such that sum(l, r) >= sum(p, l-1)
    // pR[l][r] = farthest we can merge rightwards such that sum(l, r) > sum(r+1, p)
    vector pL(N, vector<int>(N)), pR(N, vector<int>(N));
    for (int l = 0; l < N; l++) for (int r = l, p = l; r < N; r++) {
        while (p > 0 && sum(l, r) >= sum(p - 1, l - 1)) p--;
        pL[l][r] = p;
    }
    for (int r = 0; r < N; r++) for (int l = r, p = r; l >= 0; l--) {
        while (p + 1 < N && sum(l, r) > sum(r + 1, p + 1)) p++;
        pR[l][r] = p;
    }
    dp[0][N - 1] = inv.at(N - 1);
    sL[0].set(N - 1, dp[0][N - 1]), sR[N - 1].set(0, dp[0][N - 1]);
    for (int l = 0; l < N; l++) for (int r = N - 1; r >= l; r--) {
        if (l == 0 && r == N - 1) continue;
        add_self(dp[l][r], sR[r].sum(pL[l][r], l - 1)); // last merge was leftwards
        add_self(dp[l][r], sL[l].sum(r + 1, pR[l][r])); // last merge was rightwards
        if (l < r) dp[l][r] = mul(dp[l][r], inv.at(r - l)); // probability of fixed last merge
        sL[l].set(r, dp[l][r]), sR[r].set(l, dp[l][r]);
    }
    for (int c = 0; c < N; c++) cout << dp[c][c] << '\n';
}
```

The code above uses around $5N^2$ memory. It's possible to optimize this to
$N^2$ by just storing the DP array while maintaining $\texttt{pL}[l][r]$,
$\texttt{pR}[l][r]$, and the range sums in 1D arrays. 

Benjamin Qi's code:

```cpp
#include <bits/stdc++.h>
using namespace std;
 
const int MOD = 1e9 + 7;
 
int main() {
    int N; cin >> N;
    vector<int> S(N);
    for (auto &x : S) cin >> x;
    vector<int> cum_S{0};
    for (auto t : S) cum_S.push_back(cum_S.back() + t);
    vector win_prob(N + 1, vector<int>(N + 1));
    vector<int> invs(N);
    invs[1] = 1;
    for (int i = 2; i < N; i++) invs[i] = (long long) invs[MOD % i] * (MOD - MOD / i) % MOD;
    win_prob.at(0).at(N) = 1;
    vector<int> l_sum(N + 1), r_sum(N + 1), l_ptr(N + 1), r_ptr(N + 1, N);
    for (int l = 0; l <= N; l++) for (int r = N; r > l; r--) {
        while (!(cum_S.at(l) - cum_S.at(l_ptr[r]) <= cum_S.at(r) - cum_S.at(l))) {
            l_sum[r] -= win_prob.at(l_ptr[r]).at(r);
            if (l_sum[r] < 0) l_sum[r] += MOD;
            ++l_ptr[r];
        }
        while (!(cum_S.at(r) - cum_S.at(l) > cum_S.at(r_ptr[l]) - cum_S.at(r))) {
            r_sum.at(l) -= win_prob.at(l).at(r_ptr[l]);
            if (r_sum[l] < 0) r_sum[l] += MOD;
            --r_ptr[l];
        }
        win_prob.at(l).at(r) += l_sum.at(r);
        win_prob.at(l).at(r) += r_sum.at(l);
        if (win_prob[l][r] >= MOD) win_prob[l][r] -= MOD;
        if (r > l + 1) win_prob.at(l).at(r) = (long long) win_prob[l][r] * invs.at(r - l - 1) % MOD;
        l_sum[r] += win_prob.at(l).at(r);
        if (l_sum[r] >= MOD) l_sum[r] -= MOD;
        r_sum[l] += win_prob.at(l).at(r);
        if (r_sum[l] >= MOD) r_sum[l] -= MOD;
    }
    for (int i = 0; i < N; i++) cout << win_prob.at(i).at(i + 1) << '\n';
}
```
