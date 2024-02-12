# Notebook

$\mathcal O(N\log N)$

$$
\texttt{dp}[l][r] =
\sum_{\substack{p\in [1, l) \\ \texttt{sum}(l, r)\geq \texttt{sum}(p, l-1)}}{\frac{1}{r-p}\cdot \texttt{dp}[p][r]} +
\sum_{\substack{p\in (r, N] \\ \texttt{sum}(l, r) > \texttt{sum}(r+1, p)}}{\frac{1}{p-l}\cdot \texttt{dp}[l][p]}
$$
