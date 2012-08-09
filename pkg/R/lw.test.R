# Ledoit-Wolf sphericity test for high-dimensional data
lw.test = function(X) {
	n = nrow(X) - 1   # Fujikoshi et al., p. 219
	p = ncol(X)
	S = cov(X)
	trS = sum(diag(S))
	trS2 = sum(S^2)
	stat = (n*p*trS2/trS^2 - n - p - 1) / 2
	pvalue = 2 * pnorm(-abs(stat))
	list(stat=stat, pvalue=pvalue)
}
