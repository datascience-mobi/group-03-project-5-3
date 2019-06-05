# plot(density(log10(coverage_all$Coverage)), main = "Distribution of coverage values")
# We set the lower threshold for the coverage by 30, based on literature,
# and on distribution of coverage in our dataset. We will loose 4 % of our dataset this way.
# Function for drawing quantiles of istribution (log scale plot): P = seq(0,1, 0.01); qq = quantile(g_AML.covmean$rowMeans.g_AMLpat.cov., P);  qq10 = log10(qq); abline(v = qq10, col = "red")

P = seq(0,1, 0.001)
Q = ecdf(g_AML.covmean$rowMeans.g_AMLpat.cov.) (P)
Y = seq(0, 200000, 1)
plot(Y, Q, type = "n", main = "quantiles for coverage values")
lines(Y, Q)