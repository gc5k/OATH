set.seed(2016)

M = 3 #covariates
N = 200 #sample size
COL = c("red", "pink", "blue", "green")
#simu cov and phe
X = matrix(rnorm(M * N), nrow = N, ncol = M) #cov
Y = matrix(rnorm(N), nrow = N, ncol = 1) #phenotype
#standardization
X = apply(X, MARGIN = 2, scale)
Y = scale(Y)
mn=c("s={1}", "s={2}", "s={3}", "s={1,2}", "s={1,3}", "s={2,3}", "s={1,2,3}")

#simu loci
LC = 1
frq = runif(LC, 0, 0.5)
mk = matrix(0, nrow = N, ncol = LC)
for (i in 1:LC)
{
  mk[,i] = rbinom(N, 2, frq[i])
}

layout(matrix(c(1:8), 2, 4))
par(mai=c(0.75,0.75,0.3,0.3), ps=10)

#OATH vs sgwas
for (IDX in 1:LC)
{
  #make Phi matrix (Naive summary statistics)
  Phi = cov(cbind(Y,mk[,IDX], X))
  varMK = var(mk[,IDX])

  #locus-specific nss
  mmod = lm(Y ~ mk[,IDX])

  #A
  b_ = c(summary(mmod)$coefficient[2,1],
        Phi[1,3] / Phi[3,3], Phi[1,4] / Phi[4,4], Phi[1,5] / Phi[5,5])

  #Omega matrix
  Omega = Phi[-1, -1]
  #Lambda matrix
  Lambda = diag(diag(Omega), nrow=nrow(Omega), ncol=ncol(Omega))

  cnt=1
  for (i in 1:M)
  {
    md = combn(M, i)
    for (j in 1:ncol(md))
    {
      xx = X[, md[,j]]

      #individual-level data model
      smod = lm(Y ~ mk[,IDX] + xx)
      Sreg = summary(smod)$coefficients[2:(length(md[,j]) + 2),1] #partial reg coef
      Ssd = summary(smod)$coefficients[2:(length(md[,j]) + 2),2] #partial reg coef sd

      #joint estimate via OATH
      Lambda_s = Omega[c(1, md[,j]+1), c(1, md[,j]+1)]
      Omega_s = Lambda[c(1, md[,j]+1), c(1, md[,j]+1)]
      b_s = b_[c(1, md[,j]+1)]
      B = solve(Lambda_s) %*% Omega_s %*% b_s
      SigmaB = sqrt(diag((Phi[1,1] - t(B) %*% Omega_s %*% b_s)[1,1] * solve(Lambda_s) / (N - nrow(Lambda_s))))

      ####data illustration
      #compare coef
      plot(B, Sreg, xlab = "OATH", ylab = "Individal-level data", bty =
          'n', main = mn[cnt], col=COL[c(1, md[,j]+1)], pch =
          16, cex = 2, xlim=c(-0.2, 0.2), ylim=c(-0.2, 0.2))

      #compare coef sd
      for (k in 1:nrow(B))
      {
        lines(x = c(B[k,1] - SigmaB[k], B[k,1] + SigmaB[k]), y = c(B[k,1], B[k,1]))
      }
      for (k in 1:nrow(B))
      {
        lines(y = c(Sreg[k] - Ssd[k], Sreg[k] + Ssd[k]), x = c(Sreg[k], Sreg[k]))
      }

      #reference line
      abline(a = 0, b = 1, lty = 2)
      cnt=cnt+1
    }
  }
}
plot.new()
legend("topleft", legend = c("Locus", "Cov 1", "Cov 2", "Cov 3"), col = COL, pch =16, bty = 'n')
