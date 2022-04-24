ours_general = function(kernel, kappa=4){
  
  # for any general kernel matrix, output results from our method
  
  library(psych)
  library(rARPACK)
  
  n = nrow(kernel)
  an = ceiling(n*0.05)
  bn = floor(n*0.95)
  ours = rep(0, n)
  
  cusum = matrix(NA, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (j == 1 && i == 1){cusum[i, j] = kernel[i, j]}
      else if (j > 1 && i == 1){cusum[i, j] = cusum[i, j - 1] + kernel[i, j]}
      else if (j == 1 && i > 1){cusum[i, j] = cusum[i - 1, j] + kernel[i, j]}
      else {cusum[i, j] = cusum[i, j - 1] + cusum[i - 1, j]  - cusum[i - 1, j - 1] + kernel[i, j]}
    }
  }
  
  for (k in 2:n){
    for ( j in 1:(k-1) ){
      
      within = 0.5 * cusum[j, j] / j^2 + 0.5 * (cusum[k, k] - cusum[j, k] - cusum[k, j] + cusum[j, j]) / (k - j)^2
      between = (cusum[j, k] - cusum[j, j]) / (j * (k - j))
      
      tmp = 2 * (k-j)^2 * (j)^2 / (k^2 * n^2) * (within - between)
      ours[k] = max( ours[k], tmp )
      
    }
  }
  
  H = diag(1, nrow=n, ncol=n) - matrix(1, nrow=n, ncol=1) %*% matrix(1, nrow=1, ncol=n) / n
  A = H %*% kernel %*% H / n
  variance = eigs_sym(A, 1, which='LA')$values
  
  tmp = (1:n)/n * variance/2 * max(log(n)/kappa, log(50))
  n0_hat = sum(n*ours <= tmp)
  ours_threshold = mean(n*ours <= tmp) 
  ours_argmax = which.max(tmp-n*ours) / n
  
  
  return( list( ours_threshold=ours_threshold, ours_argmax=ours_argmax, ours_stat=ours*n ) )
}



matteson_general = function(d){
  
  # for any general kernel matrix, output results from matteson (tuning power of Euclidean dist)
  
  n = nrow(d)
  an = ceiling(n*0.05)
  bn = floor(n*0.95)
  matteson = rep(0, n)
  
  cusum = matrix(NA, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (j == 1 && i == 1){cusum[i, j] = d[i, j]}
      else if (j > 1 && i == 1){cusum[i, j] = cusum[i, j - 1] + d[i, j]}
      else if (j == 1 && i > 1){cusum[i, j] = cusum[i - 1, j] + d[i, j]}
      else {cusum[i, j] = cusum[i, j - 1] + cusum[i - 1, j]  - cusum[i - 1, j - 1] + d[i, j]}
    }
  }
  
  for (i in an:bn){
    within = 0.5 * cusum[i, i] / i^2 + 0.5 * (cusum[n, n] - cusum[i, n] - cusum[n, i] + cusum[i, i]) / (n - i)^2
    between = (cusum[i, n] - cusum[i, i]) / (i * (n - i))
    
    matteson[i] = i * (n - i) / n * (between - within)
  }
  
  return( list( matteson = which.max(matteson) / n ) )
}



KCPA_general = function(kernel, gamma = 1e-5){
  
  # for any general kernel matrix, output results from KCPA.
  
  
  n = nrow(kernel)
  an = ceiling(n*0.05)
  bn = floor(n*0.95)
  res = rep(NA, n)
  
  H = diag(1, nrow=n, ncol=n) - matrix(1, nrow=n, ncol=1) %*% matrix(1, nrow=1, ncol=n) / n
  A = H %*% kernel %*% H / n
  decomp = eigen(A, symmetric = TRUE)
  
  mapped = t(decomp$vectors) %*% kernel # the i-th column denote embedding for the i-th data point
  
  for (i in an:bn){
    first = apply(mapped[, 1:i], 1, mean)
    second = apply(mapped[, (i+1):n], 1, mean)
    cov_first = cov(t(mapped[, 1:i]))
    cov_second = cov(t(mapped[, (i+1):n]))
    tmp = i / n * cov_first + (n - i) / n * cov_second
    tmp_inv = solve(tmp + gamma * diag(1, n, n))
    kfdr = i * (n - i) / n * (second - first) %*% tmp_inv %*% (second - first)
    tr1 = tr(tmp_inv %*% tmp)
    tr2 = tr(tmp_inv %*% tmp_inv %*% tmp %*% tmp)
    res[i] = (kfdr - tr1) / (sqrt(2) * tr2)
  }
  
  return(list(kfdr = which.max(res) / n))
  
}
