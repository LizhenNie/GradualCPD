ours_general = function(kernel, kappa=4){
  
  # for any general kernel matrix, output results from the proposed thresholding estimator and max-gap estimator
  
  library(psych)
  library(rARPACK)
  
  n = nrow(kernel)
  an = ceiling(n*0.05)
  bn = floor(n*0.95)
  ours = rep(0, n)
  
  # calculate a CUSUM matrix, where the element in i-th row, j-th columns contains the
  # sum of all elements to the top-left block of the corresponding position in the kernel matrix
  cusum = matrix(NA, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      if (j == 1 && i == 1){cusum[i, j] = kernel[i, j]}
      else if (j > 1 && i == 1){cusum[i, j] = cusum[i, j - 1] + kernel[i, j]}
      else if (j == 1 && i > 1){cusum[i, j] = cusum[i - 1, j] + kernel[i, j]}
      else {cusum[i, j] = cusum[i, j - 1] + cusum[i - 1, j]  - cusum[i - 1, j - 1] + kernel[i, j]}
    }
  }
  
  # calaculate the scan statistic Dhat
  for (k in 2:n){
    for ( j in 1:(k-1) ){
      
      within = 0.5 * cusum[j, j] / j^2 + 0.5 * (cusum[k, k] - cusum[j, k] - cusum[k, j] + cusum[j, j]) / (k - j)^2
      between = (cusum[j, k] - cusum[j, j]) / (j * (k - j))
      
      tmp = 2 * (k-j)^2 * (j)^2 / (k^2 * n^2) * (within - between)
      ours[k] = max( ours[k], tmp )
      
    }
  }
  
  # estimate lambda1
  H = diag(1, nrow=n, ncol=n) - matrix(1, nrow=n, ncol=1) %*% matrix(1, nrow=1, ncol=n) / n
  A = H %*% kernel %*% H / n
  variance = eigs_sym(A, 1, which='LA')$values
  
  # calculate the proposed thresholding estimator and max-gap estimator
  tmp = (1:n)/n * variance/2 * max(log(n)/kappa, log(50))
  n0_hat = sum(n*ours <= tmp)
  ours_threshold = mean(n*ours <= tmp) 
  ours_argmax = which.max(tmp-n*ours) / n
  
  return( list( ours_threshold=ours_threshold, ours_argmax=ours_argmax, ours_stat=ours*n ) )
  
}
