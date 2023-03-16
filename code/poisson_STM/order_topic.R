order_topic = function(ref,est){
  n = nrow(ref)
  p = ncol(est)
  ref = apply(ref,2,function(z){z/norm(z,'2')})
  est = apply(est,2,function(z){z/norm(z,'2')})
  dist_mat = matrix(nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      dist_mat[i,j] = sqrt(mean((est[,i]-ref[,j])^2))
    }
  }
  order_output = c()
  for(i in 1:p){
    order_output[i] = which.min(dist_mat[,i])
    rm_idx = c()
    while(i > 1 & order_output[i]%in%order_output[-i]){
      rm_idx = c(rm_idx,order_output[i])
      temp = dist_mat[,i]
      temp[rm_idx] = Inf
      order_output[i] = which.min(temp)
    }
  }
  order_output
}
