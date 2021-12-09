model {
  
  # - Priors ---------------------------------------------------------------
  
  # - +Probability of  mortality ----
  p.m ~ dbeta(1,1)
  
  # - +Probability of  germination ----
  p.g ~ dbeta(1,1)
  
  # - Transformed parameters ---------------------------------------------------------------
  
  # - +Mortality history ----
  p.s[1] = (1-p.m)
  p.s[2] = (1-p.m)*(1-p.m)
  p.s[3] = (1-p.m)*(1-p.m)*(1-p.m)
  
  # - +Germination history ----
  theta_c[1] = 1
  theta_c[2] = (1-p.g)
  theta_c[3] = (1-p.g)^2
  
  # - Likelihoods ---------------------------------------------------------------
  
  # - +Intact seed observations ----
  for(i in 1:n_s){
    
    # - +Deterministic survival model ----
    # multiplies mortality and germination histories
    mu[i] <- theta_c[t_obs[i]]*p.s[t_obs[i]]
    
    # - +Likelihood ----
    y.s_obs[i] ~ dbinom(mu[i], n.s_obs[i])
    
  }
  
  # - +Seedling observations ----
  for(i in 1:n_g){
    
    # - +Likelihood ----
    y.g_obs[i] ~ dbinom(p.g, n.g_obs[i])
    
  }
  
}
