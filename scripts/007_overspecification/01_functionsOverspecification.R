####
####
# Functions used to evaluate estimability
# follows Pappalardo et al. 2020 Methods in Ecology & Evolution
# DOI: 10.1111/2041-210X.13445
# simple modifications to 01_functionsEstimability to match the parameters
# for analyzing consequences of fitting the A/C model to C/C data
####
####

# - Function to calculate credible interval/quantile of posterior and return matrix ----
credible.interval.fun = function(posterior.list, par="p.m", ci = c(.025, .975)){

  # get column index for parameter
  col.index.tmp = colnames(posterior.list[[1]]) == par;
  col.index = (1:dim(posterior.list[[1]])[2])[col.index.tmp];

  # make list of posteriors into matrix
  posterior.matrix = simplify2array(posterior.list)[,col.index,];

  # calculate 95% CI
  credible.interval = apply(posterior.matrix,2,quantile,ci);

  # return credible interval matrix
  return(credible.interval);
}

# - Function to calculate posterior mode ----

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

# - Posterior mode to matrix ----

posterior.fun <- function(posterior.list, par="p.m"){

  # get column index for parameter
  col.index.tmp = colnames(posterior.list[[1]]) == par;
  col.index = (1:dim(posterior.list[[1]])[2])[col.index.tmp];

  # make list of posteriors into matrix
  posterior.matrix = simplify2array(posterior.list)[,col.index,];

  # calculate 95% CI
  point.est = apply(posterior.matrix,2,posterior.mode);

  # return posterior mode matrix
  return(point.est);

}



# - Function to calculate 95% CI based on t-distribution ----
# code from Pappalardo et al. 2020 Methods in Ecology & Evolution
# DOI: 10.1111/2041-210X.13445
t.95CI = function(vec){

  # remove any NAs from vector
  x = na.omit(vec);

  # calculate length of vector
  n = length(x);

  # critical value: quantile at .975 for n-1 degrees of freedom
  t = qt(c(0.975), df = n - 1);

  # standard error as SD/sqrt(sample size)
  se <- sd(x)/sqrt(n);

  # 95% CI is quantile * standard error
  ci <- t * se;

  return(ci)

}

# - Function to calculate coverage of posterior ----
# takes output of credible.interval.fun()
# depends on package("binom") to calculate confidence intervals
coverage.fun = function(credible.interval.mat, par = "p.m", true.pars ){

  # get first 3 positions in par string
  par = substr(par, 1, 3)

  # reorient matrix
  credible.interval.mat = t(credible.interval.mat);

  # calculate coverage as binary
  coverage.binary = apply(credible.interval.mat,1,function(x) ifelse(x[1]<=true.pars[par] & x[2]>=true.pars[par],1,0));

  # calculate coverage proportion
  coverage.prop = sum(coverage.binary)/length(coverage.binary);

  # calculate confidence intervals on estimate
  coverage.ci = binom::binom.confint(x=sum(coverage.binary),n=length(coverage.binary),method="wilson");
  coverage.ci = coverage.ci[c("mean","lower","upper")]
  names(coverage.ci) = c("coverage.mean","coverage.lo","coverage.hi")

  # return coverage summary
  return(coverage.ci);

}

# - Function to calculate width of posterior ----
# takes output of credible.interval.fun()
width.fun = function(credible.interval.mat, par = "p.m" ){

  # reorient matrix
  credible.interval.mat = t(credible.interval.mat);

  # calculate width
  width.vec = credible.interval.mat[,2]-credible.interval.mat[,1];

  # calculate CI
  CI.up = t.95CI(width.vec)

  # create summary vector
  width.summary = c(width.mean = mean(width.vec), width.ci95 = CI.up)

  # return width summary
  return(width.summary);

}

# - Function to calculate bias of point estimate ----
# takes output of credible.interval.fun(ci=.5)
# ie. median as the point estimate
bias.fun = function(point.estimate.vec , par = "p.m" , true.pars){

  # get first 3 positions in par string
  par = substr(par, 1, 3)

  # calculate bias
  bias.vec = point.estimate.vec - true.pars[par];

  # calculate CI
  CI.up = t.95CI(bias.vec)

  # create summary vector
  bias.summary = c(bias.mean = mean(bias.vec), bias.ci95 = CI.up)

  # return bias summary
  return(bias.summary);

}

# - Function to calculate RMSE of point estimate ----
# takes output of bias.fun with median as point estimate
rmse.fun = function(bias.vec , par = "p.m" , true.pars){

  # square bias and sum (numerator of RMSE)
  numerator.rmse = sum(bias.vec^2);

  # get number of replicates
  denominator.rmse = length(bias.vec);

  # calculate rmse
  rmse = sqrt(numerator.rmse/denominator.rmse);
  names(rmse) = "rmse"

  # return rmse
  return(rmse);
}


# - Wrapper ----
# function to complete all evaluations and place in a data frame

wrapper.fun = function(posterior.list , par = "p.m", parm.df ){

  # calculate credible interval
  ci.out = credible.interval.fun(posterior.list , par = par, ci = c(.025, .975))

  # calculate point estimate
  point.out = posterior.fun(posterior.list = posterior.list, par = par)

  # construct vector of true parameters
  pars.current = as.numeric(parm.df[,c("p.m","p.g")])
  names(pars.current) = c("p.m","p.g")

  # calculate coverage
  coverage.out = coverage.fun(credible.interval.mat = ci.out, par = par, true.pars = pars.current)

  # calculate credible interval width
  width.out = width.fun(ci.out, par = par)

  # calculate bias
  bias.out = bias.fun(point.estimate.vec = point.out, par = par, true.pars = pars.current)

  # calculate rmse
  rmse.out = rmse.fun(bias.vec = bias.out, par = par, true.pars = pars.current)

  # construct dataframe
  df.out = data.frame(c(parm.df,coverage.out,width.out,bias.out,rmse.out))

  # return dataframe
  return(df.out)
}
