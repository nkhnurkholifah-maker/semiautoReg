#Identification for Partially Linear Regression Model with Autoregressive errors
#----------------------------------------------------------------------
#Internal : SCAD Penalty Value and Derivate
#----------------------------------------------------------------------
scad_deriv=function(t,lambda,a=3.7){
  t=abs(t)
  out=numeric(length(t))
  out[t<=lambda]=lambda
  idx=(t>lambda)&(t<=a*lambda)
  out[idx]=(a*lambda-t[idx])/(a-1)
  out
}
scad_lqa_w=function(t,lambda,a=3.7,eps=1e-8){
  scad_deriv(pmax(abs(t),eps),lambda,a)/pmax(abs(t),eps)
}

#----------------------------------------------------------------------
#Internal: centred B-spline design matrix
#Return n x (df-1) matrix with column means=0
#----------------------------------------------------------------------
make_Zj=function(x,df=6,degree=3){
  B=splines::bs(x,df=df,degree=degree,intercept=FALSE)
  col_means=colMeans(B)
  B_c=sweep(B,2,col_means,"-")
  attr(B_c,"col_means")=col_means
  B_c
}

#----------------------------------------------------------------------
#Internal: Group SCAD threshold (block coordinate update for group j)
#Implements closed form group SCAD thresholding
#----------------------------------------------------------------------

gscad_threshold=function(u_j,norm_u,lambda,a=3.7){
  if(norm_u==0||norm_u<=lambda){
    return(rep(0,length(u_j)))
  } else if (norm_u<=a*lambda){
    factor=(a*norm_u-a*lambda)/((a-1)*norm_u)
    return (factor*u_j)
  }else{
    return(u_j) #oracle:not penalised
  }
}

#----------------------------------------------------------------------
#Internal: Build projection - purged Y* and Z*
#----------------------------------------------------------------------
purge_linear=function(Y,X,Z){
  QR=qr(X)
  Ystar=Y-qr.fitted(QR,Y)
  Zstar=Z-qr.fitted(QR,Z)
  list(Ystar=Ystar,Zstar=Zstar,QR=QR)
}

#----------------------------------------------------------------------
#Internal: AR design matrices from residual vector
#----------------------------------------------------------------------
ar_matrices=function(e,q){
  n=length(e)
  idx=(q+1):n
  V=e[idx]
  W=matrix(NA_real_,nrow=length(idx),ncol=q)
  for (i in seq_len(q)) W[,i]=e[idx-i]
  colnames(W)=paste0("lag",seq_len(q))
  list(V=V,W=W)
}

#--------------------------------------------------------------------------
#Internal: One full BCD+ SCAD-LQA iteration
#Step1: Update theta via Group SCAD Block Coordinate Descent
#Step2: Update gamma via Group SCAD LQA (scalar SCAD on each AR coefficient)
#--------------------------------------------------------------------------
bcd_step=function(Ystar,Zstar,theta,gamma,lambda,mu,
                  n_groups,Kbar,n,q,a=3.7){
  idx=(q+1):n
  n_eff=length(idx)
  Ytilde=Ystar[idx]
  Ztilde=Zstar[idx,,drop=FALSE]
  for(i in seq_len(q)){
    Ytilde=Ytilde-gamma[i]*Ystar[idx-i]
    Ztilde=Ztilde-gamma[i]*Zstar[idx-i, , drop=FALSE]
  }

  #Step1: Block Coordinate Descent over Groups
  theta_new=theta
  for(j in seq_len(n_groups)){
    cols_j=((j-1)*Kbar+1):(j*Kbar)
    cols_mj=setdiff(seq_len(ncol(Ztilde)),cols_j)
    #Partial residual : remove contribution of all other groups
    r_j=Ytilde-Ztilde[,cols_mj,drop=FALSE]%*%theta_new[cols_mj]
    #Unconstrained block OLS estimate
    Zj=Ztilde[,cols_j,drop=FALSE]
    ZjtZj=crossprod(Zj)/n_eff
    ZjtRj=crossprod(Zj,r_j)/n_eff
    u_j=tryCatch(as.vector(solve(ZjtZj,ZjtRj)),
                 error=function(e)as.vector(MASS::ginv(ZjtZj)%*%ZjtRj)
    )
    norm_u=sqrt(sum(u_j^2))
    #Group SCAD thresholding
    theta_new[cols_j]=gscad_threshold(u_j,norm_u,lambda,a)
  }
  #Step2: Update gamma with scalar SCAD-LQA
  e_new=Ystar-Zstar%*%theta_new
  am=ar_matrices(e_new,q)
  V_ar=am$V
  W_ar=am$W
  n_ar=nrow(W_ar)
  WtW_s=crossprod(W_ar)/n_ar
  WtV_s=crossprod(W_ar,V_ar)/n_ar
  lqa_g=scad_lqa_w(gamma,lambda=mu,a=a)
  A_g=WtW_s+diag(lqa_g,nrow=q)
  b_g=WtV_s
  gamma_new=tryCatch(as.vector(solve(A_g,b_g)),
                     error=function(e)as.vector(MASS::ginv(A_g)%*%b_g)
                     )
  list(theta=theta_new,gamma=gamma_new)
}

#--------------------------------------------------------------------------
#Internal: Fit for a FIXED (lambda, mu,q)
#returns list with theta, gamma, beta, convergence info
#--------------------------------------------------------------------------
fit_fixed=function(Y,X,Z_list,q,lambda,mu,max_iter=100,tol=1e-6,a=3.7){
  p=length(Z_list)
  Z=do.call(cbind,Z_list)
  n=length(Y)

  #Purge linear part
  pg=purge_linear(Y,X,Z)
  Ystar=pg$Ystar
  Zstar=pg$Zstar
  #---Initial theta:ridge (small ridge) to get stable starting point
  #using ridge instead of raw OLS gives better-scaled initial norms,
  #which helps SCAD-LQA identify the correct groups from the first iteration
  ZtZ0=crossprod(Zstar)
  ridge_pen=max(diag(ZtZ0))*0.01 #1% of max eigenvalue
  theta_init=tryCatch(as.vector(solve(ZtZ0+diag(ridge_pen,ncol(ZtZ0)),crossprod(Zstar,Ystar))),
                      error=function(e)as.vector(MASS::ginv(ZtZ0)%*%crossprod(Zstar,Ystar))
                      )
  #initial gamma:OLS on residual
  e0=Ystar-Zstar%*%theta_init
  am0=ar_matrices(e0,q)
  gamma_init=tryCatch(as.vector(solve(crossprod(am0$W),crossprod(am0$W,am0$V))),
                      error=function(e)rep(0,q)
                      )
  theta=theta_init
  gamma=gamma_init

  #Iterative LQA
  Kbar=ncol(Z_list[[1]])
  converged=FALSE
  iter_count=0L
  for(iter in seq_len(max_iter)){
    res=bcd_step(Ystar,Zstar,theta,gamma,
                 lambda,mu,n_groups=length(Z_list),
                 Kbar=Kbar,n,q,a)
    diff=max(abs(res$theta-theta),abs(res$gamma-gamma))
    theta=res$theta
    gamma=res$gamma
    iter_count=iter
    if(diff<tol) {converged=TRUE;break}
  }
  #Hard-threshold near zero to exactly zero
  theta[abs(theta)<1e-8]=0
  gamma[abs(gamma)<1e-6]=0
  #Recover beta
  beta=as.vector(solve(crossprod(X),crossprod(X,Y-Z%*%theta)))

  #Fitted values and residuals (full model)
  fitted=as.vector(X%*%beta+Z%*%theta)
  resid=Y-fitted
  list(theta=theta,gamma=gamma,beta=beta,
         fitted=fitted,residuals=resid,
         converged=converged,iterations=iter_count,
         Y_star=Ystar,Zstar=Zstar,Z=Z)
}

#--------------------------------------------------------------------------
#Internal: CV Loss for a given (lambda,mu,q) - 5 fold by default
#note:this only called when tuning_method="CV" is explicitly chosen
#Default tuning uses BIC(faster). CV is kept for users who prefer
#prediction - based tuning or when BIC/AIC result unstable
#--------------------------------------------------------------------------
cv_loss=function(Y,X,Z_list,q,lambda,mu,nfolds=5,max_iter=100,tol=1e-6){
  n=length(Y)
  folds=sample(rep(seq_len(nfolds),length.out=n))
  losses=numeric(nfolds)
  for(k in seq_len(nfolds)){
    test=which(folds==k)
    train=which(folds!=k)
    Y_tr=Y[train];X_tr=X[train, , drop=FALSE]
    Z_tr=lapply(Z_list,function(Zj) Zj[train, , drop=FALSE])
    Y_te=Y[test];X_te=X[test, , drop=FALSE]
    Z_te=do.call(cbind,lapply(Z_list,function(Zj) Zj[test, , drop=FALSE]))
    fit_k=tryCatch(fit_fixed(Y_tr,X_tr,Z_tr,q,lambda,mu,max_iter,tol),
                   error=function(e)NULL)
    if(is.null(fit_k)) { losses[k]=Inf;next}
    pred_te=as.vector(X_te%*%fit_k$beta+Z_te%*%fit_k$theta)
    losses[k]=mean((Y_te-pred_te)^2)
  }
  ok=is.finite(losses)
  if(!any(ok)) return(Inf)
  mean(losses[ok])
}

#--------------------------------------------------------------------------
#Internal: BIC/AIC Criterion
#--------------------------------------------------------------------------
ic_loss=function(fit,n,type="BIC"){
  rss=sum(fit$residuals^2)
  df=sum(fit$theta!=0)+sum(fit$gamma!=0)+length(fit$beta)
  sigma2=rss/n
  loglik=-n/2*log(2*pi*sigma2)-rss/(2*sigma2)
  if(type=="BIC")return(-2*loglik+log(n)*df)
  if(type=="AIC")return(-2*loglik+2*df)
  stop("type must be 'BIC' or 'AIC'.")
}

#--------------------------------------------------------------------------
#Internal: Auto-Select q via BIC or AR fit of initial residuals
#--------------------------------------------------------------------------
select_q=function(e,q_max){
  bic_vals=numeric(q_max)
  n=length(e)
  for(q in seq_len(q_max)){
    am=ar_matrices(e,q)
    if(nrow(am$W)<q+2){bic_vals[q]=Inf;next}
    fit=tryCatch(lm.fit(am$W,am$V),error=function(e)NULL)
    if(is.null(fit)){bic_vals[q]=Inf;next}
    rss=sum(fit$residuals^2)
    sigma2=rss/nrow(am$W)
    ll=-nrow(am$W)/2*log(2*pi*sigma2)-rss/(2*sigma2)
    bic_vals[q]=-2*ll+log(n)*q
  }
  which.min(bic_vals)
}

#--------------------------------------------------------------------------
#Main exported function
#--------------------------------------------------------------------------
#' Fit Semiparametric Partially Linear Model with Autoregressive Errors
#'
#' @description
#' Implements the penalised identification procedure of Kazemi et al. (2021)
#' for partially linear regression models with AR(q) errors. Simultaneously
#' identifies linear and nonlinear predictors and selects the AR order using
#' SCAD penalisation via Local Quadratic Approximation (LQA).
#'
#' @param prepared An object returned by  \code{semiauto_prepare()}. Required.
#' @param q Integer.AR order. If \code{NULL} (default), selected automatically
#'  via BIC on initial residuals up to \code{q_max}.
#' @param q_max Integer. Maximum AR order considered during auto-selection.
#'  Default is \code{5}
#' @param lambda Numeric. SCAD penalty for spline(nonlinear) coefficients.
#'  If \code{NULL} (default), selected via \code{tuning_method}.
#' @param mu Numeric. SCAD penalty for AR coefficients. If \code{NULL}
#'  (default), selected via \code{tuning_method}.
#' @param lambda_grid Numeric vector. Grid for \code{lambda} search. Default
#' is \code{exp(seq(log(0.05),log(3),length.out=20))}.
#' @param mu_grid Numeric vector. Grid for \code{mu} search. Default is \code{exp(seq(log(0.1),log(1),length.out=15))}.
#' @param tuning_method Character. Method for selecting tuning parameters
#'  \code{lambda} and \code{mu}.One of:
#'  \itemize{
#'    \item \code{"BIC"}-Bayesian Information Criterion(default, fast, recommended for most cases)
#'    \item \code{"AIC"}-Akaike Information Criterion (slightly less sparse than BIC, also fast)
#'    \item \code{"CV"}-5-fold cross - validation (more accurate but significantly slower;only recommended when BIC/AIC results look
#'     unstable or for small grids)
#'}
#' @param nfolds Integer. Number of CV folds. Default is \code{5}.
#' @param df Integer. Degrees of freedom for each B-spline basis(per predictor).
#' Default is \code{6} (cubic spline with 3 internal knots).
#' @param degree Integer. Degree of B-spline polynomial. Default is \code{3} (cubic).
#' @param a Numeric. SCAD shape parameter. Default is \code{3.7} (Fan & Li 2001).
#' @param max_iter Integer. Maximum LQA iterations. Default is \code{100}.
#' @param tol Numeric. Convergence tolerance. Default is \code{1e-6}.
#' @param verbose Logical. Print progress. Default is \code{TRUE}.
#'
#' @return An object of class \code{"semiauto_fit"} containing:
#' \describe{
#'  \item{beta}{Estimated linear coefficients (named vector).}
#'  \item{theta}{Estimated spline coefficients (named vector).}
#'  \item{gamma}{Estimated AR coefficients (named vector, zeros omitted).}
#'  \item{q_selected}{AR order used.}
#'  \item{lambda}{Tuning parameter used for spline penalty.}
#'  \item{mu}{Tuning parameter used for AR penalty.}
#'  \item{linear_vars}{Character vector of predictors identified as linear.}
#'  \item{nonlinear_vars}{Character vector of predictors identified as nonlinear.}
#'  \item{active_lags}{Integer vector of significant AR lags.}
#'  \item{fitted_values}{Fitted values (length n).}
#'  \item{residuals}{Residuals (length n).}
#'  \item{converged}{Logical. Did LQA converge?}
#'  \item{iterations}{Number of LQA iterations used.}
#'  \item{tuning_method}{Tuning method used.}
#'  \item{call}{The matched call.}
#' }
#' @references
#' Kazemi, M., Shahsvani, D., Arashi, M., & Rodrigues,P.C. (2021).
#' Identification for Partially Linear Model with Autoregressive Errors
#'
#'
#' @importFrom splines bs
#' @importFrom MASS ginv
#' @importFrom stats lm.fit complete.cases
#'
#' @examples
#' \dontrun{
#' prep <- semiauto_prepare(
#' data = my_data,
#' y    = "response",
#' x    = c("x1","x2","x3","x4"),
#' time = "quarter"
#' )
#'
#' fit <- semiauto_fit(prep)
#' print(fit)
#' summary(fit)
#' }
#'
#' @export

semiauto_fit=function(prepared,
                      q=NULL,
                      q_max=5L,
                      lambda=NULL,
                      mu=NULL,
                      lambda_grid=NULL,
                      mu_grid=NULL,
                      tuning_method=c("BIC","AIC","CV"),
                      nfolds=5L,
                      df=6L,
                      degree=3L,
                      a=3.7,
                      max_iter=100L,
                      tol=1e-6,
                      verbose=TRUE
){
  #Input checks
  if(!inherits(prepared,"semiauto_prepare"))
    stop("'prepared' must be an object returned by semiauto_prepare().")
  tuning_method=match.arg(tuning_method)
  model_data=prepared$model_data
  y_name=prepared$response
  x_names=prepared$predictors
  n=nrow(model_data)
  p=length(x_names)
  if(p<1) stop("At least one predictor is required.")
  if(n<2*p*df)
    warning("Sample size may be too small relative to number of basis functions.")
  #Build Matrices
  Y=as.numeric(model_data[[y_name]])
  # X:design matrix for linear part(with intercept)
  X=stats::model.matrix(stats::as.formula(paste("~",paste(x_names,collapse="+"))),
                        data=model_data)
  #Z_list:one centred B-spline matrix per predictor
  Z_list=vector("list",p)
  names(Z_list)=x_names
  for(j in seq_len(p)){
    xj=as.numeric(model_data[[x_names[j]]])
    #normalize to [0,1] for stable spline evaluation
    xj_scaled=(xj-min(xj))/(max(xj)-min(xj)+.Machine$double.eps)
    Z_list[[j]]=make_Zj(xj_scaled,df=df,degree=degree)
  }
  #Auto select q if not supplied
  if(is.null(q)){
    if(verbose) cat("Selecting AR order (q) Via BIC on initial residual...\n")
  #Quick OLS for initial residuals
    Z_all=do.call(cbind,Z_list)
    W_all=cbind(X,Z_all)
    coef0=tryCatch(lm.fit(W_all,Y)$coefficients,error=function(e)rep(0,ncol(W_all)))
    coef0[is.na(coef0)]=0
    e0=Y-W_all%*%coef0
    q=select_q(e0,q_max=min(q_max,floor((n-1)/3)))
    if(verbose)cat("AR order selected:q=",q,"\n")
  }
  q=as.integer(q)
  if(q<1) q=1L
  #Default tuning grids
  if(is.null(lambda_grid))
    lambda_grid=exp(seq(log(0.05),log(3),length.out=20))
  if(is.null(mu_grid))
    mu_grid=exp(seq(log(0.1),log(1),length.out=15))
  #Tuning Parameter Selection
  if(is.null(lambda)||is.null(mu)){
    if(tuning_method=="CV"){
      message(
        "Note:Tuning_Method='CV'selected.\n",
        "5-fold CV runs",length(lambda_grid)*length(mu_grid)*5L,
        "model fit(",length(lambda_grid),"x",length(mu_grid),
        "grid x 5 folds).\n",
        "This may be slow. Consider tuning_method = 'BIC' for faster results."
      )
    }
    if(verbose)cat("Selecting tuning parameters (Lambda,mu) via", tuning_method,"\n")
    best_score=Inf
    best_lambda=lambda_grid[1]
    best_mu=mu_grid[1]
    total_combos=length(lambda_grid)*length(mu_grid)
    combo_count=0L
    for(lam in lambda_grid){
      for(m_val in mu_grid){
        combo_count=combo_count+1L
        score=tryCatch({
          if(tuning_method=="CV"){
            cv_loss(Y,X,Z_list,q,lam,m_val,nfolds,max_iter,tol)
          }else{
            fit_tmp=fit_fixed(Y,X,Z_list,q,lam,m_val,max_iter,tol)
            ic_loss(fit_tmp,n,type=tuning_method)
          }
        },error=function(e)Inf)
        if(is.finite(score)&&score<best_score){
          best_score=score
          best_lambda=lam
          best_mu=m_val
        }
      }
      if(verbose){
        pct=round(100*combo_count/total_combos)
        cat(sprintf("Tuning:%d%% complete (best score so far:%.4f)\r",pct,best_score))
      }
    }
    if(verbose)cat("\n lambda=",round(best_lambda,5),
                   "mu=",round(best_mu,5),"\n")
    if(is.null(lambda))lambda=best_lambda
    if(is.null(mu)) mu=best_mu
  }
  #Final fit with selected(Lambda, mu, q)
  if(verbose) cat("Fitting final model----")
  fit=fit_fixed(Y,X,Z_list,q,lambda,mu,max_iter,tol,a)
  #identify linear vs nonlinear predictors
  #recompute Kbar from actual Z dimensions
  Kbar=ncol(Z_list[[1]])
  linear_vars=character(0)
  nonlinear_vars=character(0)
  theta_list=vector("list",p)
  for(j in seq_len(p)){
    cols_j=((j-1)*Kbar+1):(j*Kbar)
    theta_j=fit$theta[cols_j]
    theta_list[[j]]=theta_j
    norm_j=sqrt(sum(theta_j^2))
    if(norm_j<1e-8){
      linear_vars=c(linear_vars,x_names[j])
    }else{
      nonlinear_vars=c(nonlinear_vars,x_names[j])
    }
  }
  names(theta_list)=x_names
  #Identify Significant AR lags
  active_lags=which(abs(fit$gamma)>1e-6)
  gamma_named=fit$gamma
  names(gamma_named)=paste0("AR(",seq_len(q),")")
  #Name beta
  beta_named=fit$beta
  names(beta_named)=colnames(X)
  #Build Result
  result=list(
    call=match.call(),
    y=y_name,
    x=x_names,
    beta=beta_named,
    theta=fit$theta,
    theta_list=theta_list,
    gamma=gamma_named,
    q_selected=q,
    lambda=lambda,
    mu=mu,
    linear_vars=linear_vars,
    nonlinear_vars=nonlinear_vars,
    active_lags=active_lags,
    fitted_values=fit$fitted,
    residuals=fit$residuals,
    converged=fit$converged,
    iterations=fit$iterations,
    tuning_method=tuning_method,
    nfolds=if(tuning_method=="CV") nfolds else NA_integer_,
    df=df,
    degree=degree,
    a=a,
    n=n,
    p=p,
    Z_list=Z_list,
    X=X,
    model_data=model_data
  )
  class(result)="semiauto_fit"
  if(verbose){
    cat("Done.\n")
    print(result)
  }
  invisible(result)
}

#--------------------------------------------------------------------------
#S3 Methods
#--------------------------------------------------------------------------
#' Print method for semiauto_fit
#' @param x A \code{semiauto_fit} object.
#' @param ... Ignored.
#' @export

print.semiauto_fit=function(x, ...){
  cat("\n")
  cat("semiautoReg based on Kazemi et al. (2021)\n")
  cat(strrep("-",45),"\n")
  cat("Response         :",x$y,"\n")
  cat("n observations   :", x$n,"\n")
  cat("Predictors       :",paste(x$x,collapse=", "),"\n")
  cat("AR order (q)     :",x$q_selected,"\n")
  cat("Tuning method    :",x$tuning_method,"\n")
  cat("lambda           :",round(x$lambda,6),"\n")
  cat("mu               :",round(x$mu,6),"\n")
  cat("Iteration        :",x$iterations,
      if(x$converged)"(converged)"else"(NOT converged-increase max_iter)","\n\n")
  cat("Structure identification:\n")
  if(length(x$linear_vars)>0){
    cat("Linear:",paste(x$linear_vars,collapse = ","),"\n")
  }else {
    cat("Linear: (none)\n")
  }
  if(length(x$nonlinear_vars)>0){
    cat("Nonlinear:",paste(x$nonlinear_vars,collapse = ","),"\n")
  }else {
    cat("Nonlinear:(none)\n\n")
  }
  cat("AR structure:\n")
  if(length(x$active_lags)>0){
    cat("Active Lags:",paste0("AR(",x$active_lags,")",collapse = ","),"\n")
    cat("gamma  :",paste0(round(x$gamma[x$active_lags],4),collapse = ","),"\n")
  }else{
    cat("No significant AR lags detected.\n\n")
  }
  cat("Linear coefficients (beta):\n")
  print(round(x$beta,4))
  invisible(x)
}
#' Summary Method for semiauto_fit
#' @param object A \code{semiauto_fit} object.
#' @param ... Ignored
#' @export

summary.semiauto_fit=function(object,...){
  x=object
  cat("\n")
  cat("====================================================\n\n")
  cat("SemiautoReg Model Summary\n")
  cat("====================================================\n\n")
  cat("Call:\n");print(x$call)
  cat("\nData:\n")
  cat(" Response      :",x$y,"\n")
  cat(" n             :",x$n,"\n")
  cat(" Predictors    :",paste(x$x,collapse=","),"\n")
  cat("\nModel Fit:\n")
  rss=sum(x$residuals^2)
  tss=sum((x$model_data[[x$y]]-mean(x$model_data[[x$y]]))^2)
  r2=1-rss/tss
  rmse=sqrt(mean(x$residuals^2))
  cat("R-squared      :",round(r2,4),"\n")
  cat("RMSE           :",round(rmse,4),"\n")
  cat("Residual SS    :",round(rss,4),"\n")
  cat("\n Structure identification:\n")
  cat("Linear         :",if(length(x$linear_vars)>0)paste(x$linear_vars,collapse = ",")else
    "(none)","\n")
  cat("Nonlinear      :",if(length(x$nonlinear_vars)>0)paste(x$nonlinear_vars,collapse = ",")else
    "(none)","\n")
  cat("\nLinear coefficients (beta):\n")
  print(round(x$beta,6))
  cat("\n AR coefficients (gamma):\n")
  if(length(x$active_lags)>0){
    print(round(x$gamma[x$active_lags],6))
  }else{
    cat(" (all shrunk to zero) \n")
  }
  cat("\nTuning:\n")
  cat(" Method      :",x$tuning_method,"\n")
  cat(" lambda      :",round(x$lambda,6),"\n")
  cat(" mu          :",round(x$mu,6),"\n")
  cat(" AR order q  :",x$q_selected,"\n")
  cat(" Spline df   :",x$df,"degree:",x$degree,"\n")
  cat(" Converged   :",x$converged,paste0("(",x$iterations,"iterations)"),"\n")
  invisible(x)
}

#' Residuals method for semiauto_fit
#' @param object A \code{semiauto_fit} object.
#' @param ... Ignored.
#' @export

residuals.semiauto_fit=function(object, ...)object$residuals

#' Fitted values method for semiauto_fit
#' @param object A \code{semiauto_fit} object.
#' @param ... Ignored.
#' @export
fitted.semiauto_fit=function(object, ...) object$fitted_values

#' Predict method for semiauto_fit
#'
#' @param object A \code{semiauto_fit} object.
#' @param newdata Optional data frame with new observations. If \code{NULL},
#'   return fitted values on training data.
#' @param ... Ignored.
#' @export
predict.semiauto_fit=function(object,newdata=NULL,...){
  if(is.null(newdata)) return (object$fitted_values)
  x_names=object$x
  missing_vars=setdiff(x_names,names(newdata))
  if(length(missing_vars)>0)
    stop ("newdata is missing variables:",paste(missing_vars,collapse=","))
  #Build X_new
  X_new=stats::model.matrix(
    stats::as.formula(paste("~",paste(x_names,collapse="+"))),
    data=newdata
  )
  #Build Z_new using same basis as training
  Z_new_list=vector("list",length(x_names))
  for(j in seq_along(x_names)){
    xj=as.numeric(newdata[[x_names[j]]])
    #Re-scale using training min/max
    xj_tr=as.numeric(object$model_data[[x_names[j]]])
    lo=min(xj_tr);hi=max(xj_tr)
    xj_sc=(xj-lo)/(hi-lo+.Machine$double.eps)
    xj_sc=pmin(pmax(xj_sc,0),1)
    B_new=splines::bs(xj_sc,df=object$df,
                     degree=object$degree,intercept=FALSE)
    #Apply same centering as training
    col_means=attr(object$Z_list[[j]],"col_means")
    if(is.null(col_means))
      col_means=colMeans(object$Z_list[[j]])
    Z_new_list[[j]]=sweep(B_new,2,col_means,"-")
  }
  Z_new=do.call(cbind,Z_new_list)
  as.vector(X_new%*%object$beta+Z_new%*%object$theta)
}




