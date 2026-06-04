# =============================================================================
# semiauto_plot
utils::globalVariables(c(
  "lag", "acf", "sig",
  "x_val", "g_val",
  "index", "value", "Series",
  "fitted", "resid", "x", "y",
  "variable", "type",
  "gamma", "active"
))
# =============================================================================

sa_theme=function(){
  ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(
      plot.title = ggplot2::element_text(face="bold",size=12),
      plot.subtitle = ggplot2::element_text(colour="grey40",size=9),
      axis.title = ggplot2::element_text(size=10),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill="#f0f4f8"),
      strip.text=ggplot2::element_text(face="bold")
    )
}
SA_COL=list(
  blue="#2C7BB6", red="#D7191C",
  green="#1A9641",orange="#FDAE61",
  grey="#636363",light="#D9EDF7"
)

#P1:Nonlinear component curves
plot_nonlinear=function(x){
  if(length(x$nonlinear_vars)==0){
    return(
      ggplot2::ggplot()+
        ggplot2::annotate("text",x=0.5,y=0.5,
                          label="No nonlinear components identified",size=5,colour="grey50")+
        ggplot2::theme_void()+
        ggplot2::labs(title="Nonlinear Components")
    )
  }
  Kbar=ncol(x$Z_list[[1]])
  plot_df=do.call(rbind,lapply(x$nonlinear_vars,function(vname){
    j=which(x$x==vname)
    cols_j=((j-1)*Kbar+1):(j*Kbar)
    theta_j=x$theta[cols_j]
    xj_raw=as.numeric(x$model_data[[vname]])
    lo=min(xj_raw);hi=max(xj_raw)
    xj_seq=seq(0,1,length.out=200)
    xj_orig=lo+xj_seq*(hi-lo)
    B_seq=splines::bs(xj_seq,df=x$df,degree=x$degree,intercept=FALSE)
    col_means=colMeans(x$Z_list[[j]])
    Bc_seq=sweep(B_seq,2,col_means,"-")
    gj_seq=as.vector(Bc_seq%*%theta_j)
    data.frame(variable=vname,x_val=xj_orig,g_val=gj_seq,
               stringsAsFactors = FALSE)
  }))
  ggplot2::ggplot(plot_df,ggplot2::aes(x=x_val,y=g_val))+
    ggplot2::geom_hline(yintercept =0, linetype="dashed",
                        colour="grey60",linewidth=0.5)+
    ggplot2::geom_line(colour=SA_COL$blue,linewidth=1.1)+
    ggplot2::facet_wrap(~variable,scales="free_x")+
    ggplot2::labs(
      title="Estimated Nonlinear Component Functions",
      subtitle = "hat(g_j)(x) for variables in hat(S_2)",
      x="x",y="hat(g_j)(x)"
    )+ sa_theme()
}

#P2:Residuals vs Fitted
plot_resid_fitted=function(x){
  df=data.frame(fitted=x$fitted_values,resid=x$residuals)
  sm=as.data.frame(stats::lowess(df$fitted,df$resid))
  ggplot2::ggplot(df,ggplot2::aes(x=fitted,y=resid))+
    ggplot2::geom_hline(yintercept=0,linetype="dashed",
                        colour=SA_COL$red,linewidth=0.7)+
    ggplot2::geom_point(colour=SA_COL$grey,alpha=0.6,size=1.4)+
    ggplot2::geom_line(data=sm,ggplot2::aes(x=x,y=y),
                       colour=SA_COL$blue,linewidth=1.1)+
    ggplot2::labs(
      title="Residuals vs Fitted",
      subtitle = "Lowess smoother in blue - should be flat near zero",
      x="Fitted values",y="Residuals"
    )+sa_theme()
}

#P3: ACF of residuals
plot_acf=function(x){
  lag_max=min(30L,floor(x$n/4L))
  acf_obj=stats::acf(x$residuals,lag.max=lag_max,plot=FALSE)
  ci=stats::qnorm(0.975)/sqrt(x$n)
  acf_df=data.frame(lag=as.numeric(acf_obj$lag[-1]),
                    acf=as.numeric(acf_obj$acf[-1])
  )
  acf_df$sig=abs(acf_df$acf)>ci
  ggplot2::ggplot(acf_df,ggplot2::aes(x=lag,y=acf))+
    ggplot2::geom_hline(yintercept=0,colour="grey50")+
    ggplot2::geom_hline(yintercept=ci,linetype="dashed",
                        colour=SA_COL$blue,linewidth=0.7)+
    ggplot2::geom_hline(yintercept=-ci,linetype="dashed",
                        colour=SA_COL$blue,linewidth=0.7)+
    ggplot2::geom_segment(ggplot2::aes(xend=lag,yend=0),
                          colour=SA_COL$grey,linewidth=0.8)+
    ggplot2::geom_point(ggplot2::aes(colour=sig),size=2.2)+
    ggplot2::scale_color_manual(values=c("FALSE"=SA_COL$grey,"TRUE"=SA_COL$red),
                                guide="none")+
    ggplot2::labs(
      title="ACF of Residuals",
      subtitle=paste0("95% CI: \u00b1",round(ci,3),
                      "- red points exceed bounds"),
      x="Lag",y="Autocorrelation"
    ) + sa_theme()
}

#P4: Observed vs Fitted time series
plot_obs_fitted=function(x){
  yobs=as.numeric(x$model_data[[x$y]])
  df_long=rbind(
    data.frame(index=seq_along(yobs),value=yobs, Series="Observed"),
    data.frame(index=seq_along(yobs),value=x$fitted_values, Series="Fitted")
  )
  ggplot2::ggplot(df_long,
                 ggplot2::aes(x=index,y=value,colour=Series,linewidth=Series))+
    ggplot2::geom_line()+
    ggplot2::scale_colour_manual(
      values=c("Observed"=SA_COL$grey,"Fitted"=SA_COL$red))+
    ggplot2::scale_discrete_manual("linewidth",values = c("Observed"=0.6,"Fitted"=1.1))+
    ggplot2::labs(
      title="Observed vs Fitted Time Series",
      subtitle = paste0("Response:", x$y),
      x="Index",y=x$y, colour=NULL,linewidth=NULL
    )+ sa_theme()+
    ggplot2::theme(legend.position = "top")
}

#P5:Structure identification bar chart
plot_structure=function(x){
  p=length(x$x)
  Kbar=ncol(x$Z_list[[1]])
  var_df=do.call(rbind,lapply(seq_len(p),function(j){
    cols_j=((j-1)*Kbar+1):(j*Kbar)
    data.frame(
      variable=x$x[j],
      norm=sqrt(sum(x$theta[cols_j]^2)),
      type=ifelse(x$x[j]%in%x$linear_vars,"Linear","Nonlinear"),
      stringsAsFactors = FALSE
    )
      }))
  var_df$variable=factor(var_df$variable,
                         levels=var_df$variable[order(var_df$norm)])
  ggplot2::ggplot(var_df,
                  ggplot2::aes(x=variable,y=norm,fill=type))+
    ggplot2::geom_col(width=0.6)+
    ggplot2::geom_hline(yintercept=0,colour="grey30")+
    ggplot2::coord_flip()+
    ggplot2::scale_fill_manual(
      values=c("Linear"=SA_COL$blue,"Nonlinear"=SA_COL$orange),
      name="Identified as")+
    ggplot2::labs(
      title="Structure Identification Summary",
      subtitle="||theta_j||=0->Linear| ||theta_j||>0 ->Nonlinear",
      x=NULL,y="||hat(theta_j)||"
    )+ sa_theme()+
  ggplot2::theme(legend.position = "top")
}

#P6 : AR coefficient bar chart
plot_ar_coefs=function(x){
  q=x$q_selected
  gamma_df=data.frame(
    lag=factor(paste0("AR(",seq_len(q),")"),
               levels=paste0("AR(",seq_len(q),")")),
    gamma=as.numeric(x$gamma),
    active=seq_len(q) %in% x$active_lags
  )
  ggplot2::ggplot(gamma_df,
                  ggplot2::aes(x=lag,y=gamma,fill=active))+
    ggplot2::geom_col(width=0.55)+
    ggplot2::geom_hline(yintercept = 0,colour="grey30")+
    ggplot2::scale_fill_manual(values=c("TRUE"=SA_COL$green,"FALSE"="grey80"),
                               labels=c("TRUE"="Significant","FALSE"="Shrunk to zero"),
                               name=NULL)+
    ggplot2::labs(
      title="Estimated AR Coefficients",
      subtitle = paste0("q = ", q,"| Active lags:",
                        if(length(x$active_lags)>0)
                           paste0("AR(",x$active_lags,")",collapse=",")
                           else "none"),
                        x=NULL,y="hat(gamma_i)"
    )+sa_theme()+
    ggplot2::theme(legend.position = "top")
}

#' Plot diagnostics for a semiauto_fit object (ggplot2)
#'
#' @description
#' Produces up to 6 ggplot2 panels:
#' \enumerate{
#'  \item Estimated nonlinear component functions
#'  \item Residuals vs Fitted Values
#'  \item ACF of residuals
#'  \item Observed vs Fitted time series
#'  \item Structure identification summary
#'  \item AR coefficient estimates }
#' @param x A\code{semiauto_fit} object.
#' @param which Integer vector. Which panels to draw. Default \code{1:6}.
#' @param combine Logical.
#' If TRUE, panels are combine using patchwork.
#' IF FALSE, a named list of ggplot objects is returned.
#' Default is FALSE.
#' @param title Character. Optional overall title.
#' @param ... Ignored.
#'
#' @return A \code{patchwork} figure (if \code{combine=TRUE}) or a named list
#' of \code{ggplot} objects.
#'
#' @examples
#' \dontrun{
#' data(macro_turkey)
#'
#' prep=semiauto_prepare(
#'  data = macro_turkey,
#'  y = "inflasi",
#'  x = c("kurs", "m3", "industri", "kredit", "brent"),
#'  time = "Tarih"
#'  )
#'
#'  fit = semiauto_fit(prepared=prep)
#'
#'  #Return diagnostic plots as a list
#'  plots = semiauto_plot(fit)
#'
#'  names(plots)
#'
#'  plots$nonlinear
#'  plots$resid_fitted
#'  plots$acf
#'  plots$obs_fitted
#'  plots$structure
#'  plots$ar_coefs
#'
#'  #Combine all panels
#'  semiauto_plot(fit,combine=TRUE)
#' }
#'
#' @export

semiauto_plot=function(x,which=1:6, combine=FALSE, title=NULL, ...){
  if(!inherits(x,"semiauto_fit"))
    stop("'x' must be a semiauto_fit object.")
  if(!requireNamespace("ggplot2",quietly=TRUE))
    stop("Install ggplot2: install.packages('ggplot2')")
  if(combine&&!requireNamespace("patchwork",quietly=TRUE))
    stop("Install patchwork: install.packages ('patchwork') ")
  builders=list(
    nonlinear=plot_nonlinear,
    resid_fitted=plot_resid_fitted,
    acf= plot_acf,
    obs_fitted= plot_obs_fitted,
    structure=plot_structure,
    ar_coefs=plot_ar_coefs
  )
  panels=list()
    for(i in seq_along(builders)){
      if(i %in% which)
        panels[[names(builders)[i]]]=builders[[i]](x)
    }
    if(!combine||length(panels)==0) return(panels)
    #Assemble: P1 full width, rest 2-per-row
    n_nl=length(x$nonlinear_vars)
    nl_rows=max(1L,ceiling(n_nl/3L))
    others=panels[names(panels)!="nonlinear"]
    if("nonlinear" %in% names(panels) &&length (others)>0){
      fig=panels$nonlinear/
        patchwork::wrap_plots(others,ncol=2)+
        patchwork::plot_layout(heights = c(2.5,3))
    }else if ("nonlinear" %in% names(panels)){
      fig=panels$nonlinear
    }else{
      fig=patchwork::wrap_plots(panels,ncol=2)
    }
    ann_title=if(!is.null(title))title else
      paste0("semiautoReg Diagnostics - response:",
             x$y," [n=", x$n, ", q=", x$q_selected, "]")
    fig+patchwork::plot_annotation(
      title=ann_title,
      subtitle = paste0("Linear:",paste0(x$linear_vars, collapse=", ")," | ",
                        "Nonlinear:",paste0(x$nonlinear_vars, collapse=", ")," | ",
                        "Active AR:", if (length(x$active_lags)>0)
                          paste0("AR(",x$active_lags,")", collapse = ",")
                        else "none"),
      theme=ggplot2::theme(
        plot.title=ggplot2::element_text(face="bold",size=14),
        plot.subtitle = ggplot2::element_text(colour="grey40",size=10)
      )
    )
}
