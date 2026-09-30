# ============================================================
# oracle_region.R
# True conditional HPD benchmark for the five simulation DGPs
#
# Integrates over the new subject's random coefficients. Gaussian cases use
# analytic HPD intervals; other cases solve a density level set numerically.
# All disconnected components are retained. Length is their total measure.
# No training data or Monte Carlo draws are needed.
# ============================================================

hpd_control <- function() list(alpha=.1,tail_tol=1e-11,quad_rel=2e-10,quad_abs=2e-12,
  root_tol=1e-10,threshold_tol=1e-12,mass_tol=1e-8,asym_mass_tol=1e-6,
  boundary_abs_tol=1e-8,boundary_rel_tol=1e-6,refine_tol=1e-7,
  mesh_scale_fraction=.5,max_mesh_points=60001L,max_refinements=3L)

oracle_parameters <- function(x,scenario) {
  if(identical(scenario,'Bimo-Fix'))scenario<-'Bimo-fix'
  scenario<-match.arg(scenario,c('Homo','Heter','Asym','Bimo','Bimo-fix'))
  x<-as.numeric(x)
  if(length(x)!=4L || any(!is.finite(x)))stop('x must have four finite covariates.')
  mu<-2*sum(x);vb<-if(scenario=='Bimo-fix')0 else sum(x^2)
  if(scenario=='Asym' && vb==0)
    stop('Asym x=0 is outside the production time design; singular pure-Gamma HPD boundary requires a separate support-boundary treatment.')
  sd<-sqrt(switch(scenario,Homo=vb+1,Heter=vb+(1+3*abs(x[2]))^2,
    Asym=vb,Bimo=vb+.25,'Bimo-fix'=.25))
  list(x=x,scenario=scenario,mu=mu,vb=vb,sd=sd)
}

oracle_distribution <- function(x,scenario,control=hpd_control()) {
  par<-oracle_parameters(x,scenario)
  state<-new.env(parent=emptyenv())
  state$pdf_evaluations<-state$cdf_evaluations<-state$sf_evaluations<-state$quadratures<-0L
  state$max_quad_abs_error<-0
  caches<-list(pdf=new.env(hash=TRUE,parent=emptyenv()),cdf=new.env(hash=TRUE,parent=emptyenv()),sf=new.env(hash=TRUE,parent=emptyenv()))
  ng_integral<-function(y,kind) {
    # G=t^10 cancels G^(-.9); exact Jacobian coefficient r^a/Gamma(a+1).
    # Beyond y+12*sd the omitted lower-CDF/normal-density kernel is negligible.
    # For SF the remaining Gamma tail is added explicitly (normal SF ~ 1).
    upper<-max(1,y+12*par$sd)
    if(kind=='sf')upper<-max(upper,qgamma(1-1e-14,.1,.1))
    cuts<-sort(unique(c(0,upper,pmax(0,pmin(upper,y+seq(-10,10,2)*par$sd)))))^.1
    fun<-function(t) {
      g<-t^10;z<-(y-g)/par$sd
      kernel<-switch(kind,pdf=dnorm(z)/par$sd,cdf=pnorm(z),sf=pnorm(z,lower.tail=FALSE))
      .1^.1/gamma(1.1)*exp(-.1*g)*kernel
    }
    ans<-0;err<-0
    for(j in seq_len(length(cuts)-1L)) {
      z<-integrate(fun,cuts[j],cuts[j+1],rel.tol=control$quad_rel,
        abs.tol=control$quad_abs/max(1,length(cuts)-1L),subdivisions=300L)
      ans<-ans+z$value;err<-err+z$abs.error;state$quadratures<-state$quadratures+1L
    }
    state$max_quad_abs_error<-max(state$max_quad_abs_error,err)
    if(kind=='sf')ans<-ans+pgamma(upper,.1,.1,lower.tail=FALSE)
    ans
  }
  calculate<-function(y,kind) {
    y<-as.numeric(y)
    if(anyNA(y))stop('Nonfinite or missing outcome argument.')
    vapply(y,function(v) {
      if(!is.finite(v))return(if(kind=='pdf')0 else if(kind=='cdf')as.numeric(v>0) else as.numeric(v<0))
      key<-sprintf('%.17g',v)
      if(exists(key,caches[[kind]],inherits=FALSE))return(get(key,caches[[kind]]))
      counter<-paste0(kind,'_evaluations');state[[counter]]<-state[[counter]]+1L
      z<-v-par$mu
      ans<-if(par$scenario=='Asym') ng_integral(z,kind) else if(par$scenario %in% c('Homo','Heter')) {
        switch(kind,pdf=dnorm(z,sd=par$sd),cdf=pnorm(z,sd=par$sd),sf=pnorm(z,sd=par$sd,lower.tail=FALSE))
      } else {
        switch(kind,pdf=.5*dnorm(z,-8,par$sd)+.5*dnorm(z,8,par$sd),
          cdf=.5*pnorm(z,-8,par$sd)+.5*pnorm(z,8,par$sd),
          sf=.5*pnorm(z,-8,par$sd,lower.tail=FALSE)+.5*pnorm(z,8,par$sd,lower.tail=FALSE))
      }
      if(!is.finite(ans) || ans< -1e-12 || (kind!='pdf'&&ans>1+1e-10))stop('Invalid conditional PDF/CDF value.')
      assign(key,ans,caches[[kind]]);ans
    },numeric(1),USE.NAMES=FALSE)
  }
  list(pdf=function(y)calculate(y,'pdf'),cdf=function(y)calculate(y,'cdf'),
       sf=function(y)calculate(y,'sf'),parameters=par,scale=par$sd,center=par$mu,
       state=state,control=control)
}

oracle_conditional_pdf <- function(y,x,scenario,control=hpd_control())oracle_distribution(x,scenario,control)$pdf(y)
oracle_conditional_cdf <- function(y,x,scenario,control=hpd_control())oracle_distribution(x,scenario,control)$cdf(y)

hpd_search_domain <- function(dist,control=dist$control) {
  width<-max(1,dist$scale);L<-dist$center-width;U<-dist$center+width
  left_expansions<-right_expansions<-0L
  while(dist$cdf(L)>control$tail_tol) {
    width<-width*2;L<-dist$center-width;left_expansions<-left_expansions+1L
    if(left_expansions>60)stop('Left tail expansion failed.')
  }
  width<-max(1,dist$scale)
  while(dist$sf(U)>control$tail_tol) {
    width<-width*2;U<-dist$center+width;right_expansions<-right_expansions+1L
    if(right_expansions>60)stop('Right tail expansion failed.')
  }
  list(bounds=c(L,U),left_tail=dist$cdf(L),right_tail=dist$sf(U),
       left_expansions=left_expansions,right_expansions=right_expansions)
}

hpd_mesh <- function(dist,domain,spacing,control=dist$control) {
  n<-max(33L,ceiling(diff(domain)/spacing)+1L)
  if(n>control$max_mesh_points)stop('Mesh exceeds prespecified numerical budget; no silent coarsening.')
  grid<-seq(domain[1],domain[2],length.out=n);den<-dist$pdf(grid)
  # Discover and refine ALL sampled extrema; no scenario-specific component count.
  mid<-2:(n-1L)
  peaks<-mid[den[mid]>=den[mid-1L]&den[mid]>=den[mid+1L]&
             (den[mid]>den[mid-1L]|den[mid]>den[mid+1L])]
  valleys<-mid[den[mid]<=den[mid-1L]&den[mid]<=den[mid+1L]&
             (den[mid]<den[mid-1L]|den[mid]<den[mid+1L])]
  extra<-numeric()
  for(j in peaks)extra<-c(extra,optimize(dist$pdf,grid[c(j-1L,j+1L)],maximum=TRUE,tol=control$root_tol)$maximum)
  for(j in valleys)extra<-c(extra,optimize(dist$pdf,grid[c(j-1L,j+1L)],maximum=FALSE,tol=control$root_tol)$minimum)
  grid<-sort(unique(c(grid,extra)));den<-dist$pdf(grid)
  list(grid=grid,density=den,peak=max(den),nodes=length(grid),spacing=max(diff(grid)),
       sampled_peaks=length(peaks),sampled_valleys=length(valleys))
}

hpd_level_set <- function(cut,dist,mesh,control=dist$control) {
  empty<-matrix(numeric(),0,2,dimnames=list(NULL,c('lo','hi')))
  grid<-mesh$grid;f<-mesh$density-cut
  if(cut<=0)return(matrix(range(grid),1,dimnames=list(NULL,c('lo','hi'))))
  if(cut>=mesh$peak)return(empty)
  crossing<-which((f[-length(f)]<0 & f[-1]>0)|(f[-length(f)]>0 & f[-1]<0))
  roots<-grid[f==0]
  for(j in crossing)roots<-c(roots,uniroot(function(y)dist$pdf(y)-cut,
    grid[c(j,j+1L)],tol=control$root_tol)$root)
  bounds<-sort(unique(c(grid[1],roots,tail(grid,1))))
  if(length(bounds)<2)return(empty)
  lo<-head(bounds,-1);hi<-tail(bounds,-1)
  keep<-dist$pdf((lo+hi)/2)>=cut
  C<-cbind(lo=lo[keep],hi=hi[keep])
  if(nrow(C)>1) {
    merged<-list(C[1,]);for(i in 2:nrow(C)) {
      j<-length(merged)
      if(C[i,1]<=merged[[j]][2]+control$root_tol)merged[[j]][2]<-C[i,2] else merged[[j+1]]<-C[i,]
    };C<-do.call(rbind,merged);colnames(C)<-c('lo','hi')
  }
  C
}

hpd_set_mass <- function(C,dist) {
  if(nrow(C)==0)return(0)
  sum(dist$cdf(C[,2])-dist$cdf(C[,1]))
}

solve_hpd_on_mesh <- function(dist,mesh,alpha,control=dist$control) {
  objective<-function(z)hpd_set_mass(hpd_level_set(z*mesh$peak,dist,mesh),dist)-(1-alpha)
  root<-uniroot(objective,c(0,1),tol=control$threshold_tol)
  cut<-root$root*mesh$peak;C<-hpd_level_set(cut,dist,mesh)
  if(nrow(C)==0 || C[1,1]<=mesh$grid[1] || tail(C[,2],1)>=tail(mesh$grid,1))
    stop('HPD set touches tail-controlled domain boundary; expand/refine required.')
  list(threshold=cut,components=C,probability=hpd_set_mass(C,dist),root_iterations=root$iter)
}

hpd_generic_solver <- function(dist,alpha=.1,control=dist$control) {
  if(length(alpha)!=1 || !is.finite(alpha) || alpha<=0 || alpha>=1)stop('Invalid alpha.')
  domain<-hpd_search_domain(dist,control);previous<-NULL;history<-list();stable<-FALSE
  for(r in 0:control$max_refinements) {
    mesh<-hpd_mesh(dist,domain$bounds,dist$scale*control$mesh_scale_fraction/2^r,control)
    ans<-solve_hpd_on_mesh(dist,mesh,alpha,control)
    discrepancy<-if(!is.null(previous)&&identical(dim(previous$components),dim(ans$components)))
      max(abs(previous$components-ans$components)) else Inf
    history[[length(history)+1]]<-data.frame(refinement=r,nodes=mesh$nodes,components=nrow(ans$components),
      probability=ans$probability,endpoint_change=discrepancy)
    if(is.finite(discrepancy)&&discrepancy<control$refine_tol){stable<-TRUE;break}
    previous<-ans
  }
  if(!stable)stop('HPD geometry did not stabilize under mesh refinement.')
  if(abs(ans$probability-(1-alpha))>control$mass_tol)stop('HPD mass tolerance failed.')
  endpoints<-as.vector(ans$components);err<-abs(dist$pdf(endpoints)-ans$threshold)
  ans$diagnostics<-list(domain=domain,mesh_refinement=do.call(rbind,history),
    boundary_abs_error=max(err),boundary_rel_error=max(err)/ans$threshold,
    pdf_at_domain_edges=dist$pdf(domain$bounds),state=as.list(dist$state))
  # Reuse only within a caller; neither training data nor a response grid enters.
  ans$mesh<-mesh
  ans
}

oracle_hpd <- function(x,scenario,alpha=.1,generic=FALSE,control=hpd_control()) {
  start<-proc.time()[['elapsed']];dist<-oracle_distribution(x,scenario,control);p<-dist$parameters
  if(p$scenario %in% c('Homo','Heter') && !generic) {
    q<-qnorm(1-alpha/2)*p$sd;C<-matrix(p$mu+c(-q,q),1,dimnames=list(NULL,c('lo','hi')))
    cut<-dist$pdf(C[1,1]);domain<-hpd_search_domain(dist)
    ans<-list(components=C,threshold=cut,probability=hpd_set_mass(C,dist),
      diagnostics=list(domain=domain,boundary_abs_error=max(abs(dist$pdf(C)-cut)),
        boundary_rel_error=max(abs(dist$pdf(C)-cut))/cut,state=as.list(dist$state)))
    algorithm<-'analytic Gaussian'
  } else {ans<-hpd_generic_solver(dist,alpha,control);algorithm<-'continuous level-set roots with mesh refinement'}
  C<-ans$components
  structure(list(scenario=p$scenario,x=p$x,alpha=alpha,threshold=ans$threshold,
    components=C,n_components=nrow(C),probability=ans$probability,
    total_length=sum(C[,2]-C[,1]),envelope=c(lo=min(C[,1]),hi=max(C[,2])),
    envelope_length=max(C[,2])-min(C[,1]),algorithm=algorithm,
    diagnostics=ans$diagnostics,runtime_sec=proc.time()[['elapsed']]-start),class='hpd_oracle')
}

evaluate_hpd_union <- function(object,y) {
  C<-object$components
  if(!is.matrix(C)||ncol(C)!=2||any(!is.finite(C))||any(C[,2]<C[,1])||
     (nrow(C)>1&&any(C[-1,1]<C[-nrow(C),2])))stop('Invalid interval union.')
  if(any(!is.finite(y)))stop('True outcomes must be finite.')
  cbind(covered=vapply(y,function(v)as.numeric(any(v>=C[,1]&v<=C[,2])),numeric(1)),
    length=rep(sum(C[,2]-C[,1]),length(y)))
}


oracle_region <- function(x_test, scenario=c('Homo','Heter','Asym','Bimo','Bimo-fix'),
                          alpha=.1, control=hpd_control()) {
  if(identical(scenario,'Bimo-Fix'))scenario<-'Bimo-fix'
  scenario<-match.arg(scenario)
  if(is.vector(x_test)&&!is.list(x_test))x_test<-matrix(x_test,nrow=1L)
  x_test<-as.matrix(x_test)
  if(ncol(x_test)!=4L || nrow(x_test)<1L || any(!is.finite(x_test)))
    stop('Production Oracle requires finite four-dimensional test covariates.')
  objects<-lapply(seq_len(nrow(x_test)),function(k)
    oracle_hpd(x_test[k,],scenario,alpha=alpha,control=control))
  # No lo_hi or training-grid region: callers must use the union evaluator.
  list(oracles=objects,components=lapply(objects,`[[`,'components'),
       density_threshold=vapply(objects,`[[`,numeric(1),'threshold'),
       conditional_probability=vapply(objects,`[[`,numeric(1),'probability'),
       total_length=vapply(objects,`[[`,numeric(1),'total_length'),
       n_components=vapply(objects,`[[`,integer(1),'n_components'),
       diagnostics=lapply(objects,`[[`,'diagnostics'),definition='true conditional HPD')
}
