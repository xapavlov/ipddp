function fp = forwardpass(alg, funcs, fp, bp)
  N=fp.horizon;
  dim_x=length(fp.x(:,1));
  dim_u=length(fp.u(:,1));
  dim_c=length(fp.c(:,1));
  xold=fp.x;
  uold=fp.u;
  yold=fp.y;
  sold=fp.s;
  cold=fp.c;
  tau=max(0.99, 1-alg.mu);
  steplist=2.^linspace(0,-10, 11);
  for step=1:length(steplist)
    
    xnew=zeros(dim_x,N);
    unew=zeros(dim_u,N);
    ynew=zeros(dim_c,N);
    snew=zeros(dim_c,N);
    cnew=zeros(dim_c,N);
    qnew=zeros(1,N);
    
    failed=false;
    stepsize=steplist(step);
    xnew(:,1)=xold(:,1);
    if alg.infeas
      for i=1:N
        ynew(:,i) = yold(:,i) + stepsize*bp.ky(:,i)+bp.Ky(:,:,i)*(xnew(:,i)-xold(:,i));
        snew(:,i) = sold(:,i) + stepsize*bp.ks(:,i)+bp.Ks(:,:,i)*(xnew(:,i)-xold(:,i));
        if or(any(ynew(:,i)<(1-tau)*yold(:,i)), any(snew(:,i)<(1-tau)*sold(:,i)))
          failed=1;
          break;
        end
        unew(:,i) = uold(:,i) + stepsize*bp.ku(:,i)+bp.Ku(:,:,i)*(xnew(:,i)-xold(:,i));
        xnew(:,i+1) = funcs.f([xnew(:,i); unew(:,i)]);
      end
    else
      for i=1:N 
        snew(:,i) = sold(:,i) + stepsize*bp.ks(:,i)+bp.Ks(:,:,i)*(xnew(:,i)-xold(:,i));
        unew(:,i) = uold(:,i) + stepsize*bp.ku(:,i)+bp.Ku(:,:,i)*(xnew(:,i)-xold(:,i));
        cnew(:,i) = funcs.c([xnew(:,i); unew(:,i)]);
        if or(any(cnew(:,i)>(1-tau)*cold(:,i)), any(snew(:,i)<(1-tau)*sold(:,i)))
          failed=1;
          break;
        end
        xnew(:,i+1) = funcs.f([xnew(:,i); unew(:,i)]);
      end
    end
    
    if failed
      continue;
    else
      for i=1:N
        qnew(i) = funcs.q([xnew(:,i); unew(:,i)]);
      end
      cost=sum(qnew)+funcs.p(xnew(:,N+1));
      
      if alg.infeas
        logcost=cost-alg.mu*sum(log(reshape(ynew,1,[])));   
        for i=1:N
          cnew(:,i) = funcs.c([xnew(:,i); unew(:,i)]);
        end
        err=max(alg.tol, norm(reshape(cnew+ynew,1,[]),1));
      else
        logcost=sum(qnew)+funcs.p(xnew(:,N+1))-alg.mu*sum(log(reshape(-cnew,1,[])));
        err=0;
      end
      
      candidate=[logcost;err];
      if any(all(candidate>=fp.filter,1))
        failed=2;
        continue;
      else
        idx=all(candidate<=fp.filter,1);
        fp.filter(:,idx)=[];
        fp.filter=[fp.filter, candidate];
        break;
      end      
    end
  end
  
  if failed
      fp.failed=failed;
      fp.stepsize=0;
  else
      fp.cost=cost;
      fp.logcost=logcost;
      fp.x=xnew;
      fp.u=unew;
      fp.y=ynew;
      fp.s=snew;
      fp.c=cnew;
      fp.q=qnew;
      fp.err=err;
      fp.stepsize=stepsize;
      fp.step=step;
      fp.failed=0;
  end
end