function [fp, bp, costs, time] = ipddp(fp, bp, funcs, alg)
  tic;
  costs=[];
  
  fp = initialroll(fp, funcs);
  
  if alg.mu==0
    alg.mu=fp.cost/fp.horizon/length(fp.s(:,1));
  end
  
  fp = resetfilter(fp, alg);
  bp = resetreg(bp);
  
  for iter = 1:alg.maxiter
    while 1
      bp = backwardpass(alg, funcs, fp, bp);
      if not(bp.failed) break; end
    end
    fp = forwardpass(alg, funcs, fp, bp);
    
    time=toc;
    if mod(iter,10)==1
        fprintf('\n');
        fprintf('%-12s','Iteration','Time','mu','Cost','Opt. error','Reg. power','Stepsize');
        fprintf('\n');
    end
    fprintf('%-12d%-12.4g%-12.4g%-12.4g%-12.4g%-12d%-12.3f\n', ...
                iter, time, alg.mu, fp.cost, bp.opterr, bp.reg, fp.stepsize);
    
    costs=[costs; fp.cost];
    steps=[steps; fp.stepsize];
    
    if max(bp.opterr,alg.mu)<=alg.tol
       disp('~~~Optimality reached~~~');
       break;
    end
    
    if bp.opterr<=0.2*alg.mu
       alg.mu=max(alg.tol/10, min(0.2*alg.mu, alg.mu^1.2));
       fp = resetfilter(fp, alg);
       bp = resetreg(bp);
    end
  end
end

function fp = initialroll(fp, funcs)
  N=fp.horizon;
  for i=1:N
    x=fp.x(:,i);
    u=fp.u(:,i);
    fp.c(:,i)=funcs.c([x; u]);
    fp.q(i)=funcs.q([x; u]);
    fp.x(:,i+1) = funcs.f([x; u]);
  end
  fp.cost= sum(fp.q)+funcs.p(fp.x(:,N+1));
end

function fp = resetfilter(fp, alg)
    if alg.infeas
      fp.logcost=fp.cost-alg.mu*sum(log(reshape(fp.y,1,[])));
      fp.err=norm(reshape(fp.c+fp.y,1,[]),1);
      if fp.err<alg.tol
        fp.err=0;
      end
    else
      fp.logcost=fp.cost-alg.mu*sum(log(reshape(-fp.c,1,[])));
      fp.err=0;
    end
    fp.filter=[fp.logcost; fp.err];
    fp.step=0;
    fp.failed=0;
end  

function bp = resetreg(bp)
    bp.reg=0;
    bp.failed=0;
    bp.recovery=0;
end


