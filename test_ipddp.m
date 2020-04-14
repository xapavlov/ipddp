
alg.tol=1e-7;
alg.maxiter=1000;
alg.mu=0; % 0 for automatic selection

alg.infeas=false;
%alg.infeas=false;

costs={};
times=[];
for i=1:50
    [funcs, fp, bp] = dynamics_invpend(i);
    %[funcs, fp, bp] = dynamics_car(i);
    %[funcs, fp, bp] = dynamics_concar(i);
    %[funcs, fp, bp] = dynamics_arm(i);
    
    [fp, bp, trace, time] = ipddp(fp, bp, funcs, alg);
    
    costs{i}=trace;
    times(i)=time;
end
