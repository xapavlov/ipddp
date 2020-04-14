function [funcs, fp, bp] = dynamics_invpend(seed)
    if nargin==0
        rng(0);
    else
        rng(seed);
    end
    N=500;
    dim_x=2;
    dim_u=1;    
    
    Q=eye(dim_x);
    R=eye(dim_u);
    
    h=0.05;
    P=10*eye(dim_x);
    
    x = sym('x_%d', [dim_x 1], 'real');
    u = sym('u','real');
    
    f = [x(1)+h*x(2); x(2)+h*sin(x(1))+h*u];
    
    fx = cat(2, diff(f,x(1)), diff(f,x(2)));
    fu = cat(2, diff(f,u));
    fxx= cat(3, diff(fx,x(1)), diff(fx,x(2)));
    fxu= cat(3, diff(fx,u));
    fuu= cat(3, diff(fu,u));

    q=0.5*h*(x'*Q*x + u'*R*u);   %stage cost
    qx=gradient(q,x);
    qu=gradient(q,u);
    qxx=cat(2, diff(qx,x(1)), diff(qx,x(2)));
    qxu=cat(2, diff(qx,u));
    quu=cat(2, diff(qu,u));
    
    p=0.5*x'*P*x;
    px=gradient(p,x);
    pxx=cat(2, diff(px,x(1)), diff(px,x(2)));
    

    c=[u-0.25; -u-0.25]%; 0.5+x(2); 0.5-x(2)]; %constraints
    %c=[u-0.1; -u-0.1; 0.5+x(2); 0.5-x(2)]; %constraints
    cx= cat(2, diff(c,x(1)), diff(c,x(2)));
    cu= cat(2, diff(c,u));
    
    dim_c=size(c,1);
    
    vars=[x;u];
    
    funcs.f=matlabFunction(f, 'vars', {vars});
    funcs.fx=matlabFunction(fx, 'vars', {vars});
    funcs.fu=matlabFunction(fu, 'vars', {vars});
    funcs.fxx=matlabFunction(fxx, 'vars', {vars});
    funcs.fxu=matlabFunction(fxu, 'vars', {vars});
    funcs.fuu=matlabFunction(fuu, 'vars', {vars});

    funcs.q=matlabFunction(q, 'vars', {vars});
    funcs.qx=matlabFunction(qx, 'vars', {vars});
    funcs.qu=matlabFunction(qu, 'vars', {vars});
    funcs.qxx=matlabFunction(qxx, 'vars', {vars});
    funcs.qxu=matlabFunction(qxu, 'vars', {vars});
    funcs.quu=matlabFunction(quu, 'vars', {vars});

    funcs.p=matlabFunction(p, 'vars', {vars});
    funcs.px=matlabFunction(px, 'vars', {vars});
    funcs.pxx=matlabFunction(pxx, 'vars', {vars});
   
    if dim_c==0
      funcs.c=[];
      funcs.cx=[];
      funcs.cu=[];
    else
      funcs.c=matlabFunction(c, 'vars', {vars});
      funcs.cx=matlabFunction(cx, 'vars', {vars});
      funcs.cu=matlabFunction(cu, 'vars', {vars});
    end
    
    funcs.dim_c=dim_c;
    funcs.dim_x=dim_x;
    funcs.dim_u=dim_u;
    
    fp.x=zeros(dim_x, N+1);
    fp.x(:,1)=[-pi;0];
    fp.u=0.02*rand(dim_u, N)-0.01;
    fp.y=0.01*ones(dim_c, N);
    fp.s=0.1*ones(dim_c, N);
    fp.mu=fp.y.*fp.s;
    fp.horizon=N;
    fp.filter=[inf;0];
    bp.ku=zeros(dim_u,N);
    bp.Ku=zeros(dim_u,dim_x,N);
    bp.ky=zeros(dim_c,N);
    bp.Ky=zeros(dim_c,dim_x,N);
    bp.ks=zeros(dim_c,N);
    bp.Ks=zeros(dim_c,dim_x,N);
end

