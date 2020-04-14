function [funcs, fp, bp] = dynamics_arm(seed)
    if nargin==0
        rng(0);
    else
        rng(seed);
    end
    N=500;
    dim_x=3;
    dim_u=3;    
    
    x = sym('x_%d', [dim_x 1], 'real');
    u = sym('u', [dim_u 1],'real');
    
    % constants
    l=1;
    h  = 0.1;     % h = timestep (seconds)

    f= x+h*u;
    
 
    fx = cat(2, diff(f,x(1)), diff(f,x(2)), diff(f,x(3)));
    fu = cat(2, diff(f,u(1)), diff(f,u(2)), diff(f,u(3)));
    fxx= cat(3, diff(fx,x(1)), diff(fx,x(2)), diff(fx,x(3)));
    fxu= cat(3, diff(fx,u(1)), diff(fx,u(2)), diff(fx,u(3)));
    fuu= cat(3, diff(fu,u(1)), diff(fu,u(2)), diff(fu,u(3)));

    rx2=l*sin(x(1)) + l*sin(x(1)+x(2));
    ry2=l*cos(x(1)) + l*cos(x(1)+x(2));
    rx3=rx2 + l*sin(x(1)+x(2)+x(3));
    ry3=ry2 + l*cos(x(1)+x(2)+x(3));

    R=eye(dim_u);
    Q=eye(dim_x);

    
    q=0.5*h*u'*R*u;% + 0.5*h*[rx;ry]'*eye(2)*[rx;ry];   %stage cost
    qx=gradient(q,x);
    qu=gradient(q,u);
    qxx=cat(2, diff(qx,x(1)), diff(qx,x(2)), diff(qx,x(3)));
    qxu=cat(2, diff(qx,u(1)), diff(qx,u(2)), diff(qx,u(3)));
    quu=cat(2, diff(qu,u(1)), diff(qu,u(2)), diff(qu,u(3)));
    
    Q=diag([100; 10; 1]);
    r=x-[pi/2; 0; 0];
    p=0.5*(r'*Q*r);
    px=gradient(p,x);
    pxx=cat(2, diff(px,x(1)), diff(px,x(2)), diff(px,x(3)));
    
    c=[u-0.10; -u-0.10; ry2-1;-ry2-1;ry3-1;-ry3-1;]; %constraints
    cx= cat(2, diff(c,x(1)), diff(c,x(2)), diff(c,x(3)));
    cu= cat(2, diff(c,u(1)), diff(c,u(2)), diff(c,u(3)));
    
    cxx= cat(3, diff(cx,x(1)), diff(cx,x(2)), diff(cx,x(3)));
    cxu= cat(3, diff(cx,u(1)), diff(cx,u(2)), diff(cx,u(3)));
    cuu= cat(3, diff(cu,u(1)), diff(cu,u(2)), diff(cu,u(3)));
    
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
      funcs.cxx=matlabFunction(cxx, 'vars', {vars});
      funcs.cxu=matlabFunction(cxu, 'vars', {vars});
      funcs.cuu=matlabFunction(cuu, 'vars', {vars});
    end
    
    funcs.dim_c=dim_c;
    funcs.dim_x=dim_x;
    funcs.dim_u=dim_u;
    
    fp.x=zeros(dim_x, N+1);
    fp.x(:,1)=[-pi/2;0;0];

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

