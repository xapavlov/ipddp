function [funcs, fp, bp] = dynamics_car(seed)
    if nargin==0
        rng(0);
    else
        rng(seed);
    end
    N=500;
    dim_x=4;
    dim_u=2;    
    
    x = sym('x_%d', [dim_x 1], 'real');
    u = sym('u', [dim_u 1],'real');
    
    % constants
    d  = 2.0;      % d = distance between back and front axles
    h  = 0.03;     % h = timestep (seconds)

    f = x+[ (d + h*x(4)*cos(u(1)) - sqrt(d^2 - (h*x(4)*sin(u(1)))^2))*[cos(x(3));sin(x(3))];
        asin(sin(u(1))*h*x(4)/d);
        h*u(2)];
    
    fx = cat(2, diff(f,x(1)), diff(f,x(2)), diff(f,x(3)), diff(f,x(4)));
    fu = cat(2, diff(f,u(1)), diff(f,u(2)));
    fxx= cat(3, diff(fx,x(1)), diff(fx,x(2)), diff(fx,x(3)), diff(fx,x(4)));
    fxu= cat(3, diff(fx,u(1)), diff(fx,u(2)));
    fuu= cat(3, diff(fu,u(1)), diff(fu,u(2)));

    

    R=1e-2*diag([1 .01]);
    Q=1e-3*[1  1]; 
    
    q=u'*R*u + 1e-3*(sqrt(x(1)^2+Q(1))-Q(1))+1e-3*(sqrt(x(2)^2+Q(2))-Q(2));   %stage cost
    qx=gradient(q,x);
    qu=gradient(q,u);
    qxx=cat(2, diff(qx,x(1)), diff(qx,x(2)), diff(qx,x(3)), diff(qx,x(4)));
    qxu=cat(2, diff(qx,u(1)), diff(qx,u(2)));
    quu=cat(2, diff(qu,u(1)), diff(qu,u(2)));
    
    P=[.01 .01 .01  1];
    p=0.1*(sqrt(x(1)^2+P(1))-P(1))+0.1*(sqrt(x(2)^2+P(2))-P(2))+...
        1*(sqrt(x(3)^2+P(3))-P(3))+0.3*(sqrt(x(4)^2+P(4))-P(4));
    px=gradient(p,x);
    pxx=cat(2, diff(px,x(1)), diff(px,x(2)), diff(px,x(3)), diff(px,x(4)));
    
    c=[u(1)-.5; -u(1)-.5; u(2)-2; -u(2)-2;]; %constraints
    %c=[u-0.1; -u-0.1; 0.5+x(2); 0.5-x(2)]; %constraints
    cx= cat(2, diff(c,x(1)), diff(c,x(2)), diff(c,x(3)), diff(c,x(4)));
    cu= cat(2, diff(c,u(1)), diff(c,u(2)));
    
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
    fp.x(:,1)=[1;1;pi*3/2;0];

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

