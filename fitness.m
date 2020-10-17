function f=fitness(gene)
highomega=5;
lowomega=1;
v1=gene(1); v2=gene(2); epsilon1=gene(3); epsilon2=gene(4); alpha1=gene(5); alpha2=gene(6);
A=0.5; B1=1.5-v1; B2=1.5-v2; C=2; D=10; E=4; F=2; 
tend=300+ceil(5*rand());
tspan=[0 tend];
ICs(1)=0.99+0.01*rand();
ICs(2)=0.01*rand(); ICs(3)=0.01*rand(); ICs(4)=0.01*rand();%initial conditions
[t1,y1]=ode45(@(t,y) odefun(t,y,A,B1,B2,C,D,E,F,epsilon1,epsilon2,highomega),tspan,ICs);
[t2,y2]=ode45(@(t,y) odefun(t,y,A,B1,B2,C,D,E,F,epsilon1,epsilon2,lowomega),tspan,ICs);
out1=alpha1*y1(:,1)+alpha2*y1(:,2);
out2=alpha1*y2(:,1)+alpha2*y2(:,2);
tspannew=[tend/2:0.01:tend];
outnew1=interp1(t1,out1,tspannew);
outnew2=interp1(t2,out2,tspannew);
f1=mean(abs(outnew1-1));
f2=mean(abs(outnew2-0));
f=f1+f2;
end
function dydt=odefun(t,y,A,B1,B2,C,D,E,F,epsilon1,epsilon2,omega)
p1=y(1); p2=y(2); y1=y(3); y2=y(4);
p_square=p1.^2+p2.^2;
y_square=y1.^2+y2.^2;
p_quad=p1.^4+p2.^4;
dydt=[p1*(F*(1-p_square)+D*(p1.^2*p_square-p_quad))+E*(-y1.^2*p1*p2+y2.^2*p2.^2);
    p2*(F*(1-p_square)+D*(p2.^2*p_square-p_quad))+E*(-y2.^2*p2*p1+y1.^2*p1.^2);
    -y1*((y1.^2-1).^2+A-B1*p1.^2+C*(y_square-y1.^2))+epsilon1*sin(omega*t)*sin(omega*t);
    -y2*((y2.^2-1).^2+A-B2*p2.^2+C*(y_square-y2.^2))+epsilon2*sin(omega*t)*sin(omega*t)];
end