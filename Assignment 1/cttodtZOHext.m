function [F,G] = cttodtZOHext(A,B,h)
%Converts continuous-time equations to discrete-time with delay
nx=length(A);nu=size(B,2);

Fx=expm(A*h);
Fu=zeros(nx,nu);
G1=(expm(A*h)-eye(length(A)))*inv(A)*B;

F=[Fx Fu zeros(nx,nu); 
   zeros(nu,nx+2*nu);
   zeros(nu,nx) 1 0];
G=[G1;1;0];
end
