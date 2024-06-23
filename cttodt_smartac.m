function [F,G] = cttodt_smartac(A,B,h)
%Converts continuous-time equations to discrete-time with delay
nx=length(A);nu=size(B,2);

Fx=expm(A*h);
Gu=expm(A*h)*inv(A)*B+1/h*(eye(nx)-expm(A*h))*(A^-2)*B;
Fu=1/h*(expm(A*h)-eye(nx))*(A^-2)*B - inv(A)*B;

F=[Fx Fu Gu; 
   zeros(nu,nx+2*nu);
   zeros(nu,nx) 1 0];
G=[zeros(nx,nu);1;0];
end