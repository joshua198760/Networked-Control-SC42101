function [F,G] = ctsstodtssdel(A,B,h,tau)
%ctsstodtss Converts continuous-time equations to discrete-time with delay
Fx=expm(A*h);
Fu=(expm(A*h)-expm(A*(h-tau)))*inv(A)*B;
G1=(expm(A*(h-tau))-eye(length(A)))*inv(A)*B;

nx=length(A);nu=size(B,2);
F=[Fx Fu; zeros(nu,nx+nu)];
G=[G1;eye(nu)];
end

