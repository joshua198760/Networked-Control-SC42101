function [F,G] = Poly(h,tau,A,Q)
%H Compute vertices of polytopic overapproximation
%(sampling time,maximum delay,system mat.,eigenvector mat.)

%System parameters
l=eig(A);
l1=l(1);l2=l(2);%eigenvalues
Lh=diag(l1,l2);
%Uncertain Parameters
ap1=exp(l1*(h-tau));ap2=exp(l2*(h-tau));
%Factorization of integral in Fu
F0s=[1/l1*exp(l1*h) 0;0 1/l2*exp(l2*h)];
F1s=[-1/l1 0;0 0];
F2s=[0 0;0 -1/l2];
%Factorization of BIG F
F0=[inv(Q)*Lh*Q inv(Q)*F0s*Q(:,2);zeros(1,3)];
F1=[zeros(2,2) inv(Q)*F1s*Q(:,2);zeros(1,3)];
F2=[zeros(2,2) inv(Q)*F2s*Q(:,2);zeros(1,3)];
F=F0+ap1*F1+ap2*F2;
%Factorization of integral in G1
G0s=[-inv(l1) 0;0 -inv(l2)];
G1s=[inv(l1) 0;0 0];
G2s=[0 0;0 inv(l2)];
%Factorization of BIG G
G.G0=[inv(Q)*G0s*Q(:,2);1];
G.G1=[inv(Q)*G1s*Q(:,2);0];
G.G2=[inv(Q)*G2s*Q(:,2);0];
G=G0+ap1*G1+ap2*G2;
end

