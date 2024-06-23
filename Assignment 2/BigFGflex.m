function [Hf,Hg] = BigFGflex(h,tau,A,Q,d1,d2,d3,d4)
%(sampling time,maximum delay,system mat.,eigenvector mat.)
l=eig(A);
l1=l(1);l2=l(2);%eigenvalues
Lh=[exp(l1*h),0;0,exp(l2*h)];
%Factorization of integral in Fu
F0s=[1/l1*exp(l1*h) 0;0 1/l2*exp(l2*h)];
F1s=[-1/l1 0;0 0];
F2s=[0 0;0 -1/l2];
iQ=inv(Q);
%Factorization of BIG F
F0=[expm(A*h) Q*F0s*iQ(:,2);zeros(1,3)];
F1=[zeros(2,2) Q*F1s*iQ(:,2);zeros(1,3)];
F2=[zeros(2,2) Q*F2s*iQ(:,2);zeros(1,3)];
%Factorization of integral in G1
G0s=[-inv(l1) 0;0 -inv(l2)];
G1s=[inv(l1) 0;0 0];
G2s=[0 0;0 inv(l2)];
%Factorization of BIG G
G0=[Q*G0s*iQ(:,2);1];
G1=[Q*G1s*iQ(:,2);0];
G2=[Q*G2s*iQ(:,2);0];
%Uncertain variables
ap1mi=exp(l1*h);ap2mi=exp(l2*h);
ap1ma=exp(l1*(h-tau));ap2ma=exp(l2*(h-tau));
%F Vertices
Hf.mi_mi=F0+ap1mi*F1+ap2mi*F2;
Hf.ma_ma=F0+ap1ma*F1+ap2ma*F2;
Hf.mi_ma=F0+ap1mi*F1+ap2ma*F2;
Hf.ma_mi=F0+ap1ma*F1+ap2mi*F2;
%Extra Vertcies
ap1c=exp(l1*(h-d1*tau));ap2c=exp(l2*(h-d2*tau));
Hf.comb=F0+ap1c*F1+ap2c*F2;
ap1c2=exp(l1*(h-d3*tau));ap2c2=exp(l2*(h-d4*tau));
Hf.comb2=F0+ap1c2*F1+ap2c2*F2;
%G Vertices
Hg.mi_mi=G0+ap1mi*G1+ap2mi*G2;
Hg.ma_ma=G0+ap1ma*G1+ap2ma*G2;
Hg.mi_ma=G0+ap1mi*G1+ap2ma*G2;
Hg.ma_mi=G0+ap1ma*G1+ap2mi*G2;
%Extra Vertices
Hg.comb=G0+ap1c*G1+ap2c*G2;
Hg.comb2=G0+ap1c2*G1+ap2c2*G2;
end

