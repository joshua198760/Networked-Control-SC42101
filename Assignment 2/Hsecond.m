function [Hf,Hg] = Hsecond(h,tau,A,Q)
%H Compute vertices of polytopic overapproximation
%(sampling time,maximum delay,system mat.,eigenvector mat.)
%System parameters
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
%Construct Combination
Hf.comb=F0;
Hf.comb(1,3)=Hf.ma_ma(1,3);
Hf.comb(2,3)=Hf.mi_mi(2,3);
%G Vertices
Hg.mi_mi=G0+ap1mi*G1+ap2mi*G2;
Hg.ma_ma=G0+ap1ma*G1+ap2ma*G2;
Hg.mi_ma=G0+ap1mi*G1+ap2ma*G2;
Hg.ma_mi=G0+ap1ma*G1+ap2mi*G2;
%Construct Combination
Hg.comb=G0;
Hg.comb(1,1)=Hg.ma_ma(1,1);
Hg.comb(2,1)=Hg.mi_mi(2,1);
end

