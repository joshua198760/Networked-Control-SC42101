function [F1,F0] = cttodtq5sys2(A,B,h,s,K)
nx=length(A);nu=size(B,2);

Fx=expm(A*h);
G=(expm(A*h)-eye(length(A)))*inv(A)*B;

if s=='z'
    F0=Fx-G*K;
    F1=Fx;
elseif s=='h'
    F0=[Fx-G*K zeros(nx,nu);
        -K zeros(nu,nu)];
    F1=[Fx G;
        zeros(nu,nx) zeros(nu,nu)];
end

end
