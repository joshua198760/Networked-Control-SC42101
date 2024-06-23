function [F1,F0] = cttodt_pl(A,B,h,s,K)
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
        zeros(nu,nx) 1];
end

end

