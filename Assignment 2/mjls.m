function [Acl_rk] = mjls(h,K,K2,q)
%mjls Closed loop matrix of MJLS
% Computes the closed-loop system matrix of the joint system,
% dependent on state q, h should be h1

%% Initialize
%Matrices
A=[-3.7 -7.5; 0 1];B=[0;1];
A2=(1/3)*A;
%Sys dimensions
nx=length(A);nu=size(B,2);
%Call System Matrices for packet loss or drop
[F1,F0]=cttodt_pl(A,B,h,'h',K);%sys1
[M1,M0]=cttodt_pl(A2,B,3*h,'h',K2);%sys2 (h1=3h2)
%% Assign Matrices to call
switch q
    case {1,3}
        Acl_rk=blkdiag(F0,eye(nx+nu));
    case {2,7}
        Acl_rk=blkdiag(F1,eye(nx+nu));
    case 4
        Acl_rk=blkdiag(F0,M1);
    case 5
        Acl_rk=blkdiag(F1,M0);
    case 6
        Acl_rk=blkdiag(F1,M1);
end

