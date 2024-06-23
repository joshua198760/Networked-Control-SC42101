function [F,G] = ctsstodtss(A,B,h)
%ctsstodtss Converts continuous-time equations to discrete-time
F=expm(A*h);
G=(expm(A*h)-eye(length(A)))*inv(A)*B;
end

