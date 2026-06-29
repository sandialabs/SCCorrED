function [V, R] = polarDecompositionVR(A)
%polarDecompositionVR performs right polar decomposition
    [U,S,v] = svd(A);
    R = U*v';
    V = U*S*U';
end