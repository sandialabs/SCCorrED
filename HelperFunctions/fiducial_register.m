function [R, t] = fiducial_register(X,Y)
    dim = size(X, 1);
    if dim == 2
        X(3,:) = zeros(1, size(X,2));
        Y(3,:) = zeros(1, size(X,2));
    end
    Xc = X - repmat(mean(X,2),1,size(X,2));
    Yc = Y - repmat(mean(Y,2),1,size(X,2));
    [U,~,V] = svd(Xc*Yc');
    R = V*diag([1, 1, det(V*U)])*U';
    t = mean(Y,2) - R*mean(X,2);
end