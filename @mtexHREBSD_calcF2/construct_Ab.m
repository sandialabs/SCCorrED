%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an updated version of construct_Ab
%
% Updates:
%
%  2024-10-07: by Thomas Bennett
%   - Assigned the normalizing factor "1E-11" to a variable so it can be
%     easily modified
%   - Corrected a mathematical error in the 'traction free' case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,b] = construct_Ab(obj, refPattern, rs, qs)
    r1 = rs(1,:)';
    r2 = rs(2,:)';
    r3 = rs(3,:)';
    q1 = qs(1,:)';
    q2 = qs(2,:)';
    zerovec = zeros(size(r1));
    A1 = [r1.*r3, r2.*r3, r3.*r3, zerovec, zerovec, zerovec, -(r1.*r1 + q1.*r1), -(r1.*r2 + q1.*r2), -(r1.*r3 + q1.*r3)];
    A2 = [zerovec, zerovec, zerovec, r1.*r3, r2.*r3, r3.*r3, -(r2.*r1 + q2.*r1), -(r2.*r2 + q2.*r2), -(r2.*r3 + q2.*r3)];
    b1 = q1.*r3;
    b2 = q2.*r3;
    C_rotate = obj.ft.Qp2s'*refPattern.g';
    C_p = rotate4thorder(refPattern.C.M*1E9, C_rotate);
    n = obj.ft.Qp2s'*[0;0;1];
    factor = 1e-11; % This factor affects the relative weight of the 'free-surface'/'traction-free' constraint
    switch obj.assumption
        case 'free-surface'
            A5=[ C_p(1,1,1,1)*n(1)+C_p(1,2,1,1)*n(2)+C_p(1,3,1,1)*n(3) C_p(1,1,1,2)*n(1)+C_p(1,2,1,2)*n(2)+C_p(1,3,1,2)*n(3)    C_p(1,1,1,3)*n(1)+C_p(1,2,1,3)*n(2)+C_p(1,3,1,3)*n(3) C_p(1,1,1,2)*n(1)+C_p(1,2,1,2)*n(2)+C_p(1,3,1,2)*n(3) C_p(1,1,2,2)*n(1)+C_p(1,2,2,2)*n(2)+C_p(1,3,2,2)*n(3)  C_p(1,1,2,3)*n(1)+C_p(1,2,2,3)*n(2)+C_p(1,3,2,3)*n(3)  C_p(1,1,1,3)*n(1)+C_p(1,2,1,3)*n(2)+C_p(1,3,1,3)*n(3)      C_p(1,1,2,3)*n(1)+C_p(1,2,2,3)*n(2)+C_p(1,3,2,3)*n(3) C_p(1,1,3,3)*n(1)+C_p(1,2,3,3)*n(2)+C_p(1,3,3,3)*n(3)]*factor;
            A6=[ C_p(2,1,1,1)*n(1)+C_p(2,2,1,1)*n(2)+C_p(2,3,1,1)*n(3) C_p(2,1,1,2)*n(1)+C_p(2,2,1,2)*n(2)+C_p(2,3,1,2)*n(3)    C_p(2,1,1,3)*n(1)+C_p(2,2,1,3)*n(2)+C_p(2,3,1,3)*n(3) C_p(2,1,1,2)*n(1)+C_p(2,2,1,2)*n(2)+C_p(2,3,1,2)*n(3) C_p(2,1,2,2)*n(1)+C_p(2,2,2,2)*n(2)+C_p(2,3,2,2)*n(3)  C_p(2,1,2,3)*n(1)+C_p(2,2,2,3)*n(2)+C_p(2,3,2,3)*n(3)  C_p(2,1,1,3)*n(1)+C_p(2,2,1,3)*n(2)+C_p(2,3,1,3)*n(3)      C_p(2,1,2,3)*n(1)+C_p(2,2,2,3)*n(2)+C_p(2,3,2,3)*n(3) C_p(2,1,3,3)*n(1)+C_p(2,2,3,3)*n(2)+C_p(2,3,3,3)*n(3)]*factor;
            A7=[ C_p(3,1,1,1)*n(1)+C_p(3,2,1,1)*n(2)+C_p(3,3,1,1)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3)    C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3) C_p(3,1,1,2)*n(1)+C_p(3,2,1,2)*n(2)+C_p(3,3,1,2)*n(3) C_p(3,1,2,2)*n(1)+C_p(3,2,2,2)*n(2)+C_p(3,3,2,2)*n(3)  C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3)  C_p(3,1,1,3)*n(1)+C_p(3,2,1,3)*n(2)+C_p(3,3,1,3)*n(3)      C_p(3,1,2,3)*n(1)+C_p(3,2,2,3)*n(2)+C_p(3,3,2,3)*n(3) C_p(3,1,3,3)*n(1)+C_p(3,2,3,3)*n(2)+C_p(3,3,3,3)*n(3)]*factor;
            A = [A1; A2; A5; A6; A7];
            b = [b1; b2;  0;  0;  0];
        case 'trace=0'
            A7=[1 0 0 0 1 0 0 0 1];
            b7 = 0;
            A = [A1; A2; A7];
            b = [b1; b2; b7];
        case 'traction-free'
            % Set up the condition that the traction normal to the sample 
            % surface is zero, i.e., t_normal = n_i * sigma_ij * n_j = 0.
            % Note that sigma_ij = C_ijkl * epsilon_kl
            %                    = C_ijkl * beta_kl (by symmetry in k,l)
            A7 = [n'*C_p(:,:,1,1)*n  n'*C_p(:,:,1,2)*n  n'*C_p(:,:,1,3)*n  n'*C_p(:,:,2,1)*n  n'*C_p(:,:,2,2)*n  n'*C_p(:,:,2,3)*n  n'*C_p(:,:,3,1)*n  n'*C_p(:,:,3,2)*n  n'*C_p(:,:,3,3)*n]*factor;
            A = [A1; A2; A7];
            b = [b1; b2; 0];
    end
end
