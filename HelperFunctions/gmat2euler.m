function [phi1, PHI, phi2]=gmat2euler(g)
%gmat2euler - retrieves the euler angles from a g-matrix
%   according to bunge for phi1,PHI,phi2 in radians
%
% Corrected by TJB, 2024-11-15
    tol = 1e-10;
    if g(3,3) > 1-tol
        PHI=0.0;
        if g(1,1) > 1-tol
            phi1=0.0;
        elseif g(1,1) < -1.0
            phi1 = pi;
        else
            if g(1,1) > 1
                temp=1.0;
            elseif g(1,1) < -1
                temp=-1.0;
            else
                temp=g(1,1);
            end
            phi1=acos(temp);
        end
        if g(1,2) < 0.0
            phi1 = 2*pi-phi1;
        end
        phi2=0.0;
    elseif g(3,3) < -1+tol
        PHI=pi;
        if g(1,1) > 1-tol
            phi1 = 0.0;
        elseif g(1,1) < -1.0
            phi1 = pi;
        else
            if g(1,1) > 1
                temp=1.0;
            elseif g(1,1) < -1
                temp=-1.0;
            else
                temp=g(1,1);
            end
            phi1=acos(temp);
        end
        if g(1,2) < 0.0
            phi1 = 2*pi-phi1;
        end
        phi2=0.0; % corrected, TJB 2024-11-15
    else
        if g(3,3) > 1
            temp=1.0;
        elseif g(3,3) < -1
            temp=-1.0;
        else
            temp=g(3,3);
        end
        PHI=acos(temp);
	    phi1 = atan2(g(3,1),-g(3,2));
	    phi2 = atan2(g(1,3), g(2,3));
        if phi1 < 0.0
            phi1 =phi1+2*pi;
        end
        if phi2 < 0.0
            phi2 =phi2+2*pi;
        end
    end
end