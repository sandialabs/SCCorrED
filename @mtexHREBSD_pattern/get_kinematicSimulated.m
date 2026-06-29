function I = get_kinematicSimulated(obj)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    g = obj.g;
    xstar = obj.patternCenter(1);
    ystar = obj.patternCenter(2);
    zstar = obj.patternCenter(3);
%     pixsize = mtexHREBSD.scan.pixelSize;
%     Av = mtexHREBSD.scan.KV*1E3;
%     sampletilt = mtexHREBSD.scan.sampleTilt*pi/180;
%     elevang = mtexHREBSD.scan.cameraElevation*pi/180;
%     material = mtexHREBSD.readMaterial(obj.material);


    pixsize = obj.SimData.pixelSize;
    Av = obj.SimData.Av;
    sampletilt = obj.SimData.sampletilt;
    elevang = obj.SimData.elevang;
    material = obj.SimData.material;

    Fhkl = material.Fhkl;
    dhkl = material.dhkl;
    hkl = material.hkl;
    axs = material.axs;
    lattice = material.lattice;
    alattice = material.a1;
    blattice = material.b1;
    clattice = material.c1;
    F = obj.SimData.F;
    useeuler = 0;

    sFhkl = Fhkl.^2;
    Wa = 6.626e-34/sqrt(2*1.602e-19*9.109e-31*Av+(1.602e-19*Av/2.998e8)^2); %calculate wavelength
    
    simpat = zeros(pixsize,pixsize);
    
    
    %Coordinate frame transformation
    Qvp=[-1 0 0; 0 -1 0; 0 0 1];
    Qps = obj.SimData.ft.Qp2s;
    
    PC = [xstar*pixsize;(1-ystar)*pixsize;zstar*pixsize];
    % Sample to Crystal
    if length(g(:)) < 9
        phi1=g(1);
        PHI=g(2);
        phi2=g(3);
        Qsc=euler2gmat(phi1,PHI,phi2); %rotation sample to crystal
    else
        Qsc=g;
        [R U] = poldec(Qsc);
        if sum(sum(U-eye(3)))>1e-10
            error('g must be a pure rotation')
        end
    end
    Qcs=Qsc';
    
    % Phospher to sample
    % I wrote frameTransforms.phosphorToSample to consolidate this
    % functionality into one place, but there isn't an easy way to call that
    % function from here. Eventually, I want to use that function here, but it
    % is a larger refactor than I want to deal with right now. --Zach Clayburn
%     if ~useeuler
%         alpha = pi/2 - sampletilt + elevang;
%         % Phospher to sample
%         Qps=[0 -cos(alpha) -sin(alpha);...
%             -1     0            0;...
%             0   sin(alpha) -cos(alpha)];
%     else
%         Qmp = euler2gmat(camphi1,camPHI,camphi2);
%         Qmi = [0 -1 0;1 0 0;0 0 1];
%         Qio = [cos(sampletilt) 0 -sin(sampletilt);0 1 0;sin(sampletilt) 0 cos(sampletilt)];
%         Qpo = Qio*Qmi*Qmp'*[-1 0 0;0 1 0;0 0 -1];
%         
%         Qps = Qpo;
%     end
    
    UsePermHKL = 0; %can't get this working for anything non-cubic, so use old 
    %method for other symmetries. PermuteHKL does speed up cubic symmetries
    %significantly and removes double counting. Some day someone should get it
    %working for other symmetries as well and then use genEBSDVect.m
    if strcmp(lattice,'hexagonal') == 1
        SymOps = gensymopsHex;
        numsyms = 12;
    elseif strcmp(lattice,'tetragonal')
        SymOps = gensymopsTet(axs);
        numsyms = 8;
    elseif strcmp(lattice,'cubic')
        UsePermHKL = 1;
    end
%     warning('off', 'MATLAB:colon:operandsNotRealScalar')
    for i = 1:length(dhkl)
        
        if UsePermHKL
            NewHKLList = PermuteHKL(hkl(i,:),lattice);
    %         for j = 1:size(NewHKLList,1)%numsyms
                numsyms = size(NewHKLList,1);
        end
        for j = 1:numsyms
            if UsePermHKL
                      eco3 = NewHKLList(j,:)';
            else
                eco3=hkl(i,:);
            end
            if eco3(1) ~= 0
                eco3(1) = 1/eco3(1);
                eco3(1) = eco3(1)*alattice;
                %             eco3(1) = eco3(1)*1/alattice;
                eco3(1) = 1/eco3(1);
            end
            if eco3(2) ~= 0
                eco3(2)= 1/eco3(2);
                eco3(2) = eco3(2)*blattice;
                %             eco3(2) = eco3(2)*1/blattice;
                eco3(2)= 1/eco3(2);
            end
            if eco3(3) ~= 0
                eco3(3) = 1/eco3(3);
                eco3(3) = eco3(3)*alattice; % ****this was clattice but is now alattice due to changes in hkl values for HCP, may need to fix for tetragonal
                %             eco3(3) = eco3(3)*1/clattice;
                eco3(3) = 1/eco3(3);
            end
            if ~UsePermHKL
                eco3 = squeeze(SymOps(j,1:3,1:3))*eco3';
            end
            %         N=eco3/norm(eco3);
         
            C=eco3/norm(eco3)*dhkl(i);% normal to hkl plane, with length equal to distance between planes
            if C(3) == 0%my stuff
                A = [0 0 1]';% normal vector to C
            else%my stuff
                A = [1 1 -(C(1)+C(2))/C(3)]';% normal vector to C
            end%my stuff
     
            B = cross2(A,C);% normal vector to A and C
            NPreNorm = cross2(A,B);% this should be in direction of C (normal to hkl plane) - can't remember why all the work with A and B
            N = NPreNorm/norm(NPreNorm);
            n = norm(NPreNorm)/norm(cross2(F*A,F*B))*det(F)*inv(F)'*N;% deformed normal to hkl plane in crystal frame
            tdhkl = abs(xdotyMex((F*C)',n));% deformed distance between planes????
            theta=real(asin(Wa/2/tdhkl));
            %         eco3=Qcs*n';
            eco3 = Qcs*n; % normal to hkl plane in sample frame; ecoi are the axes of the reference frame associated with the cone for this hkl frame
            %Find the in plane bases in the crystal frame
            %         eco3=F*eco3;
            %         eco3=Qcs*eco3;
            eco3=eco3';
            eco3=eco3/norm(eco3);
            eco2=[0 0 0];
            if eco3(1)==0
                eco2(1)=1;
            elseif eco3(2)==0
                eco2(2)=1;
            elseif eco3(3)==0
                eco2(3)=1;
            else
                eco2(1)=1;
                eco2(2)=1;
                eco2(3)=(-eco3(1)-eco3(2))/eco3(3);
                eco2=eco2/norm(eco2);
            end
            eco1=cross(eco2,eco3);
            
            Qsco=[eco1; eco2; eco3];
            
            Q = Qsco*Qps*Qvp;
            
            %Equation of intersection
            t=-Q*PC;
            
            ts=tan(theta)^2;
    %         keyboard
            a=ts*(Q(1,1)^2+Q(2,1)^2)-Q(3,1)^2;%x^2
            b=ts*(Q(1,2)^2+Q(2,2)^2)-Q(3,2)^2;%y^2
            c=ts*(2*Q(1,1)*Q(1,2)+2*Q(2,1)*Q(2,2))-2*Q(3,1)*Q(3,2);%xy
            d=ts*(2*Q(1,1)*t(1)+2*Q(2,1)*t(2))-2*Q(3,1)*t(3);%x
            e=ts*(2*Q(1,2)*t(1)+2*Q(2,2)*t(2))-2*Q(3,2)*t(3);%y
            f=ts*(t(1)^2+t(2)^2)-t(3)^2;%1
            %Choose y and solve for
            y=0:pixsize-1;
            qa=a*ones(size(y));
            qb=(c*y+d*ones(size(y)));
            qc=(b*y.^2+e*y+f*ones(size(y)));
            
            xp=((-qb+sqrt(qb.^2-4*qa.*qc))./qa*.5);
            xm=((-qb-sqrt(qb.^2-4*qa.*qc))./qa*.5);
            
            %If necessary choose x and solve for y
            if sum(abs(imag(xp)))>0
                x=0:pixsize-1;
                qa=b*ones(size(x));
                qb=(c*x+e*ones(size(x)));
                qc=(a*x.^2+d*x+f*ones(size(x)));
                
                yp=((-qb+sqrt(qb.^2-4*qa.*qc))./qa*.5);
                ym=((-qb-sqrt(qb.^2-4*qa.*qc))./qa*.5);
                
                %sort to find the high and low vals
                ymin=ceil(min([yp;ym]));
                ymax=floor(max([yp;ym]));
    %             ymin=round(min([yp;ym]));
    %             ymax=round(max([yp;ym]));
                %max sure they fall on the screen
                x(ymin>pixsize)=[];
                ymax(ymin>pixsize)=[];
                ymin(ymin>pixsize)=[];
                x(ymax<1)=[];
                ymin(ymax<1)=[];
                ymax(ymax<1)=[];
                ymin(ymin<1)=1;
                ymax(ymax>pixsize)=pixsize;
                for ind=1:length(x)
                    simpat((ymin(ind)):(ymax(ind)),x(ind)+1)=simpat((ymin(ind)):(ymax(ind)),x(ind)+1)+sFhkl(i);
                end
            else
                %sort to find the high and low vals
                xmin=ceil(min([xp;xm]));
    %             xmin=round(min([xp;xm]));
                xmax=floor(max([xp;xm]));
    %             xmax=round(max([xp;xm]));
                %make sure they fall on the screen
                y(xmin>pixsize)=[];
                xmax(xmin>pixsize)=[];
                xmin(xmin>pixsize)=[];
                y(xmax<1)=[];
                xmin(xmax<1)=[];
                xmax(xmax<1)=[];
                xmin(xmin<1)=1;
                xmax(xmax>pixsize)=pixsize;
                for ind=1:length(y)
                    % applies a guassian distribution to the main bands
                    
                    %                 if i < 4
                    %
                    %                     le = xmax(ind)-xmin(ind)+1;
                    %                     sincer = 1:le;
                    %                     sincer = sincer-length(sincer)*.5;
                    %                     sincer = exp(-sincer.^2/(2*(le*.5)^2));
                    %                     simpat(y(ind)+1,(xmin(ind)):(xmax(ind))) = simpat(y(ind)+1,(xmin(ind)):(xmax(ind)))+sFhkl(i)*sincer;
                    %
                    %                 else
                    %                     if ind == 1
                    %                     xmin(ind)
                    %                     xmax(ind)
                    %                     keyboard
                    %                     end
                    simpat(y(ind)+1,(xmin(ind)):(xmax(ind))) = simpat(y(ind)+1,(xmin(ind)):(xmax(ind))) + sFhkl(i);
                    
                    %                 end
                end
            end 
        end 
    end
%     warning('on', 'MATLAB:colon:operandsNotRealScalar')
    I = single(simpat);
end


function [NewP] = PermuteHKL(hkl,lattice)
    % Adaptation of Stuart Wright's TSL code for removing redundant symmetries
    % for a given hkl.
    % permutation of Miller indices according to symmetry
    % The permutations correspond only to TRUE rotations
    % no inversions - rotations as given in SymElements.h
    % This means the sum of the # of permutations and # of
    % negations must be even.
    % clear all; clc;
    % lattice = 'cubic'
    % hkl = [1 -1 -1]
    
    if ~strcmp(lattice,'cubic') && ~strcmp(lattice,'hexagonal')
        
        error('Permutations not supported for this crystal lattice type');
    end
    
    in = 1;
    
    %        int i,j,k,h,l,ii,jj,kk,ll,in=0,numpossible,tag;
    
    %        h = hkl(1); k = hkl(2); l = hkl(3); i = -(h+k);
    switch (lattice)
        
        case 'cubic'
            for l = 0:5
                
                if l < 3
                    i = l;
                    j = mod(i+1,3);
                    k = mod(i+2,3);
                    ll = 1;
                else
                    j = l-3;
                    i = mod(j+1,3);
                    k = mod(j+2,3);
                    ll = -1;
                end
                
                for ii = 1:-2:-1
                    for jj = 1:-2:-1
                        for kk = 1:-2:-1
                            
                            if ii*jj*kk*ll < 0
                                continue
                            end
                            p(in,1) = ii * hkl(i+1);
                            p(in,2) = jj * hkl(j+1);
                            p(in,3) = kk * hkl(k+1);
                            in = in + 1;
                            
                        end
                    end
                end
            end
            NewP = unique(p,'rows');
            KeepList = [];
            if strcmp(lattice,'cubic')
                keepcnt = 1;
                for gg = 1:size(NewP,1)
                    
                    BadEgg = 0;
                    for hh = gg+1:size(NewP,1)
                        %               [NewP(gg,:) ; NewP(hh,:)]
                        
                        if (NewP(gg,:) == -NewP(hh,:))
                            BadEgg = 1;
                        end
                    end
                    if ~BadEgg
                        KeepList(keepcnt) = gg;
                        keepcnt = keepcnt + 1;
                    end
                end
            end
            NewP(1:size(KeepList,1),:);
            
            
            
        case 'hexagonal'
            h = hkl(1); k = hkl(2); l = hkl(3); i = -(h+k);
            p(in,1) =  h; p(in,2) =  k; p(in,3) =  l; in = in + 1;
            p(in,1) =  k; p(in,2) =  i; p(in,3) =  l; in = in + 1;
            p(in,1) =  i; p(in,2) =  h; p(in,3) =  l; in = in + 1;
            p(in,1) = -h; p(in,2) = -k; p(in,3) =  l; in = in + 1;
            p(in,1) = -k; p(in,2) = -i; p(in,3) =  l; in = in + 1;
            p(in,1) = -i; p(in,2) = -h; p(in,3) =  l; in = in + 1;
            p(in,1) =  k; p(in,2) =  h; p(in,3) = -l; in = in + 1;
            p(in,1) =  i; p(in,2) =  k; p(in,3) = -l; in = in + 1;
            p(in,1) =  h; p(in,2) =  i; p(in,3) = -l; in = in + 1;
            p(in,1) = -k; p(in,2) = -h; p(in,3) = -l; in = in + 1;
            p(in,1) = -i; p(in,2) = -k; p(in,3) = -l; in = in + 1;
            p(in,1) = -h; p(in,2) = -i; p(in,3) = -l; in = in + 1;
            
            
            SymOps = gensymopsHex;
           
            
            for jj = 1:size(p,1)
                NewP(jj,:) = squeeze(SymOps(jj,1:3,1:3))*p(jj,:)';
            end
            
            
    end
end