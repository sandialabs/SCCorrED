classdef PatternSimulator < dynamicprops
    %PatternSimulator samples individual patterns from a given master
    %pattern
    %   Detailed explanation goes here
    %   
    % Will Gilliland 2025-1-27

    properties
        delta
        detectorEuler
        filter
        h5file
        phase
        mLPNH
        mLPSH
        patDims
        pixelSize
        Qpc2spec
        sampletilt
    end

    methods
        function obj = PatternSimulator(input, Phase)
            arguments
                input mtexHREBSD
                Phase {...
                    mustBeA(...
                        Phase,...
                        ["string", "char"]...
                    )} = []
            end
            obj.delta = input.scan.delta;
            if isprop(input.scan,'detectorOrientation')
                obj.detectorEuler = input.scan.detectorOrientation;
            else
                obj.detectorEuler = [0,input.scan.cameraElevation+90,0]'*pi/180;
            end
            obj.filter = input.patterns.filter;
            obj.phase = Phase;
            obj.h5file = fullfile( ...
                input.paths.EMdataPath,sprintf('%s_EBSDmaster.h5',Phase) ...
                );
            obj.patDims = [input.scan.pixelSize,input.scan.pixelSize]; %input.patterns.imSize;
            obj.pixelSize = input.scan.pixelSize;
            obj.sampletilt = input.scan.sampleTilt;
            obj.arm_dynamic_pattern;
        end


        function obj = arm_dynamic_pattern(obj)
            % =========================================================== %
            % This function returns simulated EBSD pattern (dynamic scattering)
            %
            % Inputs
            %   h5file : name of the master EBSD pattern file (h5 format)
            %   detectorEuler : array [phi1,Phi,phi2], Bunge Euler angles 
            %                   providing the EBSD detector orientation 
            %                   (Oxford convention)
            %
            % Outputs
            %   mLPNH : master pattern Northern hemisphere (modified Lambert)
            %   mLPSH : master pattern Southern hemisphere (modified Lambert)
            %   Qpc2spec : precomputed coordinate transformation from EMsoft PC 
            %              frame to specimen frame
            %
            % Written by Thomas Bennett, 2025-01-24
            %
            % Updates
            %   Will Gilliland, 2025-01-27 built class around function
            %
            % =========================================================== %
        
            % Read the master patterns (modified Lambert projection)
            obj.mLPNH = h5read(obj.h5file,"/EMData/EBSDmaster/mLPNH");
            obj.mLPSH = h5read(obj.h5file,"/EMData/EBSDmaster/mLPSH");
        
            % Add together individual master patterns for complex crystals
            if length(size(obj.mLPNH)) == 4
                obj.mLPNH = sum(obj.mLPNH,4);
                obj.mLPSH = sum(obj.mLPSH,4);
            end
        
            % Sum the master patterns from each electron energy bin
            obj.mLPNH = sum(obj.mLPNH,3);
            obj.mLPSH = sum(obj.mLPSH,3);
            obj.mLPNH = rescale(obj.mLPNH(:,:,end));
            obj.mLPSH = rescale(obj.mLPSH(:,:,end));
            
            % Determine the coordinate transformation from the EMsoft PC 
            % frame to the specimen frame
            % EMsoft PC frame to detector frame 2
            R1 = rotation.byEuler([pi,pi,0]);
            % detector frame 1 to detector frame 2 (invert later)
            R2 = rotation.byEuler(obj.detectorEuler');
            % detector frame 1 to EMsoft specimen frame
            R3 = rotation.byEuler([0,obj.sampletilt,-90]*degree);    
            % Computer the overall overall transformation 
            obj.Qpc2spec = R3'*R2*R1';
        end


        function EBSP = get_dynamic_pattern(obj,pattern)
            % =================================================================== %
            % This function returns simulated EBSD pattern (dynamic scattering)
            %
            % Inputs
            %   mLPNH : master pattern Northern hemisphere (modified Lambert)
            %   mLPSH : master pattern Southern hemisphere (modified Lambert)
            %   Qpc2spec : precomputed coordinate transformation from EMsoft PC 
            %              frame to specimen frame
            %   pcx : pattern center x (pixels, EMsoft)
            %   pcy : pattern center y (pixels, EMsoft)
            %   L : detector distance (micrometers)
            %   delta : pixel width on the EBSD detector screen (micrometers)
            %   patDims : array, [pattern width, pattern height]
            %   phi1, Phi, phi2 : Bunge Euler angles defining crystal orientation
            %                     in degrees
            %
            % Outputs
            %   EBSP : simulated pattern image
            %
            % Written by Thomas Bennett, 2025-01-24
            %
            % Updates
            %
            % To-do
            %   - Add pattern deformation
            %   - Verify the correctness of the detectorEuler implementation
            %
            % =================================================================== %
        
            % Quick and dirty info from pattern class
            pcx = -(pattern.patternCenter(1) - 0.5) * obj.pixelSize;
            pcy = (pattern.patternCenter(2) - 0.5) * obj.pixelSize;
            L = pattern.patternCenter(3) * obj.pixelSize * obj.delta;
%             phi1 = pattern.rotations(1);
%             Phi = pattern.rotations(2);
%             phi2 = pattern.rotations(3);
            q_emsoft = [0,-1,0; -1,0,0; 0,0,-1];
            [phi1, Phi, phi2] = gmat2euler(pattern.g * q_emsoft);
%             simulatedPattern = pattern; % when this isn't here, it slowed down for some reason?

            % Get the size of the master pattern
            [nh_MP,nw_MP] = size(obj.mLPNH);
        
            % Calculate the transformation matrix (EMsoft PC to crystal frame)
            g = rotation.byEuler([phi1,Phi,phi2]); % radians
            Qtotal = g'*obj.Qpc2spec;
            
            % Number of pixels in the pattern
            n = prod(obj.patDims);
            % Construct phosphor screen pixel coordinates
            % (EMsoft PC frame, micrometers)
            [xs,ys] = meshgrid(0:1:(obj.patDims(1)-1),0:1:(obj.patDims(2)-1));
            xs = ((obj.patDims(1) - 1) / 2 - pcx - xs(:))*obj.delta;
            ys = ((obj.patDims(2) - 1) / 2 - pcy - ys(:))*obj.delta;
        
            % Normalize the coordinates and combine into one array
            D = sqrt(xs.^2 + ys.^2 + L^2);
            xyz = [xs ./ D,ys ./ D,-L ./ D];
            % Perform the coordinate transformation
            xyzNew = (xyz * Qtotal.matrix');
            x = xyzNew(:,1);
            y = xyzNew(:,2);
            z = xyzNew(:,3);
        
            % Allocate memeory for mapped coordinates
            X = zeros([n,1]);
            Y = zeros([n,1]);
            % Distinguish between northern and southern hemispheres
            isNorthernHemisphere = true([n,1]);
        
            % Map cartesian coordinates to the modified Lambert projection
            % (see P. Callahan and M. de Graef, Microsc. Microanal. 2013)
            c1 = 0.5 * sqrt(pi); % a constant
            for i = 1:n
                z_tmp = z(i); % store the z coordinate to avoid repeated indexing
                if z_tmp < 0 % make corrections if z is negative
                    z_tmp = -z_tmp;
                    isNorthernHemisphere(i) = false;
                end
                if z_tmp >=1; z_tmp = 1; end % z <= 1
                c2 = sqrt(2*(1-z_tmp)); % another constant
                if abs(y(i)) <= abs(x(i)) % case 1 (Eq 10a)
                    if x(i) < 0; c2 = -c2; end
                    X(i) = c2 * c1;
                    if x(i) == 0 % special case if x = y = 0
                        Y(i) = c2 * c1;
                    else
                        Y(i) = c2 * atan(y(i)/x(i)) / c1;
                    end
                else                      % case 2 (Eq 10b)
                    if y(i) < 0; c2 = -c2; end
                    X(i) = c2 * atan(x(i)/y(i)) / c1;
                    Y(i) = c2 * c1;
                end
            end
        
            % convert modified Lambert coordinates to fractional array indices
            X2 = (nw_MP-1) * 0.5 * (X / sqrt(0.5*pi) + 1);
            Y2 = (nh_MP-1) * 0.5 * (Y / sqrt(0.5*pi) + 1);
        
            % Allocate memory for the EBSD pattern
            EBSP = zeros(n,1);
        
            % Use the master pattern as a look-up table 
            % (This is slightly faster than MATLAB's interp2() for patterns under
            % 4 million pixels)
            for i=1:n
                x_tmp = X2(i);
                y_tmp = Y2(i);
                i_x = floor(x_tmp);
                i_y = floor(y_tmp);
                % Check bounds
                if i_x < 0
                    i_x = 0;
                elseif i_x > nw_MP-2
                    i_x = nw_MP-2;
                end
                if i_y < 0
                    i_y = 0;
                elseif i_y > nh_MP-2
                    i_y = nh_MP-2;
                end
                dx = x_tmp - i_x;
                dy = y_tmp - i_y;
                % Get intensities at each of the four nearest pixels
                if isNorthernHemisphere(i)
                    I_11 = obj.mLPNH(i_x+1,i_y+1);
                    I_12 = obj.mLPNH(i_x+1,i_y+2);
                    I_21 = obj.mLPNH(i_x+2,i_y+1);
                    I_22 = obj.mLPNH(i_x+2,i_y+2);
                else
                    I_11 = obj.mLPSH(i_x+1,i_y+1);
                    I_12 = obj.mLPSH(i_x+1,i_y+2);
                    I_21 = obj.mLPSH(i_x+2,i_y+1);
                    I_22 = obj.mLPSH(i_x+2,i_y+2);
                end
                % Bilinear interpolation
                EBSP(i) = (1-dx) * (I_11 * (1-dy) + I_12 * dy) + dx * ...
                        (I_21 * (1-dy) + I_22 * dy);
            end
            EBSP = imadjust(reshape(EBSP,fliplr(obj.patDims)),[],[],0.7);
%             EBSP = reshape(EBSP,fliplr(obj.patDims));
            EBSP = obj.filter.filterImage(EBSP);
            
%             simulatedPattern.image = EBSP;
        end
    end
end