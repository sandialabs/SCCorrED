classdef frameTransformations
    %frameTransformations is a class containing static methods for
    %performing frame transformations during HREBSD analysis! 
    properties
        m
        Qp2s % transformation from phosphor to sample frames
        Qiv2pf = [-1, 0, 0; 0, -1, 0; 0, 0, 1];
%         Qiv2pf = [0, -1, 0; -1, 0, 0; 0, 0, 1];
    end
    
    methods 
        function obj = frameTransformations(mtexHREBSD)
            obj.m = mtexHREBSD.roi.pixelSize;
            obj.Qp2s = obj.get_Qps(mtexHREBSD);
        end

        function phosphorVec = imageVec2phosphorFrame(obj, imageVec)
            temp = imageVec;
            temp(3,:) = zeros(1,size(imageVec,2));
            phosphorVec = (obj.Qiv2pf*temp)./obj.m;
        end

        function imageVec = phosphorFrame2imageVec(obj, phosphorVec)
            temp = (obj.Qiv2pf'*phosphorVec);
            imageVec = temp(1:2,:).*obj.m;
        end
    end


    methods(Static)
        function Qps = get_Qps(mtexHREBSD)
            sampleTilt = mtexHREBSD.scan.sampleTilt * pi/180;
            if isprop(mtexHREBSD.scan, "detectorOrientation")
                Qmp = euler2gmat(mtexHREBSD.scan.detectorOrientation(1), ...
                                 mtexHREBSD.scan.detectorOrientation(2), ...
                                 mtexHREBSD.scan.detectorOrientation(3));
                Qmi = [ 0, -1, 0; 
                        1,  0, 0;
                        0,  0, 1];
                c = cos(sampleTilt);
                s = sin(sampleTilt);
                q_emsoft = [0,-1,0; -1,0,0; 0,0,-1];
                Qio = [c, 0, -s;
                       0, 1,  0;
                       s, 0,  c];
                Qps = double(q_emsoft'*Qio*Qmi*Qmp'*[-1, 0, 0; 0, 1, 0; 0, 0, -1]);
            else
                cameraElevation = mtexHREBSD.scan.cameraElevation * pi/180;
                alphaRotation = pi/2 - sampleTilt + cameraElevation;
                c = cos(alphaRotation);
                s = sin(alphaRotation);
                q_emsoft = [0,-1,0; -1,0,0; 0,0,-1];
                Qps = double(q_emsoft'*[ 0, -c, -s;...
                                 -1,  0,  0;...
                                  0,  s, -c]);
            end
        end


        function Qps = phosphor2sample(scanSettings)
            % add options for camera Euler angles
            sampleTilt = scanSettings.sampleTilt * pi/180;
            cameraElevation = scanSettings.cameraElevation * pi/180;
            alphaRotation = pi/2 - sampleTilt + cameraElevation;
            c = cos(alphaRotation);
            s = sin(alphaRotation);
            Qps = [ 0, -c, -s;...
                   -1,  0,  0;...
                    0,  s, -c];
        end


%         function phosphorVec = imageVec2phosphorFrame(imageVec)
%             phosphorVec = [-imageVec(2), -imageVec(1), 0];
%         end


%         function imageVec = phosphorFrame2imageVec(phorsphorVec)
%             imageVec = [-phorsphorVec(2), -phorsphorVec(1)];
%         end


        function Mi = rotate2imageFrame(Mp)
            Q = [0, -1, 0;
                 -1, 0, 0;
                 0, 0, 1];
            Mi = Q*Mp*Q';
        end
        

        function Mp = rotate2phosphorFrame(Mi)
            Q = [0, -1, 0;
                 -1, 0, 0;
                 0, 0, 1];
            Mp = Q'*Mi*Q;
        end


        function Pp = pc2phosphorFrame(P)
            Pp = [-P(1); -(1-P(2)); P(3)];
        end


        function P_EM = pc2EMsoft(P,m,n)
            P_EM = [m*(P(1) - 0.5), ...
                    m*P(2) - n/2, ...
                    m*P(3)];
        end


        function Pi = EMsoft2image(P,m,n)
            Pi = [P(2) + n/2, ...
                  P(1) + m/2, ...
                  P(3)];
        end
    end
end