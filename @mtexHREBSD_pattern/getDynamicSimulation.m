function image = getDynamicSimulation(obj, Ftensor, applyDeformation)
%gen_emsoftPattern Summary of this function goes here
%   Detailed explanation goes here
%     if ~exist('Ftensor', 'var')
%         Ftensor = '1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,';
%     end
    % set these paths to correct local path for EMsoft!!
    pc_p = obj.patternCenter;
    EMdataPath = obj.SimData.EMdataPath;
    EMsoftPath = obj.SimData.EMsoftPath;
    ImageInd = obj.scanIndex;
    if strcmp(obj.material, 'Ni')
        Material = 'Nickel';
    else
        Material = obj.material;
    end
    L = pc_p(3) * obj.SimData.pixelSize * obj.SimData.delta;
    thetac = obj.SimData.cameraElevation;
    delta = obj.SimData.delta;
    numsx = obj.SimData.pixelSize;
    numsy = numsx;
    xpc = -(pc_p(1) - 0.5) * obj.SimData.pixelSize;
    ypc = (pc_p(2) - 0.5) * obj.SimData.pixelSize;
    omega = 0;
    alphaBD = 0;
    energymin = 5.0;
    energymax = 20.0;
    includebackground = 'n';
    anglefile = ['OpenXY_euler' num2str(ImageInd) '.txt']; 
    eulerconvention = 'tsl';
    if ~exist(fullfile(EMdataPath,sprintf('%s_EBSDmaster.h5',Material)),'file')
        masterfile = (sprintf('%s_EBSDmasterout.h5',Material));
    else
        masterfile = (sprintf('%s_EBSDmaster.h5',Material)); 
    end
    energyfile = masterfile;
    datafile = ['EBSDout_' num2str(ImageInd) '.h5'];
    bitdepth = '8bit';
    beamcurrent = 15; 
    dwelltime = 100;   
    binning = 1;       
%     applyDeformation = 'n';
%     Ftensor = '1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,';
    scalingmode = 'gam';
    gammavalue = 0.4;   
    maskpattern = 'n';
    nthreads = 1;
    q_emsoft = [0,-1,0; -1,0,0; 0,0,-1];
    [phi1, PHI, phi2]=gmat2euler(obj.g * q_emsoft);
    datafilepath = fullfile(EMdataPath,datafile);
    inputfile = fullfile(EMdataPath,['OpenXY_' num2str(ImageInd) '.nml']);
    %Write testeuler.txt file
    fid=fopen(fullfile(EMdataPath,anglefile),'w');
    
    fprintf(fid,'eu\n');
    fprintf(fid,'1\n');
    fprintf(fid,'%g,%g,%g\n',phi1*180/pi,PHI*180/pi,phi2*180/pi);% in degrees
    fclose(fid);
    
    cleanupAngles = onCleanup(@() delete(fullfile(EMdataPath,anglefile)));

    %Write EMEBSDexample.nml file
    fid=fopen(inputfile,'w');
    formatString = [
    ' &EBSDdata\n'...
    ... template file for the EMEBSD program
    ...
    ... distance between scintillator and illumination point [microns]
    ' L = %g,\n'...
    ... tilt angle of the camera (positive below horizontal, [degrees])
    ' thetac = %g,\n'...
    ... CCD pixel size on the scintillator surface [microns]
    ' delta = %g,\n'...50.0
    ... number of CCD pixels along x and y
    ' numsx = %g,\n'...
    ' numsy = %g,\n'...
    ... pattern center coordinates in units of pixels
    ' xpc = %g,\n'...
    ' ypc = %g,\n'...
    ... angle between normal of sample and detector
    '! omega = %g,\n'...
    ... transfer lens barrel distortion parameter
    ' alphaBD = %g,\n'...0.0
    ... energy range in the intensity summation [keV]
    ' energymin = %g,\n'...5.0
    ' energymax = %g,\n'...20.0
    ... include a realistic intensity background or not ...
    ' includebackground = ''%s'',\n'...'y'
    ... name of angle file (euler angles or quaternions); path relative to EMdatapathname
    ' anglefile = ''%s'',\n'...'testeuler.txt'
    ... does this file have only orientations ('orientations') or does it also have pattern center and deformation tensor ('orpcdef')
    ... if anglefiletype = 'orpcdef' then each line in the euler input file should look like this: (i.e., 15 floats)
    ...   55.551210  58.856774  325.551210  0.0  0.0  15000.0  1.00 0.00 0.00 0.00 1.00 0.00 0.00 0.00 1.00
    ...   <-   Euler angles  (degrees)  ->  <- pat. ctr.   ->  <- deformation tensor in column-major form->
    ' anglefiletype = ''orientations'',\n'...
    ... 'tsl' or 'hkl' Euler angle convention parameter
    ' eulerconvention = ''%s'',\n'...'tsl'
    ... name of EBSD master output file; path relative to EMdatapathname
    ' masterfile = ''%s'',\n'...'master.h5'
    ...
    ' energyfile = ''%s'',\n'...'master.h5'
    ... name of output file; path relative to EMdatapathname
    ' datafile = ''%s'',\n'...'EBSDout.h5'
    ... bitdepth '8bit' for [0..255] bytes; 'float' for 32-bit reals; '##int' for 32-bit integers with ##-bit dynamic range
    ... e.g., '9int' will get you 32-bit integers with intensities scaled to the range [ 0 .. 2^(9)-1 ];
    ... '17int' results in the intensity range [ 0 .. 2^(17)-1 ]
    ' bitdepth = ''%s'',\n'...'8bit'
     ... incident beam current [nA]
    ' beamcurrent = %g,\n'...'150.0'
    ... beam dwell time [micro s]
    ' dwelltime = %g,\n'...'100.0'
    ... include Poisson noise ? (y/n) (noise will be applied *before* binning and intensity scaling)
    ' poisson = ''n'',\n'...
    ... binning mode (1, 2, 4, or 8)
    ' binning = %u,\n'...'1'
    ... should we perform an approximate computation that includes a lattice distortion? ('y' or 'n')
    ... This uses a polar decomposition of the deformation tensor Fmatrix which results in
    ... an approcimation of the pattern for the distorted lattice; the bands will be very close
    ... to the correct position in each pattern, but the band widths will likely be incorrect.
    ' applyDeformation = ''%s'',\n'...'n'
    ... if applyDeformation='y' then enter the 3x3 deformation tensor in column-major form
    ... the default is the identity tensor, i.e., no deformation
    ' Ftensor = %s\n'...1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,
    ... intensity scaling mode 'not' = no scaling, 'lin' = linear, 'gam' = gamma correction
    ' scalingmode = ''%s'',\n'...'not',
    ... gamma correction factor
    ' gammavalue = %g,\n'...1.0,
    ... if the 'makedictionary' parameter is 'n', then we have the normal execution of the program
    ... if set to 'y', then all patterns are pre-processed using the other parameters below, so that
    ... the resulting dictionary can be used for static indexing in the EMEBSDDI program...
    ... these parameters must be taken identical to the ones in the EMEBSDDI.nml input file to have
    ... optimal indexing...
    ' makedictionary = ''n'',\n'...
    ... should a circular mask be applied to the data? 'y', 'n'
    ' maskpattern = ''%s'',\n'...'n',
    ... mask radius (in pixels, AFTER application of the binning operation)
    ' maskradius = %u,\n'...
    ... hi pass filter w parameter; 0.05 is a reasonable value
    ' hipassw = 0.05,\n'...
    ... number of regions for adaptive histogram equalization
    ' nregions = 10,\n'...
    ... number of threads (default = 1)
    ' nthreads = %u,\n'...1
    ' /\n'];
    fprintf(fid,formatString,L,thetac,delta,numsx,numsy,xpc,ypc,omega,...
        alphaBD,energymin,energymax,includebackground,anglefile,...
        eulerconvention,masterfile,masterfile,datafile,bitdepth,...
        beamcurrent,dwelltime,binning,applyDeformation,Ftensor,...
        scalingmode,gammavalue,maskpattern,floor(numsx/2),nthreads);
    fclose(fid);
    cleanupNamelist = onCleanup(@() delete(inputfile));
    currentDirectory = pwd;
    cd(EMdataPath);
%     disp(['"' fullfile(EMsoftPath,'EMEBSD') '" ' inputfile])
    [status, cmdout] = system(['"' fullfile(EMsoftPath,'EMEBSD') '" ' inputfile]);
    cd(currentDirectory);
    cleanupDataFile = onCleanup(@() delete(fullfile(EMdataPath,datafile)));
%     if status
%         disp(cmdout)
%     end
    %generate pic
%     disp(datafilepath)
    h5infostruct=h5info(datafilepath);
    data1=h5read(h5infostruct.Filename,'/EMData/EBSD/EBSDPatterns');
%     data1 = reshape(data1, [numsx, numsy]);
    image=zeros(numsx,numsy);
    image(:,:)=data1(:,:,1);
    image=(image');
end


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