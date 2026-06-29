function iniWrite(obj, kwargs)
% iniWrite writes a .ini config file to save mtexHREBSD settings
    %
    % Syntax
    %   mtexHREBSD.iniWrite(Name,Value)
    %
    % Name-Value Arguments
    %   File - where to write .ini file
    %       [] (default)|char|string
    %       Will enumerate filename is File already exists. If empty, 
    %       saves as "<pwd>/hrebsdConfig.ini"
    arguments
        obj mtexHREBSD
        kwargs.File {...
                    mustBeA(...
                        kwargs.File,...
                        ["string", "char", "double"]...
                    )} = []
    end
    filename = getFile(kwargs.File);
    disp("Writing .ini file to " + filename);
    fid=fopen(filename,'w');
    cleanupfid = onCleanup(@() fclose(fid));
    writeFirstBlock(obj, fid);
    writeBlock(obj,fid, "analysis");
    writeBlock(obj,fid, "scan");
    writeBlock(obj,fid, "ft");
    writeBlock(obj,fid, "roi");
    writeBlockPatterns(obj,fid);
end


function writeBlockPatterns(obj, fid)
    fprintf(fid, "[patterns]\n");
    fprintf(fid, "doCropSquare = %d\n", obj.patterns.doCropSquare);
    fprintf(fid, "\n");
    fprintf(fid, "[patterns.filter]\n");
    fprintf(fid, "lowerRadius = %d\n", obj.patterns.filter.lowerRadius);
    fprintf(fid, "upperRadius = %d\n", obj.patterns.filter.upperRadius);
    fprintf( ...
        fid, ...
        "lowerSmoothing = %d\n", ...
        obj.patterns.filter.lowerSmoothing ...
        );
    fprintf( ...
        fid, ...
        "upperSmoothing = %d\n", ...
        obj.patterns.filter.upperSmoothing ...
        );
    fprintf(fid, "doFilter = %d\n",obj.patterns.filter.doFilter);
    fprintf("\n")
end


function writeBlock(obj, fid, property) 
    blacklist = {'centers', 'custfilt', 'windowfunc'};
    keys = properties(obj.(property));
    fprintf(fid, "["+property+"]\n");
    for i = 1:length(keys)
        if ~any(strcmp(keys{i}, blacklist))
            fmtstr = string(keys{i}) + " = ";
            fmtstr = addFormatting(fmtstr, obj.(property).(keys{i}));
            fprintf( ...
                fid, ...
                fmtstr, ...
                obj.(property).(keys{i})(:) ...
            );
        end
    end
    fprintf(fid, "\n");
end


function fmtType = checkValueType(values)
    valueType = class(values);
    if any(strcmp(valueType, {'char', 'string'}))
        fmtType = "%s";
    else 
        fmtType = "%d";
    end
end


function fmtstr = addFormatting(fmtstr,values)
    fmtType = checkValueType(values); 
    if fmtType == "%s"
        values = string(values);
    end
    if length(values(:)) == 1
        fmtstr = fmtstr + fmtType + "\n";
    else
        fmtstr = fmtstr + "[";
        for i = 1:length(values(:))-1
            fmtstr = fmtstr + fmtType + ", ";
        end
        fmtstr = fmtstr + fmtType + "]\n";
    end
end


function writeFirstBlock(obj, fid)
    fprintf(fid,"[mtexHREBSD]\n");
    fprintf(fid, "version = ""%s""\n", obj.version);
    if isempty(obj.patternCenterOffset)
        fprintf(fid, "patternCenterOffset = []\n");
        fprintf(fid, "patternCenterOffsetMatrix = []\n");
    else
        fprintf( ...
            fid, ...
            "patternCenterOffset = [%d, %d, %d]\n", ...
            obj.patternCenterOffset ...
            );
        fprintf( ...
            fid, ...
            "patternCenterOffsetMatrix = " +...
            "[%d, %d, %d; %d, %d, %d; %d, %d, %d]\n", ...
            obj.patternCenterOffsetMatrix ...
            );
    end
    fprintf(fid, "\n");
end


function filename = getFile(File)
    if ~isempty(File)
        filename = File;
    else 
        filename = fullfile(pwd, "hrebsdConfig.ini");
    end
    filename = enumerateFilenames(filename);
end


function newFilename = enumerateFilenames(filename)
    if isfile(filename)
        [path,name,ext] = fileparts(filename);
        for i = 1:100
            newname = name + "("+num2str(i)+")"+ext;
            newFilename = fullfile(path, newname);
            if ~isfile(newFilename)
                break
            end
        end
    else
        newFilename = filename;
    end
end