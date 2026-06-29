classdef H5oinaPatternProvider < patterns.PatternProvider
    %H5PATTERNPROVIDER Provides patterns from an HDF5 file
    
    properties
        imSize double
    end
    
    properties(Access = private)
        patternPath(1,:) char
        patternSize
    end
    
    methods
        function obj = H5oinaPatternProvider(fileName)
            fileInfo = h5info(fileName);
%             assert(length(fileInfo.Groups) == 1,...
%                 ['More than one root group in the HDF5 file %s,\n'...
%                 'could not determine correct path to patterns.'], fileName);
            rootGroupName = fileInfo.Groups(1).Name;
            dataPath = [rootGroupName '/EBSD/Data'];
            patternPath = [dataPath '/Processed Patterns'];
            patternInfo = h5info(fileName, patternPath);
            patternSize = patternInfo.ChunkSize;
            
            obj@patterns.PatternProvider(fileName, min(patternSize(1:2)));
            obj.patternPath = patternPath;
            obj.patternSize = patternSize;
            obj.imSize = patternSize(1:2);
        end
        
        function sobj = saveobj(obj)
            sobj = saveobj@patterns.PatternProvider(obj);
            sobj.patternPath = obj.patternPath;
            sobj.patternSize = obj.patternSize;
        end
    end
    
    methods (Access = protected)
        function pattern = getPatternData(obj, ind)
            start = ones(size(obj.patternSize));
            start(obj.patternSize == 1) = ind;
            count = obj.patternSize;
            pattern = h5read(obj.fileName, obj.patternPath, start, count)';
        end
    end
    
    methods (Static)
        function obj = restore(loadStruct)
            obj = patterns.H5PatternProvider(loadStruct.fileName);
            obj.patternPath = loadStruct.patternPath;
            obj.patternSize = loadStruct.patternSize;
            obj.hexOffset = loadStruct.hexOffset;
        end  
    end
end