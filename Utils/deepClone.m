function [ clone ] = deepClone( dataStruct )
    % deepClone. Make deep clone of input structure.
    %
    %   For more details see: http://stackoverflow.com/questions/247430/matlab-copy-constructor
    %
    %   [ clone ] = deepClone( dataStruct )
    %
    %   INPUT
    %       dataStruct  input structure.
    %
    %   OUTPUT
    %       clone   cloned structure with same data.
    %
    fn = [];
    
    try
        fn = [fn; fieldnames(dataStruct)];
    catch MEstruct
        throw(MEstruct)
    end
    
    if length(fn) ~= length(unique(fn))
        error('[ deepClone::FieldsNotUniqu ] Field names must be unique');
    end
    
    c = [];
    try
        c = [c ; struct2cell(dataStruct)];
    catch MEdata
        throw(MEdata);
    end
    
    clone = cell2struct(c, fn, 1);
end
