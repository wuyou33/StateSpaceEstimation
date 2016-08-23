function [ clone ] = deepClone( dataStruct )
    % deepClone
    %   Make deep clone of input structure
    %   http://stackoverflow.com/questions/247430/matlab-copy-constructor
    
    fn = [];
    
    try
        fn = [fn; fieldnames(dataStruct)];
    catch MEstruct
        throw(MEstruct)
    end    
    
    if length(fn) ~= length(unique(fn))
        error('deepClone:FieldsNotUnique', 'Field names must be unique');
    end
    
    c = [];
    try
        c = [c ; struct2cell(dataStruct)];
    catch MEdata
        throw(MEdata);
    end
    
    clone = cell2struct(c, fn, 1);
end
