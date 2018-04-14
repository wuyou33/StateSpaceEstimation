function [ clone ] = deep_clone(data_struct)
    % deep_clone. Make deep clone of input structure.
    %
    %   For more details see: http://stackoverflow.com/questions/247430/matlab-copy-constructor
    %
    %   [ clone ] = deep_clone( data_struct )
    %
    %   INPUT
    %       data_struct  input structure.
    %
    %   OUTPUT
    %       clone   cloned structure with same data.
    %
    fn = [];
    
    try
        fn = [fn; fieldnames(data_struct)];
    catch MEstruct
        throw(MEstruct)
    end
    
    if length(fn) ~= length(unique(fn))
        error('[ deep_clone::FieldsNotUniqu ] Field names must be unique');
    end
    
    c = [];
    try
        c = [c ; struct2cell(data_struct)];
    catch MEdata
        throw(MEdata);
    end
    
    clone = cell2struct(c, fn, 1);
end
