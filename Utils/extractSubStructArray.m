function [ subStruct ] = extractSubStructArray( initialStruct, indexes )
    % extractSubStructArray. For every field of struct (which should be one dimensional array) extract sub array with specified indexes.
    %
    %   [ subStruct ] = extractSubStructArray( initialStruct, indexes )
    %
    %   INPUT
    %       initialStruct 	struct with initial data;
    %       indexes       	array of indexes, which sould be extracted.
    %
    %   OUPUT
    %       subStract  	struct with fields, which corresponging initialStruct, but each field is a subarray with specified indexes
    %
    narginchk(2, 2);
    
    fields = fieldnames(initialStruct);
    
    for i = 1:length(fields)
        fieldValue = initialStruct.(fields{i});
        subStruct.(fields{i}) = fieldValue(indexes);
    end
end
