function match = string_match(leftStr, rightStr)
    % string_match. Check that two string is match.
    %
    %   Returns match > 0 if string1 and string2 match (string2 can be a cell array of strings). 
    %   match = 0 is returned if not match is found.
    %
    %   INPUT
    %       leftStr     first string;
    %       rightStr    second string.
    %
    %   OUTPUT
    %       match   result of matching.
    %
    match = sum(strcmp(leftStr, rightStr));
end
