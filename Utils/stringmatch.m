function match = stringmatch(leftStr, rightStr)
% STRINGMATCH  Returns match > 0 if string1 and string2 match (string2 can be a cell array of strings). match = 0 is returned if not match is found.
    match = sum(strcmp(leftStr, rightStr));
end