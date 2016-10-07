function  res = getfielddefault (inStruct, fieldname, defaultval)
    % query non recursive structure array with default values.
    %
    % Input:
    % ------
    % inStruct,     enqueried struct
    % fieldname,    field asked for
    % defaultval,   default value to be returned
    %
    % Return:
    % -------
    % res,          enqueried value/defaultval
    %
    % contact:      yunpeng.wng@gmail.com
    %
    res = 0;
    f = fieldnames(inStruct);
    for i = 1:length(f)
        if (strcmp(f{i}, strtrim(fieldname)))
            res = getfield(inStruct, fieldname);
            return;
        else
            res = defaultval;
        end
    end
