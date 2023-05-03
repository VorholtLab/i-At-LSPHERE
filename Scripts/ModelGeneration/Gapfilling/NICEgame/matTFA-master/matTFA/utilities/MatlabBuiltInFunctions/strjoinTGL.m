function joinedStr = strjoinTGL(str, delimiter)
%STRJOIN  Join cell array of strings into single string
%   S = STRJOIN(C) constructs the string S by linking each string within
%   cell array of strings C together with a space.
%
%   S = STRJOIN(C, DELIMITER) constructs S by linking each element of C
%   with the elements of DELIMITER. DELIMITER can be either a string or a
%   cell array of strings having one fewer element than C.
%
%   If DELIMITER is a string, then STRJOIN forms S by inserting DELIMITER
%   between each element of C. DELIMITER can include any of these escape
%   sequences:
%       \\   Backslash
%       \0   Null
%       \a   Alarm
%       \b   Backspace
%       \f   Form feed
%       \n   New line
%       \r   Carriage return
%       \t   Horizontal tab
%       \v   Vertical tab
%
%   If DELIMITER is a cell array of strings, then STRJOIN forms S by
%   interleaving the elements of DELIMITER and C. In this case, all
%   characters in DELIMITER are inserted as literal text, and escape
%   characters are not supported.
%
%   Examples:
%
%       c = {'one', 'two', 'three'};
%
%       % Join with space.
%       strjoin(c)
%       % 'one two three'
%
%       % Join as a comma separated list.
%       strjoin(c, ', ')
%       % 'one, two, three'
%
%       % Join with a cell array of strings DELIMITER.
%       strjoin(c, {' + ', ' = '})
%       % 'one + two = three'
%
%   See also STRCAT, STRSPLIT.

%   Copyright 2012-2014 The MathWorks, Inc.

    narginchk(1, 2);
    
    strIsString  = ischar(str);
    strIsCellstr = iscellstr(str);
    
    % Check input arguments.
    if ~strIsCellstr && ~strIsString
        error(message('MATLAB:strjoin:InvalidCellType'));
    end

    numStrs = numel(str);

    if nargin < 2
        delimiter = {' '};
    elseif ischar(delimiter)
        delimiter = {strescape(delimiter)};
    elseif iscellstr(delimiter) || ischar(delimiter)
        numDelims = numel(delimiter);
        if numDelims ~= 1 && numDelims ~= numStrs-1 
            error(message('MATLAB:strjoin:WrongNumberOfDelimiterElements'));
        elseif strIsCellstr && ischar(delimiter)
            delimiter = cellstr(delimiter);
        end
        delimiter = reshape(delimiter, numDelims, 1);
    else
        error(message('MATLAB:strjoin:InvalidDelimiterType'));
    end

    str = reshape(str, numStrs, 1);
    
    if strIsString
        joinedStr = join(str, delimiter);
    elseif numStrs == 0
        joinedStr = '';
    else
        joinedCell = cell(2, numStrs);
        joinedCell(1, :) = str;
        joinedCell(2, 1:numStrs-1) = delimiter;

        joinedStr = [joinedCell{:}];
    end
end
