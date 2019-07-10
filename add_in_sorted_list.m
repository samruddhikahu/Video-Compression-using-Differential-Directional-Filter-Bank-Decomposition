function list = add_in_sorted_list(list,value,column)

%     This file is part of SPECK_codec.
% 
%     SPECK_codec is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     SPECK_codec is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with SPECK_codec.  If not, see <http://www.gnu.org/licenses/>.
%
% Xavier Alameda Pineda (xavi.alameda at gmail.com)
% Image and Video Processing Group
% Technical University of Catalonia
% http://gps-tsc.upc.es/imatge

if nargin < 3
    disp('Usage: list = add_in_sorted_list(list,value,column)');
end

% If it's an empty list --> add to the list
if size(list,1)==0
    list = value;
    return;
end

% If it's below all values in the list
if value(column) <= list(1,column)
    list = cat(1,value,list);
    return;
% If it's over all values in the list
elseif value(column) > list(end,column)
    list = cat(1,list,value);
    return;
end

% Application of Bolzano's theorem.
% In this way we can add an element in time O(log(n))
I = floor(size(list,1)/2);
index = I;
while 1
    i = rem(I,2);
    I = floor(I/2);
    % If it's over the low-point
    if value(column) > list(index,column)
        % If it's below the next point
        if value(column) <= list(index+1,column)
            list = cat(1,list,zeros(size(value)));
            list(index+2:end,:) = list(index+1:end-1,:);
            list(index+1,:) = value;
            break;
        else
            index = index + 1;
        end
        index = index + I;
        if index > size(list,1)
            index = size(list,1);
        end
    elseif value(column) < list(index,column)
        if value(column) >= list(index-1,column)
            list = cat(1,list,zeros(size(value)));
            list(index+1:end,:) = list(index:end-1,:);
            list(index,:) = value;
            break;
        else
            index = index - 1;
        end
        index = index - I;
        if index < 1
            index = 1;
        end
    else
        list = cat(1,list,zeros(size(value)));
        list(index+2:end,:) = list(index+1:end-1,:);
        list(index+1,:) = value;
        break;
    end
end
