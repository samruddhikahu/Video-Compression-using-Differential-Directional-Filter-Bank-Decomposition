function [bit_plane_max, bitstream, Stat] = SPECK_encoder1_fast2(im,i_sizes)

%
%
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
% SPECK ALGORITHM (ENCODER)
% 
% Xavier Alameda Pineda (xavi.alameda at gmail.com)
% Image and Video Processing Group
% Technical University of Catalonia
% http://gps-tsc.upc.es/imatge
% 
% Syntaxis
% ========
% SPECK_encoder(options,im,sizes)
% 
% Input
% =====
% options ~ SPECK codec options. Type 'help SPECK_codec_options' for help.
% im ~ input image
% sizes ~ wavelet subband sizes (see wavedec2 of MATLAB Wavelet Toolbox
% for details)
%
% File Output
% ===========
% output_data/filename_speck_data.mat :
%       bp_max ~ maximum bit plane
%       output ~ SPECK-encoded image (bitstream)
%       sizes ~ a copy of the input argument
% output_data/filename_speck_additional_data.mat :
%       labels ~ (if desired) one label for each bit, for the arithmetic coding
% output_data/filename_speck_options.mat :
%       options ~ a copy of the input argument
% raw_data/filename_speck_raw.txt :
%       the raw data is written inside
%
% see also SPECK_encoder

% <><><><><><>
% INITIALITZE
% <><><><><><>

% Setting global variables
% <><><><><><><><><><><><>
clear global
global LSP LNP LIP LIS LISp output image 

LSP = [];
LNP = [];
LIP = [];
LIS = [];
LISp = [];
image = im;
output = [];
%options = opt;

% LIS (LIP) is initialitzed with the S set (pixel)
if i_sizes(1,1)*i_sizes(1,2)==1
    % [row, column, parent_bit_plane]
    LIP = [1,1,-1];
else
    % [row, column, height, width, num_elements, parent_bit_plane]
    LIS = [1,1,i_sizes(1,1),i_sizes(1,2),i_sizes(1,1)*i_sizes(1,2),-1];
end

% When does set I start?
i_row = i_sizes(1,1);
i_col = i_sizes(1,2);

% Higher bit plane
if(max(abs(image(:)))==0)
    bit_plane_max = 0;
else
    bit_plane_max = floor(log2(max(abs(image(:)))));
end

% <><><><><>
% ALGORITHM
% <><><><><>

% For each bit plane
bit_plane = bit_plane_max;
pass = 1;

while bit_plane >= 0
    % Sorting pass
    % <><><><><><>

    % For each not yet significant pixel:-
    if(~isempty(LIP))
        LIP = Process_pix(LIP,bit_plane);
    end

    % For each not yet sifnificant set in size increasing order
    if(~isempty(LIS))
        Process_Set(LIS, bit_plane);
    end
%     for n = 1:size(LIS,1),
%         ProcessS(LIS(n,:),bit_plane);
%     end
    
    % If I ~= empty set, ProcessI
    if i_row < i_sizes(end,1) && i_col < i_sizes(end,2)
        [i_row, i_col] = ProcessI(bit_plane,i_row,i_col,i_sizes);
    end
        
    % Restore LIS and provisional LIS
    LIS = LISp;
    LISp = [];
            
    % Refinement pass
    % <><><><><><><><>
    
    % For each pixel in the LSP code its bit_plane-th bit
    if(~isempty(LSP))
        %for m=1:size(LSP,1)
        clear Lind coeffx Val coeff;
        Lind = sub2ind(size(image),LSP(:,1),LSP(:,2));
        Val = image(Lind);
        coeffx = dec2bin(abs(Val));
        coeff = coeffx(:,end-bit_plane)';
        %coeff = (dec2bin(abs(image(LSP(m,1),LSP(m,2)))));
        output = [output coeff];
        %labels = [labels get_label(LSP(m,:),label_mode,3)];
        %end
    end
    
    % Restore LSP and LNP
    LSP = cat(1,LSP,LNP);
    LNP = [];
    
    Stat(pass,1) = pass;
    Stat(pass,2) = 2^bit_plane;
    Stat(pass,3) = length(output);
    
    pass = pass + 1;
    % Quantization step
    % <><><><><><><><><>
    bit_plane = bit_plane - 1;
end

bitstream = output;

function LIP = Process_pix(LIP, bit_plane)
    global image output LNP;
    Lind = sub2ind(size(image),LIP(:,1),LIP(:,2));
    Val = image(Lind);
    op = abs(Val)>=(2^bit_plane);
    Val_sgn = (op==1).*Val;
    ops = (Val_sgn>0);
    op_out = op + ops;
    op_out = ((op_out==1).*2) + op_out;
    
    for ix = 1:length(op_out)
        output = [output dec2bin(op_out(ix))];
    end
    if(any(op))
        ax1 = LIP(:,1).*op;
        ax2 = LIP(:,2).*op;
        ax3 = LIP(:,3).*op;
        
        ax(:,1) = ax1(ax1~=0);
        ax(:,2) = ax2(ax2~=0);
        ax(:,3) = ax2(ax3~=0);
        
        LNP = [LNP; ax];
    end
    
    if(any(~op))
        LIP1 = LIP(:,1).*(~op);
        LIP2 = LIP(:,2).*(~op);
        LIP3 = LIP(:,3).*(~op);
        LIP1 = LIP1(LIP1~=0);
        LIP2 = LIP2(LIP2~=0);
        LIP3 = LIP3(LIP3~=0);
        
        LIP = [];
        LIP(:,1) = LIP1;
        LIP(:,2) = LIP2;
        LIP(:,3) = LIP3;
    else
        LIP = [];
    end
    

function Process_Set(LIS, bit_plane)

global LISp image output

idx = find(LIS(:,5)<4);
if(~isempty(idx))
    LISs = LIS(idx,:);
    for ix = 1:size(LISs,1)
        Process_smallS(LISs(ix,:),bit_plane);
    end
    LIS(idx,:) = [];
end
xc = 0; xp = 0;
LISpx = []; LISc = [];
for ii = 1:size(LIS,1)
    clear temp_sub_image;
    temp_sub_image = image( LIS(ii,1) : LIS(ii,1) + LIS(ii,3) - 1,...
                            LIS(ii,2) : LIS(ii,2) + LIS(ii,4) - 1);
    if 2^bit_plane <= max(max(abs(temp_sub_image(:)))) % < 2^(bit_plane+1)
        % Output significance
        output = [output '1'];
        %labels = [labels get_label(position,label_mode,1)];
        % It's not necessary to remove from LIS, due to the provisional
        % LIS: LISp. If we don't add it to he LISp, it won't appear as
        % not significant in next significance scan steps.
        % However, we have to code it anyway:
        xc = xc + 1;
        LISc(xc,:) = LIS(ii,:);
        %CodeS(LIS(ii,:), bit_plane);
    % If it's not significant
    else
        % Output significance
        output = [output '0'];
        %labels = [labels get_label(position,label_mode,1)];
        % Add to the SORTED provisional LIS
        xp = xp + 1;
        LISpx(xp,:) = LIS(ii,:);
        %LISp = add_in_sorted_list(LISp,LIS(ii,:),5);
    end
end
for ic = 1:size(LISc,1)
    CodeS(LISc(ic,:),bit_plane);
end
if(~isempty(LISpx))
    for ip = 1:size(LISpx,1)
        LISp = add_in_sorted_list(LISp,LISpx(ip,:),5);
    end
end

function Process_smallS(position, bit_plane)

global LIP image LNP output

% ttx = image;
% ttx = LNP;

% for each pixel in the set
for ii=1:position(3),
    for jj=1:position(4),
        % Get its value
        pos = [position(1)+ii-1,position(2)+jj-1,position(6)];
        value = image(pos(1),pos(2));
        % If it's significant
        if abs(value)>=2^bit_plane
            % Output significance
            output = [output '1'];
            %labels = [labels get_label(pos,label_mode,1)];
            % Output its sign
            if value >= 0
                output = [output '0'];
                %labels = [labels get_label(pos,label_mode,2)];
            else
                output = [output '1'];
                %labels = [labels get_label(pos,label_mode,2)];
            end
            % Add to the LNP
            LNP = cat(1,LNP,pos);
            % Remove from LIP, if any
            if numel(LIP)~=0
                indexes = find(LIP(:,1)==pos(1) & LIP(:,2)==pos(2));
                if numel(indexes) == 1
                    LIP(indexes,:) = [];
                elseif numel(indexes) > 1
                    error('Pixel repeated in the LIP');
                end
            end
            % If it's not significant
        else
            % Output significance
            output = [output '0'];
            %labels = [labels get_label(pos,label_mode,1)];
            % Add to the LIP
            if numel(LIP)~=0
                if ~sum(double(LIP(:,1)==pos(1) & LIP(:,2)==pos(2)))
                    LIP = cat(1,LIP,pos);
                end
            else
                LIP = pos;
            end
        end
    end % jj
end % ii
    
function CodeS(position, bit_plane)

global LIP LNP output image

if numel(position) ~= 6
    disp('aturat');
end

num_pixels = position(5);

% If the S set is small
if num_pixels <= 4
    % for each pixel in the set
    for ii=1:position(3),
        for jj=1:position(4),
            pos = [position(1)+ii-1,position(2)+jj-1,bit_plane];
            % Get its value
            value = image(pos(1),pos(2));
            % If it's significant
            if abs(value)>=2^bit_plane
                % Output significance
                output = [output '1'];
                %labels = [labels get_label(pos,label_mode,1)];
                % Output its sign
                if value >= 0
                    output = [output '0'];
                    %labels = [labels get_label(pos,label_mode,2)];
                else
                    output = [output '1'];
                    %labels = [labels get_label(pos,label_mode,2)];
                end
                % Add to the LNP
                LNP = cat(1,LNP,pos);
            % If it's not significant
            else
                % Output significance
                output = [output '0'];
                %labels = [labels get_label(pos,label_mode,1)];
                % Add to the LIP
                LIP = cat(1,LIP,pos);
            end
        end % jj
    end % ii
% If the S set is big
else
    % Very excentric set
    % This set is an horitzontal rectangle. We split this set in four
    % rectangles:
    % ---------------------------------
    % |       |       |       |       |
    % |   1   |   2   |   3   |   4   |
    % |       |       |       |       |
    % ---------------------------------
    if position(4) >= 4*position(3)
        pos = zeros(4,6);
        II = rem(position(4),4);
        for ii=1:4
            pos(ii,1:2) = position(1:2) + (ii-1)*[0,floor(position(4)/4)];
            if ii<=II
                pos(ii,2) = pos(ii,2)+1;
                pos(ii,4) = floor(position(4)/4)+1;
            else
                pos(ii,4) = floot(position(4)/4);
            end
            pos(ii,3) = position(3);
            pos(ii,5) = pos(ii,3)*pos(ii,4);
            pos(ii,6) = bit_plane;
            % Process this set, analyzing its significance
        end
        ProcessS(pos, bit_plane);
    % Very excentric set
    % This is a vertical rectangle. We split this set in an analog way of
    % the last case.
    elseif position(3) >= 4*position(4);
        pos = zeros(4,5);
        II = rem(position(3),4);    
        for ii=1:4
            pos(ii,1:2) = position(1:2) + (ii-1)*[floor(position(3)/4),0];
            if ii<=II
                pos(ii,1) = pos(ii,1)+1;
                pos(ii,3) = floor(position(3)/4)+1;
            else
                pos(ii,3) = floor(position(3)/4);
            end
            pos(ii,4) = position(4);
            pos(ii,5) = pos(ii,3)*pos(ii,4);
            pos(ii,6) = bit_plane;
            % Process this set, analyzing its significance
        end
        Process_Set(pos, bit_plane);
    else
        pos = zeros(4,5);
        half_rows = position(3)/2;
        half_cols = position(4)/2;
        % Definition of the four subsets in the following way
        % ----------------
        % |       |      |
        % |   1   |   3  |
        % |       |      |
        % ----------------
        % |       |      |
        % |   2   |   4  |
        % ----------------
        pos(1,1) = position(1);
            pos(1,2) = position(2);
                pos(1,3) = ceil(half_rows);
                    pos(1,4) = ceil(half_cols);
                        pos(1,5) = pos(1,3)*pos(1,4);
        pos(2,1) = position(1) + ceil(half_rows);
            pos(2,2) = position(2);
                pos(2,3) = floor(half_rows);
                    pos(2,4) = ceil(half_cols);
                        pos(2,5) = pos(2,3)*pos(2,4);
        pos(3,1) = position(1);
            pos(3,2) = position(2) + ceil(half_cols);
                pos(3,3) = ceil(half_rows);
                    pos(3,4) = floor(half_cols);
                        pos(3,5) = pos(3,3)*pos(3,4);
        pos(4,1) = position(1) + ceil(half_rows);
            pos(4,2) = position(2) + ceil(half_cols);
                pos(4,3) = floor(half_rows);
                    pos(4,4) = floor(half_cols);
                        pos(4,5) = pos(4,3)*pos(4,4);
                        
        pos(:,6) = bit_plane;
        
        %for ii=1:4,
            Process_Set(pos, bit_plane);
        %end
    end
end

return

function [i_row, i_col] = ProcessI(bit_plane,i_row,i_col,i_sizes)

global output image

temp_im = image;
temp_im(1:i_row,1:i_col) = 0;

if 2^bit_plane <= max(max(abs(temp_im(:)))) % < 2^(bit_plane+1)
    % Output '1'
    output = [output '1'];
    %labels = [labels 1];
    [i_row, i_col] = CodeI(bit_plane, i_row, i_col, i_sizes);
else
    % Output '0'
    output = [output '0'];
    %labels = [labels 1];
end

return

function [i_row, i_col] = CodeI(bit_plane, i_row, i_col, i_sizes)

%global i_row i_col i_sizes

% In case the I set is not empty
if i_row < i_sizes(end,1)

    pos = zeros(3,4);

    next_i_row = i_sizes(find(i_sizes(:,1)==i_row)+1,1);
    next_i_col = i_sizes(find(i_sizes(:,2)==i_col)+1,2);

    % Compute the three subbands
    pos(1,1:4) = [ i_row+1 ,       1 , next_i_row-i_row ,            i_col ];
    pos(2,1:4) = [       1 , i_col+1 ,            i_row , next_i_col-i_col ];
    pos(3,1:4) = [ i_row+1 , i_col+1 , next_i_row-i_row , next_i_col-i_col ];
    pos(:,5) = pos(:,3).*pos(:,4);
    pos(:,6) = bit_plane;

    % Process these subbands
    %for ii=1:3,
        Process_Set(pos,bit_plane);
    %end

    i_row = next_i_row;
    i_col = next_i_col;

    % In case the I set is not empty, process it
	if i_row < i_sizes(end,1),
        [i_row, i_col] = ProcessI(bit_plane,i_row,i_col,i_sizes);
	end
    
end

return
