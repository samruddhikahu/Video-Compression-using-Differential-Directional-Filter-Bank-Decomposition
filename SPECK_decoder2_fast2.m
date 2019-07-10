function [rim_final, P_im, Stat] = SPECK_decoder2_fast2(br, bit_plane_max, output, i_sizes)

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
% SPECK ALGORITHM (DECODER)
% 
% Xavier Alameda Pineda (xavi.alameda at gmail.com)
% Image and Video Processing Group
% Technical University of Catalonia
% http://gps-tsc.upc.es/imatge
% 
% Syntaxis
% ========
% rim = SPECK_decoder(options)
% 
% Input
% =====
% options ~ SPECK codec options. Type 'help SPECK_codec_options' for help.
%     
% Output
% ======
% rim ~ restored image
% 
% see also SPECK_encoder

% <><><><><><>
% INITIALITZE
% <><><><><><>
clear global
computing_time = clock;

% Setting global variables
global rrim LSP LNP LIP LIS LISp input input_counter

LSP = [];
LNP = [];
LIP = [];
LIS = [];
LISp = [];
%options = opt;

% LIS (LIP) is initialitzed with the S set (pixel)
if i_sizes(1,1)*i_sizes(1,2)==1
    % [row, column, parent_bit_plane]
    LIP = [1,1,-1];
else
    % [row, column, height, width, num_elements, parent_bit_plane]
    LIS = [1,1,i_sizes(1,1),i_sizes(1,2),i_sizes(1,1)*i_sizes(1,2),-1];
end

% % i_sizes computing
% i_sizes = zeros(size(sizes)-[1,0]);
% i_sizes(1,:) = sizes(1,:);
% for ii=2:size(i_sizes,1),
%     i_sizes(ii,:) = i_sizes(ii-1,:) + sizes(ii,:);
% end

% When does set I start?
i_row = i_sizes(1,1);
i_col = i_sizes(1,2);

% Input management variables
result = 0;


input = output;
input_counter = 0;
br = 0;
% % Do we have to truncate the bitstream?
% maxbits = floor(prod(i_sizes(end,:))*br);
% if maxbits < numel(input),
%     input = input(1:maxbits);
% end
% input_counter = 0;
%rrim = zeros(sum(sizes(1:end-1,:),1));
rrim = zeros(i_sizes(end,:));

pass = 1;
% <><><><><>
% ALGORITHM
% <><><><><>

% For each bit plane
bit_plane = bit_plane_max;
while bit_plane >= 0
    % Sorting pass
    % <><><><><><>

    % For each possible newly significant pixel
    if(~isempty(LIP))
        LIP = Process_pix(LIP,bit_plane);
    end
            
    if result < 0, break; end

    % For each possible set size in increasing order
    if(~isempty(LIS))
        result = Process_Set(LIS,bit_plane);
    end
    
%     for n = 1:size(LIS,1),
%         result = ProcessS(LIS(n,:), bit_plane);
%         if result < 0, break; end
%     end
    
    if result < 0, break; end
    
    % If I ~= empty set, ProcessI
    if i_row < i_sizes(end,1) && i_col < i_sizes(end,2)
        [result, i_row, i_col] = ProcessI(bit_plane, i_row, i_col, i_sizes);
    end
    
    % Restore LIS and provisional LIS
    LIS = LISp;
    LISp = [];
    
    if result < 0, break; end

    % Refinement pass
    % <><><><><><><><>

    % For each pixel in the LSP decode it's bit_plane-th bit
    if(~isempty(LSP))
        Lind = sub2ind(size(rrim),LSP(:,1),LSP(:,2));
        Val = rrim(Lind);
        s = sign(Val);
        bits = input(1:length(Val))';
        input(1:length(Val)) = [];
        %bits = bits';
        %bitsn = (bits=='0').*(-1);
        bitsop = ((bits=='0').*(-1)) + str2num(bits);
        Val = abs(Val) + bitsop.*(2^(bit_plane-1));
        %Val = s.*Val;
        rrim(Lind) = s.*Val;
    end
        
    if result < 0, break; end
        
    % Restore LSP and LNP
    LSP = cat(1,LSP,LNP);
    LNP = [];
    
    Stat(pass,1) = pass;
    Stat(pass,2) = 2^bit_plane;
%     Stat(pass,3) = mse;
%     Stat(pass,4) = psnr;
    
    P_im(:,:,pass) = rrim;
       
    pass = pass + 1;
    % Quantization step
    % <><><><><><><><><>
    bit_plane = bit_plane - 1;
    
end

rim_final = rrim;

clear global

return

function LIP = Process_pix(LIP,bit_plane)
    global input rrim LNP
    tmp_array = 2^bit_plane : 2^(bit_plane+1)-1;
    tmp_value = mean(tmp_array);
    %[rx, cx] = size(LIP);
    bits1 = input(1:size(LIP,1))';
    input(1:size(LIP,1)) = [];
    
    x = 1; ix = 1;
    while(x<=size(LIP,1))
        if(bits1(ix)=='1')
            bits1 = [bits1; input(1)];
            input(1) = [];
            bmx(x) = bin2dec(bits1(ix:ix+1)');
            x = x + 1; ix = ix + 2;
        else
            bmx(x) = bin2dec(bits1(ix));
%             bits1 = [bits1; input(1)];
%             input(1) = [];
            x = x + 1; ix = ix + 1;
        end
    end
    
    tmpx = round(((bmx>0).*2.5) - bmx);
    Lnd = sub2ind(size(rrim),LIP(:,1),LIP(:,2));
    rrim(Lnd) = tmpx'.*tmp_value;
    
    ax1 = LIP(:,1).*(bmx>0)';
    ax2 = LIP(:,2).*(bmx>0)';
    ax3 = LIP(:,3).*(bmx>0)';
    
    if(any(ax1))
        ax(:,1) = ax1(ax1~=0);
        ax(:,2) = ax2(ax2~=0);
        ax(:,3) = ax2(ax3~=0);
        
        LNP = [LNP; ax];
    end
    
    LIP1 = LIP(:,1).*(bmx==0)';
    LIP2 = LIP(:,2).*(bmx==0)';
    LIP3 = LIP(:,3).*(bmx==0)';
    LIP1 = LIP1(LIP1~=0);
    LIP2 = LIP2(LIP2~=0);
    LIP3 = LIP3(LIP3~=0);
    
    LIP = [];
    if(~isempty(LIP1))
        LIP(:,1) = LIP1;
        LIP(:,2) = LIP2;
        LIP(:,3) = LIP3;
    end


function result = Process_Set(LIS, bit_plane)

global LISp input
%ttx = LIS;
%ttx = input;

result = 0;

idx = find(LIS(:,5)<4);
if(~isempty(idx))
    LISs = LIS(idx,:);
    for ix = 1:size(LISs,1)
        result = Process_smallS(LISs(ix,:),bit_plane);
        if result < 0
            return;
        end
    end
end
    LIS(idx,:) = [];
    if(~isempty(LIS))
        bitx = input(1:size(LIS,1))';
        input(1:size(LIS,1)) = [];
        
        LISc1 = LIS(:,1).*str2num(bitx);
        LISc2 = LIS(:,2).*str2num(bitx);
        LISc3 = LIS(:,3).*str2num(bitx);
        LISc4 = LIS(:,4).*str2num(bitx);
        LISc5 = LIS(:,5).*str2num(bitx);
        LISc6 = LIS(:,6).*str2num(bitx);
        
        if(any(LISc1))
            clear id6;
            id6 = find(LISc1>0);
            
            LISc(:,1) = LISc1(LISc1~=0);
            LISc(:,2) = LISc2(LISc2~=0);
            LISc(:,3) = LISc3(LISc3~=0);
            LISc(:,4) = LISc4(LISc4~=0);
            LISc(:,5) = LISc5(LISc5~=0);
            %LISc(:,6) = LISc6(LISc6~=0);
            LISc(:,6) = LISc6(id6,:);
            
            for ii = 1:size(LISc,1)
                result = CodeS(LISc(ii,:), bit_plane);
                if result < 0
                    return;
                end
            end
        end
        
        LISp1 = LIS(:,1).*(~str2num(bitx));
        LISp2 = LIS(:,2).*(~str2num(bitx));
        LISp3 = LIS(:,3).*(~str2num(bitx));
        LISp4 = LIS(:,4).*(~str2num(bitx));
        LISp5 = LIS(:,5).*(~str2num(bitx));
        LISp6 = LIS(:,6).*(~str2num(bitx));
        
        if(any(LISp1))
            clear id6;
            id6 = find(LISp1>0);
            
            LISpx(:,1) = LISp1(LISp1~=0);
            LISpx(:,2) = LISp2(LISp2~=0);
            LISpx(:,3) = LISp3(LISp3~=0);
            LISpx(:,4) = LISp4(LISp4~=0);
            LISpx(:,5) = LISp5(LISp5~=0);
            %LISpx(:,6) = LISp6(LISp6~=0);
            LISpx(:,6) = LISp6(id6,:);
            
            for ii = 1:size(LISpx,1)
                LISp = add_in_sorted_list(LISp,LISpx(ii,:),5);
            end
            
        end
    end

% % Fill in case of a pixel-set
% if numel(position) ~= 3 && numel(position) ~= 6
%     disp('aturat');
%     position = [position 1 1 1];
% end

% if numel(position) == 3
%     num_pixels = 1;
% else
%     num_pixels = position(5);
%     if num_pixels == 1
%         position = [position(1),position(2),position(6)];
%     end
% end
return


function result = CodeS(position, bit_plane)

%global LIP LNP rrim

result = 0;

num_pixels = position(5);

% If the S set is small
if num_pixels <= 4
    result = Process_smallS(position, bit_plane);
    if result < 0
        return;
    end
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
            Process_Set(pos(ii,:),bit_plane);
        end
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
                pos(ii,3) = floot(position(3)/4);
            end
            pos(ii,4) = position(4);
            pos(ii,5) = pos(ii,3)*pos(ii,4);
            pos(ii,6) = bit_plane;
            % Process this set, analyzing its significance
            Process_Set(pos(ii,:),bit_plane);
        end
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
        
        Process_Set(pos,bit_plane);
%         for ii=1:4,
%             Process_Set(pos(ii,:),bit_plane);
%         end
    end
end

return

function result = Process_smallS(position, bit_plane)

global LIP LNP rrim

result = 0;
% for each pixel in the set
% temporal values for reconstruction
temp_array = 2^bit_plane : 2^(bit_plane+1)-1;
temp_value = mean(temp_array);

for ii=1:position(3),
    for jj=1:position(4),
        pos = [position(1)+ii-1,position(2)+jj-1,position(6)];
        switch input_bit%(pos,1)
            % If it's significant
            case '1'
                % Decode it's sign
                switch input_bit%(pos,2)
                    case '1'
                        rrim(pos(1),pos(2)) = -temp_value;
                    case '0'
                        rrim(pos(1),pos(2)) = temp_value;
                    case 'E'
                        return;
                end
                % Add to the LNP
                LNP = cat(1,LNP,pos);
                if numel(LIP)~=0
                    indexes = find(LIP(:,1)==pos(1) & LIP(:,2)==pos(2));
                    if numel(indexes) == 1
                        LIP(indexes,:) = [];
                    elseif numel(indexes) > 1
                        error('A pixel repeated in the LIP');
                    end
                end
                % If it's not significant
            case '0'
                % Add to the LIP
                if numel(LIP)~=0
                    if ~sum(double(LIP(:,1)==pos(1) & LIP(:,2)==pos(2)))
                        LIP = cat(1,LIP,pos);
                    end
                else
                    LIP = pos;
                end
            case 'E'
                % Bit-stream finished
                result = -1;
                return;
        end
    end % jj
end % ii

return

function [result, i_row, i_col] = ProcessI(bit_plane, i_row, i_col, i_sizes)

result = 0;

switch input_bit%([],0) 
    case '1'
        % If it's significant
        % Code I set
        [result, i_row, i_col] = CodeI(bit_plane, i_row, i_col, i_sizes);
    case 'E'
        % If bit stream is finished
        % Finish code
        result = -1;
    case '0'
end

return

function [result, i_row, i_col] = CodeI(bit_plane, i_row, i_col, i_sizes)

%global i_row i_col i_sizes bit_plane

result = 0;

if i_row < i_sizes(end,1)

    pos = zeros(3,4);

    next_i_row = i_sizes(find(i_sizes(:,1)==i_row)+1,1);
    next_i_col = i_sizes(find(i_sizes(:,2)==i_col)+1,2);

    % Identify new S sets
    pos(1,1:4) = [ i_row+1 ,       1 , next_i_row-i_row ,            i_col ];
    pos(2,1:4) = [       1 , i_col+1 ,            i_row , next_i_col-i_col ];
    pos(3,1:4) = [ i_row+1 , i_col+1 , next_i_row-i_row , next_i_col-i_col ];
    pos(:,5) = pos(:,3).*pos(:,4);
    pos(:,6) = bit_plane;
    
    result = Process_Set(pos, bit_plane);

%     for ii=1:3,
%         % Process this sets
%         result = Process_Set(pos(ii,:), bit_plane);
%         if result < 0
%             return
%         end
%     end

    i_row = next_i_row;
    i_col = next_i_col;
    
    % Process bit plane
    if i_row < i_sizes(end,1),
        [result, i_row, i_col] = ProcessI(bit_plane, i_row, i_col, i_sizes);
    end
end

return

function bit = input_bit

global input

% if label_class == 0
%     label = 1;
% else
%     label = get_label(positionlabel_mode,label_class);
% end

if numel(input) > 0
    bit = input(1);
    input(1) = [];
else
    bit = 'E';
end

return
