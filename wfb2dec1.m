function DWT = wfb2dec1(x, lvl, h, g)
% WFB2DEC   2-D Wavelet Filter Bank Decomposition
%
%       y = wfb2dec(x, h, g)
%
% Input:
%   x:      input image
%   h, g:   lowpass analysis and synthesis wavelet filters
%
% Output:
%   x_LL, x_LH, x_HL, x_HH:   Four 2-D wavelet subbands

% Make sure filter in a row vector
h = h(:)';
g = g(:)';

h0 = h;
len_h0 = length(h0);
ext_h0 = floor(len_h0 / 2);

% Highpass analysis filter: H1(z) = -z^(-1) G0(-z)
len_h1 = length(g);
c = floor((len_h1 + 1) / 2); 
% Shift the center of the filter by 1 if its length is even.
if mod(len_h1, 2) == 0
    c = c + 1;
end
h1 = - g .* (-1).^([1:len_h1] - c);
ext_h1 = len_h1 - c + 1;

[m, n] = size(x);
%x_LL = x;

for i = 1:lvl
    
    clear x_L x_H x_LH x_HL x_HH;
    % Row-wise filtering
    x_L = rowfiltering(x, h0, ext_h0);
    x_L = x_L(:, 1:2:end);
    %x_L = x_L(:, 2:2:end);
    
    x_H = rowfiltering(x, h1, ext_h1);
    x_H = x_H(:, 1:2:end);
    %x_H = x_H(:, 2:2:end);
    
    % Column-wise filtering
    x_LL = rowfiltering(x_L', h0, ext_h0)';
    x_LL = x_LL(1:2:end, :);
    %x_LL = x_LL(2:2:end, :);
    
    x_LH = rowfiltering(x_L', h1, ext_h1)';
    x_LH = x_LH(1:2:end, :);
    %x_LH = x_LH(2:2:end, :);
    
    x_HL = rowfiltering(x_H', h0, ext_h0)';
    x_HL = x_HL(1:2:end, :);
    %x_HL = x_HL(2:2:end, :);
    
    x_HH = rowfiltering(x_H', h1, ext_h1)';
    x_HH = x_HH(1:2:end, :);
    %x_HH = x_HH(2:2:end, :);
    
    a = 2^i; b = 2^(i-1);
    DWT(1:(m/a),1:(n/a)) = x_LL;
    DWT(((m/a)+1):(m/b),1:(n/a)) = x_HL;
    DWT(1:(m/a),((n/a)+1):(n/b)) = x_LH;
    DWT(((m/a)+1):(m/b),((n/a)+1):(n/b)) = x_HH;
    
    clear x;
    x = x_LL;
    
end


% Internal function: Row-wise filtering with border handling 
function y = rowfiltering(x, f, ext1)
ext2 = length(f) - ext1 - 1;
x = [x(:, end-ext1+1:end) x x(:, 1:ext2)];
y = conv2(x, f, 'valid');