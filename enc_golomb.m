function [bits] = enc_golomb(symbol, sign)

% Zeroth Order Exponential Golomg Codes:-
bits = '';

% If signed_symbol flag is 1
if (sign)
    if (symbol ==0)
%         symbol = symbol;
    elseif (symbol>0)
        symbol = 2*symbol -1;
    else 
        symbol = (-2)*symbol;
    end
end

% Here code_num = symbol
% M is prefix, info is suffix
M = floor(log2(symbol + 1));
info = dec2bin(symbol + 1 - 2^M,M);

for j=1:M
    bits = [bits '0'];
end
bits = [bits '1'];
bits = [bits info];
    
