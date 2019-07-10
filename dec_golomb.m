function [symbol,i] = dec_golomb(i,bits, sign)

% i = 1;
length_M = 0;
x = 0; % x is a flag to exit when decoding of symbol is done

while x<1
    switch bits(i)
        case '1'
            if (length_M == 0)
                symbol = 0;
                i = i + 1;
                x = 1;
            else
                info = bin2dec(bits(i+1 : i+length_M));
                symbol = 2^length_M + info -1;
                i = i + length_M + 1;
                length_M = 0;
                x = 1;
            end
            
        case '0'
            length_M = length_M + 1;
            i = i + 1;
    end
end

if sign
        if symbol==0
        else
            symbol = (-1)^(symbol+1)*ceil(symbol/2);
        end
end

end