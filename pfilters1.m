function [h, g] = pfilters1(fname)
% PFILTERS    Generate filters for the Laplacian pyramid
%
%	[h, g] = pfilters(fname)
%
% Input:
%	fname:	Name of the filters, including the famous '9-7' filters
%		    and all other available from WFILTERS in Wavelet toolbox
%
% Output:
%	h, g:	1D filters (lowpass for analysis and synthesis, respectively)
%		    for seperable pyramid

switch fname
    case {'9-11', '9/11'}
        h = [0.00647208691208, -0.03125, -0.05177669529664, 0.28125];
        h = [h, 0.59060921676911, fliplr(h)];
        
        g = [0.00134041308792, 0.00647208691208, -0.00938289161544...
            , -0.05177669529664, 0.25804247852752];
        g = [g, 0.59060921676911, fliplr(g)];
    
    case {'13-19', '13/19'}
        h = [-0.00080901086401, 0, 0.01456219555218, -0.03125, -0.05096768443263, 0.28125];
        h = [h, 0.57442899948892, fliplr(h)];
        
        g = [0.00002094395450, 0, -0.00056548677147, -0.00080901086401, 0.00452389417173,...
            0.01456219555218, -0.04163820143138, -0.05096768443263, 0.28664758649661];
        g = [g, 0.57442899948892, fliplr(g)];
        
    case {'9-7', '9/7'}
	    h = [.037828455506995 -.023849465019380 -.11062440441842 ...
	         .37740285561265];	
	    h = [h, .85269867900940, fliplr(h)];
    
	    g = [-.064538882628938 -.040689417609558 .41809227322221];
	    g = [g, .78848561640566, fliplr(g)];
        
    case {'5-3', '5/3'}
	    h = [-1, 2, 6, 2, -1] / (4 * sqrt(2));
        g = [1, 2, 1] / (2 * sqrt(2));
        
    case {'Burt'}
	    h = [0.6, 0.25, -0.05];
	    h = sqrt(2) * [h(end:-1:2), h];
	
	    g = [17/28, 73/280, -3/56, -3/280];
	    g = sqrt(2) * [g(end:-1:2), g];
 	    
    case {'pkva'}	% filters from the ladder structure	
	    % Allpass filter for the ladder structure network
	    beta = ldfilter(fname);
	
	    lf = length(beta);
	    n = lf / 2;
	
	    if n ~= floor(n)
	        error('The input allpass filter must be even length');
        end
	
	    % beta(z^2)
	    beta2 = zeros(1, 2*lf-1);
	    beta2(1:2:end) = beta;
	
	    % H(z)
	    h = beta2;
	    h(2*n) = h(2*n) + 1;
	    h = h / 2;
	
	    % G(z)
	    g = -conv(beta2, h);
	    g(4*n-1) = g(4*n-1) + 1;
	    g(2:2:end) = -g(2:2:end);
	
	    % Normalize
	    h = h * sqrt(2);
	    g = g * sqrt(2);

    otherwise
	    [h, g] = wfilters(fname, 'l');
end