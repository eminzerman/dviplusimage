function P = pu_encode( L )
% Convert luminance to perceptually uniform code values that are
% backward-compatible with sRGB. 
%
% P = pu_encode( L )
%
% P - perceptually uniform code values
% L - luminance in cd/m^2 (must be in absolute, not relative units!)
%
% Refer to [1] for details. Please cite 
% the paper if you find our work useful. The paper and bib file can be
% found at the project website: 
% http://www.mpi-inf.mpg.de/resources/hdr/fulldr_extension/
%
% [1] T. O. Ayd{\i}n, R. Mantiuk, and H.-P. Seidel, "Extending quality 
% metrics to full dynamic range images," in Human Vision and Electronic "
% Imaging XIII, Proceedings of SPIE, (San Jose, USA), January 2008
% 

persistent P_lut;
persistent l_lut;

if( isempty( P_lut ) ) % caching for better performance

    data = dlmread('pu_space.csv', ',');

    l_lut = log10( data(:,1) );
    P_lut = data(:,2);  
    
end

epsilon = 1e-4;

l_min = 10^(l_lut(1)+epsilon);
l_max = 10^(l_lut(end)-epsilon);

l = log10(max(min(L,l_max),l_min));

P = interp1( l_lut, P_lut, l );
