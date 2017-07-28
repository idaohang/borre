function  [AtA,AtY] = normals(AtA,AtY,H,omc,var)
% NORMALS  Accumulates the contribution of one
%	   observation equation and adds it to
%	   the coefficient matrix AtA and the
%	   right side AtY. The accumulated result
%	   is outputted under the same name.

%   Written by Kai Borre
%   September 16, 1996

	AtY = AtY + H*omc/var;
	AtA = AtA + H*H'/var;

%%%%%%% end normals.m  %%%%%%%%%%%%%
