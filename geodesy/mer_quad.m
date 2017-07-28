function B = mer_quad(Qm,n,phi)
%MER_QUAD Computes the length of the meridian from Equator to the point
%         with latitude phi.

%Kai Borre 10-13-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0  $  $Date:1998/10/23  $

B = Qm*(phi+(-3/2*n+9/16*n^3)*sin(2*phi)...
          	  +(15/16*n^2-15/32*n^4)*sin(4*phi)...
	            -35/48*n^3*sin(6*phi)+(315/512*n^4)*sin(8*phi));
%%%%%%%%%%%%%%%%%%%%% end mer_quad.m %%%%%%%%%%%%%%
