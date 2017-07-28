function nw13(triangles,nodes)
%NW13     Computes the structure matrix for each triangle
%         in a network. The structure matrix is a generalized
%         weight matrix derived according to the continous
%         network theory. The network is described by coordinates 
%         of the nodes and a list of the triangles. The 
%         coordinates are given in the file 'nodes' and the 
%         triangles are given in the file 'triangles'

%Kai Borre May 8, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

%Originally programmed in algol 5 on March 8, 1976

if nargin == 0
   triangles = 'triangle.dat';
   nodes = 'node_utm.dat';
end

tic
fidt = fopen(triangles,'r');
fidn = fopen(nodes,'r');

u = zeros(5,1);
sumt = 0;

node = [];
for j = 1:59
   node = [node; fscanf(fidn,'%f %f %f',[1 3])];
end
fclose(fidn);

fidd = fopen('net_anal.dat','w');
fidtt = fopen('triangle.tex','w');

point = node(:,1);
N = node(:,2);
E = node(:,3);

tri = [];
for j = 1:86
   tri = [tri; fscanf(fidt,'%g %g %g',[1 3])];
end
fclose(fidt);
tri(:,4) = tri(:,1);

for j = 1:86           % runs over triangles
   w = zeros(5,1);
   for i = 1:3         %triangle under consideration
      i_point = find(point == tri(j,i));
      X(i) = node(i_point,3); % easting
      Y(i) = node(i_point,2); % northing
   end
   X(4) = X(1);
   Y(4) = Y(1);
   for i = 1:3
      side(i) = norm([X(i+1)-X(i) Y(i+1)-Y(i)]);
   end   
   % area of actual triangle
   s = sum(side)/2;
   T = sqrt(s*(s-side(1))*(s-side(2))*(s-side(3)));
   %A nice procedure for finding neighbouring triangles     
   [r1,s1] = find(tri(j,1) == tri);
   [r2,s2] = find(tri(j,2) == tri);
   [r3,s3] = find(tri(j,3) == tri);
   n1 = intersect(r1,r2);
   n2 = intersect(r2,r3);
   n3 = intersect(r3,r1);
   nn1 = setxor(n1,j);
   if (isempty(nn1) == 1), n(1) = 0; else n(1) = nn1; end;
   nn2 = setxor(n2,j);
   if (isempty(nn2) == 1), n(2) = 0; else n(2) = nn2; end;
   nn3 = setxor(n3,j);
   if (isempty(nn3) == 1), n(3) = 0; else n(3) = nn3; end;
   sumn = 0;
   
   for i = 1:3  % for every triangle side
      if n(i) == 0
         TN = 0;
      else
         for k = 1:3
            i_point = find(point == tri(n(i),k));
            XN(k) = node(i_point,3); % easting
            YN(k) = node(i_point,2); % northing
         end
         XN(4) = XN(1);
         YN(4) = YN(1);          
         for k = 1:3
            siden(k) = norm([XN(k+1)-XN(k) YN(k+1)-YN(k)]);
         end
         sn = sum(siden)/2; 
         TN = sqrt(sn*(sn-siden(1))*(sn-siden(2))*(sn-siden(3)));
      end;
      sumn = sumn+TN/3;   
      theta = atan2(Y(i+1)-Y(i),X(i+1)-X(i));
      e = 6/(9*1.e-12+0.0004/(side(i))^2); 
      w(1) = w(1)+                e/(T+TN);
      w(2) = w(2)+   e*sin(2*theta)/(T+TN);
      w(3) = w(3)+.5*e*sin(4*theta)/(T+TN);
      w(4) = w(4)+   e*cos(2*theta)/(T+TN);
      w(5) = w(5)+.5*e*cos(4*theta)/(T+TN);
   end; % side i
   u = u+w*(T+sumn);
   sumt = sumt+2*T;
   fprintf(fidd,'\n %5.0f  %10.2f %10.2f %10.2f\n',...
                           tri(j,1), w(1), w(4), w(2));
   fprintf(fidd,' %5.0f \n %5.0f %22.2f %10.2f\n',...
                           tri(j,2),tri(j,3),w(5),w(3));   
   fprintf(fidtt,' %7.0f  %7.0f  %7.0f  %7.0f\n',...
                           tri(j,1),tri(j,2),tri(j,3),w(1));
end
u = u/sumt;
fprintf(fidd,'\n mean %12.2f %10.2f %10.2f\n',...
                                   u(1), u(4), u(2));
fprintf(fidd,'\n %28.2f %10.2f\n',u(5),u(3));   
fclose('all');
toc
%%%%%%%%%%%%%%%%%%% end nw13.m  %%%%%%%%%%%%%%%%%%%   