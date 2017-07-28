% Mandelbrot fractal and different visualizations

% Written by Kai Borre
% October 4, 2000

renum = input('renum:  '); % Reads number of real points
imnum = input('imnum:  '); % Reads number of imaginary points

remin = -2; remax = 1;     % Defines which numbers to compute
immin = -1.5; immax = 1.5;

reval1 = linspace(remin,remax,renum);
imval1 = linspace(immin,immax,imnum);
[Reval, Imval] = meshgrid(reval1,imval1);
Imvalreal = Imval; Imval = Imval*i;
Cgrid = Reval+Imval;
for reind = 1:renum
   disp(['reind = ',int2str(reind)]); % Writes loop status
   %Iteration loop z(i+1) = (z(i))^2+c
   for imind = 1:imnum
      c = Cgrid(reind,imind);
      numc = 0;
      zold = 0.0+i*0.0;
      z = zold^2+c;
      while (abs(z) <= 2) & (numc < 100)
         numc = numc+1;
         zold = z;
         z = zold^2+c;
      end
      Mandelbrot(reind,imind) = numc;
   end
end

clf
subplot(2,2,1), mesh(reval1,imval1,Mandelbrot)
axis([-2 1 -1.5 1.5 0 100])

subplot(2,2,2), contour(reval1,imval1,Mandelbrot,100), grid

subplot(2,1,2), surf(Reval,Imvalreal,Mandelbrot)

view(2);
shading flat, colormap(flipud(jet))
colorbar
axis([-2 1 -1.5 1.5])

%%%%%%%%%%%%%%%%% end ex11.m  %%%%%%%%%%%%%%%%%%%%%%
