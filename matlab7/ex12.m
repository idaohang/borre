%ex12  demonstration of various figure settings
% written by Kai Borre
% October 7, 2004

figure;
axis;
grid
set((get(gcf,'children')),'color',[1 0 0],'gridlinestyle','--','linewidth',5)
set((get(gcf,'children')),'xcolor',[0 1 0],'ycolor',[0 0 1])
%%%%%%%%%%%%%%%%%%  end ex12.m %%%%%%%%%%%%%%%%%%%%%%