%ex13.m  Gift to 7th semester international, October 2003

load topo;
[x y z] = sphere(100);
s = surface(x,y,z,'facecolor','texture','cdata',topo);
set(s,'edgecolor','none','FaceAlpha','texture','AlphaData',topo)
set(s,'BackFaceLighting','unlit')
colormap(topomap1)
alpha('direct')
alphamap([.1;1])
axis off vis3d
campos([2 13 10])
camlight
lighting gouraud

%%%%%%%%% end  ex13.m  %%%%%%%%%%%% 