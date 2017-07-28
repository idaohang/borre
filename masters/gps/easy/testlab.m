% test

x=[1 1.53 4];
y=[1 2 3];
plot(x,y)
set(gca,'XTick',x)
set(gca,'XTickLabel', sprintf('%3.4f|',x))
set(gca,'YTick',y)
set(gca,'YTickLabel', sprintf('%+1.2f|',y))