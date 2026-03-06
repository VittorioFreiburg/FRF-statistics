function h=FRF_Plotbands(time,avg,band,color)
h=gcf;
plot(time,avg,'color',color/2,'linewidth',2); hold on;
plot(time,band(1,:),'color',color,'linewidth',1.5,'linestyle','-.');
plot(time,band(2,:),'color',color,'linewidth',1.5,'linestyle','-.');