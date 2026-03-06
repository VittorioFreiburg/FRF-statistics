tic
% Test Prediction
f1=[0.0500    0.1500    0.3000    0.4000    0.5500    0.7000    0.9000    1.1000    1.3500    1.7500    2.2000];
load('SET1.mat');
sf=22;

 leout=20;

SET=SET1;
SET(leout,:)=[];
X=SET1(leout,:);
[xt,t]=pseudopulse(X,f1,sf);

%% Bootstrap
[avg,sigma,band,Cp,chist,values,alpha] = FRF_MinimalPredictionBand(X,SET,f1,1/22,1000);
%[avg,sigma,band,Cp,chist,values] = FRF_PredictionBand(SET,f1,1/22,0.95,4);

%% plots
ntrials=size(SET,1);
ns=441;

yt=zeros(ntrials,ns);

for i=1:ntrials %for the plot
    [x,t]=pseudopulse(SET(i,:),f1,sf);
    yt(i,:)=x;
end


% Histogram
figure(11)
 
plot(values(1:end-1),chist)


xlabel('$C_c$','interpreter','latex')
ylabel('$\max \left( \left| \hat{x}_b - \hat{x}) \right| / \hat{\sigma_b} \right)$','interpreter','latex','Fontsize',11)

hold on;
yline(alpha,'Linewidth',1.5,'Linestyle','-.','Color',[1 0 0],'Label',['\alpha^* = ',num2str(alpha*100,4),'%'],'Fontsize',10)
hold on;
xline(Cp,'Linewidth',1.5,'Linestyle','-.','Color',[1 0 1],'Label','C_p^*','Fontsize',10,'LabelVerticalAlignment','middle')


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print(gcf,'HistogramPredictionMinimal','-dpdf','-r0')



UPPER=band(1,:);
LOWER=band(2,:);
res=zeros(1,length(UPPER));
res(UPPER<0)=UPPER(UPPER<0);
res(LOWER>0)=LOWER(LOWER>0);

st=t(2)-t(1);
ns=length(res);
fd=fft(res); %back to frequency domain
fp=(0:(ns-1))/(st*(ns-1));
RF=zeros(1,11);
for f=1:11
    RF(f)=2*fd(fp==f1(f))/ns;
end


figure(12)

plot(t,avg,'color',[0 0 0],'linewidth',1.5); hold on
plot(t,UPPER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(t,xt,'linewidth',2,'color',[1 0 0]);
plot(t,yt,'color',[0.8 0.8 0.8]); hold on
plot(t,avg,'color',[0 0 0],'linewidth',1.5);
plot(t,UPPER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(t,LOWER,'color',[0.1 0.1 1],'linewidth',1.5,'linestyle','-.');
plot(t,xt,'linewidth',2,'color',[1 0 0]);
yline(0,'k-.');

legend({'$\hat{x}(t)$',[num2str(alpha*100,4),'%'],'test','$x_i(t)$'},'interpreter','latex');


xlabel('time (s)');
ylabel('PIR (°)');

title('Minimal Prediction Band Including the Sample')

% tested samples in FRF_t,yx



set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print(gcf,'TimeDomainMinPrediction','-dpdf','-r0')
toc
% pdf estimation
[cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(X,SET,f1,1/22,1000)
