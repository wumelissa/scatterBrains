function plotResults(dcs_sensitivity,scos_sensitivity,dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2,...
   dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2,intDist,fs_bfi)
%% Plots the final reults of the simulations
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)
%% Give the settings for the plots
pdiff=0.1;

% Calcution of the x and y limits for the plots
sensLim=[min([dcs_sensitivity,scos_sensitivity])-pdiff*range([dcs_sensitivity,scos_sensitivity]),...
    max([dcs_sensitivity,scos_sensitivity])+pdiff*range([dcs_sensitivity,scos_sensitivity])];
covLim=[min([dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2])*(1-pdiff*log10(range([dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2]))),...
    max([dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2])*(1+pdiff*log10(range([dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2])))];
cnrLim=[min([dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2])*(1-pdiff*log10(range([dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2]))),...
    max([dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2])*(1+pdiff*log10(range([dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2])))];
sdsLim=[min(intDist)-pdiff*range(intDist),max(intDist)+pdiff*range(intDist)];

% Settings for the plots to match colors and symbols
dcs_colors=[0,0,0;0,0.3059,0.5961];
scos_colors=[0,0,0;.2275,.6471,.4314;0,0.3059,0.5961];
symbols={'o','s','d'};
lw=5;
xylw=3;
ms=10;
fs=15;

%% Figure 1: Comparison of the cw and pwm DCS simulations
figure(1)
subplot(131)
hold off
plot(intDist,100*dcs_sensitivity,['-',symbols{1}],'color',dcs_colors(1,:),'markerfacecolor',dcs_colors(1,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(100*sensLim);
A=gca;
A.FontSize=fs;
grid on
grid minor
A.YAxis.LineWidth=xylw;
A.XAxis.LineWidth=xylw;
A.Box='off';
ylabel('Cerebral sensitivity (%)')
xlabel('Source-detector separation (mm)')

subplot(132)
hold off
semilogy(intDist,dcs_cov_cw,['-',symbols{1}],'color',dcs_colors(1,:),'markerfacecolor',dcs_colors(1,:),...
'linewidth',lw,'markersize',ms)
hold on
semilogy(intDist,dcs_cov_pwm,['-',symbols{3}],'color',dcs_colors(2,:),'markerfacecolor',dcs_colors(2,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(covLim);
B=gca;
B.FontSize=fs;
grid on
grid minor
B.YAxis.LineWidth=xylw;
B.XAxis.LineWidth=xylw;
B.Box='off';
ylabel(sprintf('Coefficient of variation @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')
title('Comparison of CW and PWM DCS')
legend({'CW','Pulse-wave modulated'},'fontsize',fs,'box','off','location','northwest')
subplot(133)
hold off
semilogy(intDist,dcs_cnr_cw,['-',symbols{1}],'color',dcs_colors(1,:),'markerfacecolor',dcs_colors(1,:),...
'linewidth',lw,'markersize',ms)
hold on
semilogy(intDist,dcs_cnr_pwm,['-',symbols{3}],'color',dcs_colors(2,:),'markerfacecolor',dcs_colors(2,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(cnrLim);
C=gca;
C.FontSize=fs;
grid on
grid minor
C.YAxis.LineWidth=xylw;
C.XAxis.LineWidth=xylw;
C.Box='off';
ylabel(sprintf('Contrast-to-noise ratio @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')

%% Figure 2: Comparison of the results from the SCOS simulations
figure(2)
subplot(131)
hold off
plot(intDist,100*scos_sensitivity,['-',symbols{1}],'color',scos_colors(1,:),'markerfacecolor',scos_colors(1,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(100*sensLim);
D=gca;
D.FontSize=fs;
grid on
grid minor
D.YAxis.LineWidth=xylw;
D.XAxis.LineWidth=xylw;
D.Box='off';
ylabel('Cerebral sensitivity (%)')
xlabel('Source-detector separation (mm)')

subplot(132)
hold off
semilogy(intDist,scos_cov_cw,['-',symbols{1}],'color',scos_colors(1,:),'markerfacecolor',scos_colors(1,:),...
'linewidth',lw,'markersize',ms)
hold on
semilogy(intDist,scos_cov_pwm1,['-',symbols{2}],'color',scos_colors(2,:),'markerfacecolor',scos_colors(2,:),...
'linewidth',lw,'markersize',ms)
semilogy(intDist,scos_cov_pwm2,['-',symbols{3}],'color',scos_colors(3,:),'markerfacecolor',scos_colors(3,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(covLim);
E=gca;
E.FontSize=fs;
grid on
grid minor
E.YAxis.LineWidth=xylw;
E.XAxis.LineWidth=xylw;
E.Box='off';
ylabel(sprintf('Coefficient of variation @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')
title('Comparison of CW and PWM SCOS');
legend({'CW','Pulse-wave modulated (reduced max power)','Pulse-wave modulated (reduced frame rate)'},'fontsize',fs*.75,'box','off','location','northwest')
subplot(133)
hold off
semilogy(intDist,scos_cnr_cw,['-',symbols{1}],'color',scos_colors(1,:),'markerfacecolor',scos_colors(1,:),...
'linewidth',lw,'markersize',ms)
hold on
semilogy(intDist,scos_cnr_pwm1,['-',symbols{2}],'color',scos_colors(2,:),'markerfacecolor',scos_colors(2,:),...
'linewidth',lw,'markersize',ms)
semilogy(intDist,scos_cnr_pwm2,['-',symbols{3}],'color',scos_colors(3,:),'markerfacecolor',scos_colors(3,:),...
'linewidth',lw,'markersize',ms)
xlim(sdsLim)
ylim(cnrLim);
F=gca;
F.FontSize=fs;
grid on
grid minor
F.YAxis.LineWidth=xylw;
F.XAxis.LineWidth=xylw;
F.Box='off';
ylabel(sprintf('Contrast-to-noise ratio @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')

%% Figure 3: Comparing the optimal implementation of SCOS and DCS
%  Optimal is determined by the CNR

% Calcultion of the optimal operating conditions
[dcs_cnr_opt,dcs_ind]=max([dcs_cnr_cw;dcs_cnr_pwm],[],1);
dcs_cov_opt=dcs_sensitivity./dcs_cnr_opt;
[scos_cnr_opt,scos_ind]=max([scos_cnr_cw;scos_cnr_pwm1;scos_cnr_pwm2],[],1);
scos_cov_opt=scos_sensitivity./scos_cnr_opt;

dcs_color_opt=[0.5412    0.3098    0.4902];
scos_color_opt=[0.9569    0.2667    0.1804];

lw2=1;
ms2=150;

figure(3)
subplot(131)
hold off
plot(intDist,100*dcs_sensitivity,'color',dcs_color_opt,'linewidth',lw)
hold on
plot(intDist,100*scos_sensitivity,'color',scos_color_opt,'linewidth',lw)
scatter(intDist(dcs_ind==1),100*dcs_sensitivity(dcs_ind==1),ms2,dcs_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(dcs_ind==2),100*dcs_sensitivity(dcs_ind==2),ms2,dcs_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==1),100*scos_sensitivity(scos_ind==1),ms2,scos_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==2),100*scos_sensitivity(scos_ind==2),ms2,scos_color_opt,'filled',symbols{2},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==3),100*scos_sensitivity(scos_ind==3),ms2,scos_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
xlim(sdsLim)
ylim(100*sensLim);
G=gca;
G.FontSize=fs;
grid on
grid minor
G.YAxis.LineWidth=xylw;
G.XAxis.LineWidth=xylw;
G.Box='off';
ylabel('Cerebral sensitivity (%)')
xlabel('Source-detector separation (mm)')

subplot(132)
hold off
semilogy(intDist,dcs_cov_opt,'color',dcs_color_opt,'linewidth',lw)
hold on
semilogy(intDist,scos_cov_opt,'color',scos_color_opt,'linewidth',lw)
scatter(intDist(dcs_ind==1),dcs_cov_opt(dcs_ind==1),ms2,dcs_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(dcs_ind==2),dcs_cov_opt(dcs_ind==2),ms2,dcs_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==1),scos_cov_opt(scos_ind==1),ms2,scos_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==2),scos_cov_opt(scos_ind==2),ms2,scos_color_opt,'filled',symbols{2},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==3),scos_cov_opt(scos_ind==3),ms2,scos_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
xlim(sdsLim)
ylim(covLim);
H=gca;
H.FontSize=fs;
grid on
grid minor
H.YAxis.LineWidth=xylw;
H.XAxis.LineWidth=xylw;
H.Box='off';
ylabel(sprintf('Coefficient of variation @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')
title('Comparison of DCS and SCOS');
legend({'DCS','SCOS'},'fontsize',fs,'box','off','location','northwest')

subplot(133)
hold off
semilogy(intDist,dcs_cnr_opt,'color',dcs_color_opt,'linewidth',lw)
hold on
semilogy(intDist,scos_cnr_opt,'color',scos_color_opt,'linewidth',lw)
scatter(intDist(dcs_ind==1),dcs_cnr_opt(dcs_ind==1),ms2,dcs_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(dcs_ind==2),dcs_cnr_opt(dcs_ind==2),ms2,dcs_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==1),scos_cnr_opt(scos_ind==1),ms2,scos_color_opt,'filled',symbols{1},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==2),scos_cnr_opt(scos_ind==2),ms2,scos_color_opt,'filled',symbols{2},'MarkerEdgeColor','k','LineWidth',lw2)
scatter(intDist(scos_ind==3),scos_cnr_opt(scos_ind==3),ms2,scos_color_opt,'filled',symbols{3},'MarkerEdgeColor','k','LineWidth',lw2)
xlim(sdsLim)
ylim(cnrLim);
I=gca;
I.FontSize=fs;
grid on
grid minor
I.YAxis.LineWidth=xylw;
I.XAxis.LineWidth=xylw;
I.Box='off';
ylabel(sprintf('Contrast-to-noise ratio @ %0.0f Hz',fs_bfi))
xlabel('Source-detector separation (mm)')
end