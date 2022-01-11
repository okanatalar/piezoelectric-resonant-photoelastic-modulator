%% Plotting |s_11| in dB
figure;
plot(freq/1e6,20*log10(abs(s_11)),'LineWidth',4,'Color',[0,0,0.5])
hold on
%scatter(freq/1e6,20*log10(abs(s_11)),'filled', 'MarkerFaceColor', [0 0 0.5], 'Linewidth',40)
xlabel('Frequency (MHz)','FontSize',40)
ylabel('|s_{11}| (dB)','FontSize',40)
set(gca,'FontSize',60)

xlim([2.75 3.25])
ylim([-1 0.1]) 
xticks([2.8 2.9 3.0 3.1 3.2])
set(gca,'linew',4)
%set(gca,'Color',[203 226 230]/255)



