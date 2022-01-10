%% Plotting |s_11| in dB
figure;
x = [3.7238 3.7438 3.7438 3.7238];
y = [-29.9 -29.9 0.9 0.9];
patch(x,y,[203 226 230]/255,'LineStyle','none')
hold on
plot(freq/1e6,20*log10(abs(s_11)),'.-','LineWidth',4,'Color',[0,0,0.5])
xlabel('Frequency (MHz)','FontSize',40)
ylabel('|s_{11}| (dB)','FontSize',40)
set(gca,'FontSize',60)
xlim([3.4 4.2])
ylim([-30 1]) 
xticks([3.4 3.6 3.8 4.0 4.2])
yticks([-25 -20 -15 -10 -5 0])
set(gca,'Box','on');
set(gca,'linew',4)
%Take from 3.4 MHz to 4.2 MHz

%% Plotting zoomed in version 
figure;
plot(freq/1e6,20*log10(abs(s_11)),'LineWidth',4,'Color',[0,0,0.5])
hold on
%scatter(freq/1e6,20*log10(abs(s_11)),'filled', 'MarkerFaceColor', [0 0 0.5], 'Linewidth',40)
xlabel('Frequency (MHz)','FontSize',40)
ylabel('|s_{11}| (dB)','FontSize',40)
set(gca,'FontSize',60)

xlim([3.7238 3.7438])
ylim([-13 0.1]) 
xticks([3.725 3.73 3.735 3.74])
set(gca,'linew',4)
set(gca,'Color',[203 226 230]/255)



