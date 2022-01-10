%Reading txt file for broad_abs(s_11)
%data = readtable('broad_abs(s_11).txt'); % requires 2013a or later
data = readtable('broad_abs(s_11)_dB_1e-3_damping.txt'); % requires 2013a or later
freq_broad = table2array(data(:,1)); 
s_11_broad = table2array(data(:,2));

%Reading txt file for focused_abs(s_11)
%data = readtable('focused_abs(s_11).txt'); % requires 2013a or later
data = readtable('focused_abs(s_11)_dB_1e-3_damping.txt'); % requires 2013a or later
freq_focused = table2array(data(:,1)); 
s_11_focused = table2array(data(:,2));

%% Plotting |s_11| in dB (Broad)
figure;
x = [3.76 3.78 3.78 3.76];
y = [-12.965 -12.965 0.94 0.94];
patch(x,y,[203 226 230]/255,'LineStyle','none')
hold on
plot(freq_broad/1e6,s_11_broad,'.-','LineWidth',4,'Color',[0,0,0.5])

xlabel('Frequency (MHz)','FontSize',40)
ylabel('|s_{11}| (dB)','FontSize',40)
set(gca,'FontSize',60)
xlim([3.4 4.2])
ylim([-13 1]) 
set(gca,'Box','on');
set(gca,'linew',4)

%% Plotting |s_11| in dB (Focused)
figure;
plot(freq_focused/1e6,s_11_focused,'.-','LineWidth',4,'Color',[0,0,0.5])
xlabel('Frequency (MHz)','FontSize',40)
ylabel('|s_{11}| (dB)','FontSize',40)
set(gca,'FontSize',60)
xlim([3.76 3.78])
ylim([-13 1]) 
set(gca,'linew',4)
set(gca,'Color',[203 226 230]/255)

%% Strain plot along y direction (at center of wafer - (0,0) coordinates)
%Using interpolation by using 10 points 
y_vector = 1e-6*[50:50:450];
y_vector = [0.01e-6,y_vector,499.99e-6];

Syz_vector = [7.12e-6,1.36e-5,1.92e-5,2.35e-5,2.62e-5,2.71e-5,2.62e-5,2.35e-5,1.92e-5,1.36e-5,7.08e-6];


y = interpft(Syz_vector,1000);
Syz_interp = y(1:910);
figure;

%%
y_interp = [0:500e-6/909:500e-6];
plot(y_interp*1e6,Syz_interp,'.-','LineWidth',4,'Color',[0.75,0,0])

%Plotting S_yz bar
Syz_bar = ones(size(Syz_interp))*mean(Syz_interp);



xlabel('y (\mum)','FontSize',40)
ylabel('S_{yz}','FontSize',40)
set(gca,'FontSize',60)
xlim([-10 510])
ylim([-0.1e-5 3.2e-5]) 
set(gca,'linew',4)
