%% 0Vpp Processing

spectrum = fft(trace_data);

DC_0Vpp = abs(spectrum(1));

Amp_0Vpp = 4*abs(spectrum(49));

DoM_0Vpp = Amp_0Vpp/DC_0Vpp;


%% 6Vpp Processing

spectrum = fft(trace_data);

DC_6Vpp = abs(spectrum(1));

Amp_6Vpp = 4*abs(spectrum(49));

DoM_6Vpp = Amp_6Vpp/DC_6Vpp;

%% 8Vpp Processing

spectrum = fft(trace_data);

DC_8Vpp = abs(spectrum(1));

Amp_8Vpp = 4*abs(spectrum(49));

DoM_8Vpp = Amp_8Vpp/DC_8Vpp;

%% 10Vpp Processing

spectrum = fft(trace_data);

DC_10Vpp = abs(spectrum(1));

Amp_10Vpp = 4*abs(spectrum(49));

DoM_10Vpp = Amp_10Vpp/DC_10Vpp;

%% 12Vpp Processing

spectrum = fft(trace_data);

DC_12Vpp = abs(spectrum(1));

Amp_12Vpp = 4*abs(spectrum(49));

DoM_12Vpp = Amp_12Vpp/DC_12Vpp;

%% 14Vpp Processing

spectrum = fft(trace_data);

DC_14Vpp = abs(spectrum(1));

Amp_14Vpp = 4*abs(spectrum(49));

DoM_14Vpp = Amp_14Vpp/DC_14Vpp;

%% Plotting
%4000 is length of trace_data
Voltage_vector = [0 6 8 10 12 14];
figure;
plot(Voltage_vector.^2*90/12^2,[DC_0Vpp,DC_6Vpp,DC_8Vpp,DC_10Vpp,DC_12Vpp,DC_14Vpp]/4000,'o','MarkerEdgeColor',[0 0 0.5],'MarkerFaceColor',[0 0 .5],'MarkerSize',25)
%plot(Voltage_vector,20*log10(abs(s_11)),'.-','LineWidth',4,'Color',[0,0,0.5])
xlabel('RF Power (mW)','FontSize',40)
ylabel('DC (V)','FontSize',40)
set(gca,'FontSize',60)
xlim([-20 160])
ylim([-2 14]) 
%xticks([3.4 3.6 3.8 4.0 4.2])
%yticks([-25 -20 -15 -10 -5 0])
set(gca,'Box','on');
set(gca,'linew',4)

%% Plotting
%4000 is length of trace_data
Voltage_vector = [0 6 8 10 12 14];
figure;
plot(Voltage_vector.^2*90/12^2,[Amp_0Vpp,Amp_6Vpp,Amp_8Vpp,Amp_10Vpp,Amp_12Vpp,Amp_14Vpp]/4000,'o','MarkerEdgeColor',[0 0 0.5],'MarkerFaceColor',[0 0 .5],'MarkerSize',25)
%plot(Voltage_vector,20*log10(abs(s_11)),'.-','LineWidth',4,'Color',[0,0,0.5])
xlabel('RF Power (mW)','FontSize',40)
ylabel('Beat Tone (V)','FontSize',40)
set(gca,'FontSize',60)
xlim([-20 160])
ylim([-0.2 2]) 
%xticks([3.4 3.6 3.8 4.0 4.2])
%yticks([-25 -20 -15 -10 -5 0])
set(gca,'Box','on');
set(gca,'linew',4)

%% Plotting
%4000 is length of trace_data
Voltage_vector = [0 6 8 10 12 14];
figure;
plot(Voltage_vector.^2*90/12^2,100*[DoM_0Vpp,DoM_6Vpp,DoM_8Vpp,DoM_10Vpp,DoM_12Vpp,DoM_14Vpp],'o','MarkerEdgeColor',[0 0 0.5],'MarkerFaceColor',[0 0 .5],'MarkerSize',25)
%plot(Voltage_vector,20*log10(abs(s_11)),'.-','LineWidth',4,'Color',[0,0,0.5])
xlabel('RF Power (mW)','FontSize',40)
ylabel('DoM (%)','FontSize',40)
set(gca,'FontSize',60)
xlim([-20 160])
ylim([-10 85]) 
%xticks([3.4 3.6 3.8 4.0 4.2])
%yticks([-25 -20 -15 -10 -5 0])
set(gca,'Box','on');
set(gca,'linew',4)

%% Plotting
figure;
time_vector = 0:1/(40e6):200/(40e6) - 1/(40e6);
plot(time_vector*1e6,trace_data(1:200),'LineWidth',4,'Color',[0,0,0.5])
xlabel('Time (\mus)','FontSize',40)
ylabel('Voltage (V)','FontSize',40)
set(gca,'FontSize',60)
xlim([0 5])
ylim([-0.01 0.12]) 
%xticks([3.4 3.6 3.8 4.0 4.2])
%yticks([-25 -20 -15 -10 -5 0])
set(gca,'Box','on');
set(gca,'linew',4)



