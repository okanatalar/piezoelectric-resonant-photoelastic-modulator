%% 4-Frame Detection for single frequency of operation (20 MHz)
L = 2e-2;
d = 0.5:0.01:7.5;
QE = 0.7;
P = 1/5;
h = 6.63e-34;
c = 3e8;
lambda_optical = 905e-9;
Pixel_No = 1e6;
r = 0.1;
f_s = 60;
lambda = c/20e6;
Mod_Depth = 0.8;
sigma_squared = 36;%Fixed noise variance per frame

DC = 0.5.*(P.*r./Pixel_No).*pi.*L.^2./(2.*pi.*d.^2).*QE./(h.*c./lambda_optical).*1./(f_s);%Photons per second per pixel

sigma = (sqrt(DC + sigma_squared)./sqrt(2)).*lambda./(4.*pi.*DC.*Mod_Depth);%Distance error

figure;
plot(d,sigma*100,'LineWidth',4,'Color',[0 0 0.5])
xlabel('Target Distance (m)','fontweight','bold','FontSize',36)
ylabel({'Standard Deviation of';'Distance Estimate (cm)'},'fontweight','bold','FontSize',36)
title('20 MHz')
set(gca,'FontSize',36); 
set(gca,'linew',4)

%% 4-Frame Detection for single frequency of operation (3.77 MHz)
L = 2e-2;
d = 0.5:0.01:7.5;
QE = 0.7;
P = 1/5;
h = 6.63e-34;
c = 3e8;
lambda_optical = 905e-9;
Pixel_No = 1e6;
r = 0.1;
f_s = 60;
lambda = c/3.77e6;
Mod_Depth = 0.8;
sigma_squared = 36;%Fixed noise variance per frame

DC = 0.5.*(P.*r./Pixel_No).*pi.*L.^2./(2.*pi.*d.^2).*QE./(h.*c./lambda_optical).*1./(f_s);%Photons per second per pixel

sigma = (sqrt(DC + sigma_squared)./sqrt(2)).*lambda./(4.*pi.*DC.*Mod_Depth);%Distance error

figure;
plot(d,sigma*100,'LineWidth',4,'Color',[0 0 0.5])
xlabel('Target Distance (m)','fontweight','bold','FontSize',36)
ylabel({'Standard Deviation of';'Distance Estimate (cm)'},'fontweight','bold','FontSize',36)
title('3.77 MHz')
set(gca,'FontSize',36); 
set(gca,'linew',4)