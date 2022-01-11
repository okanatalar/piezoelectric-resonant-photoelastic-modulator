%% Monte Carlo Simulation for Distance Error Multi-Frequency Approach
L = 2e-2;
d = 0.5:0.01:80;
d_true = 0.5:0.1:80;
QE = 0.7;
P = 1/15;
h = 6.63e-34;
c = 3e8;
lambda_optical = 905e-9;
Pixel_No = 1e6;
r = 0.1;%Target reflectivity
f_s = 60;
Mod_Depth = 0.8;
sigma_squared = 36;%Fixed noise variance per frame


%% This assumes single camera with time multiplexing

%% This assumes camera number equal to mod.frequency with no time multiplexing
DC = 0.5.*(P.*r./Pixel_No).*pi.*L.^2./(2.*pi.*d_true.^2).*QE./(h.*c./lambda_optical).*1./(f_s);%Photons per frame per pixel - 0.5 due to polarizer
sigma = (sqrt(DC + sigma_squared)./sqrt(2))./(DC*Mod_Depth);%Phase error

c = 3e8;
lambda_1 = c/(11e6);
lambda_2 = c/(13e6);
lambda_3 = c/(20e6);
N = 3;%no frequencies used

estimation_standard_deviation = zeros(size(d_true));
repeat = 1000;

for ite = 1:length(d_true)
    for re=1:repeat   
        if N == 3
            phi_1_tilde = wrapTo2Pi(4*pi/lambda_1*d_true(ite) + sigma(ite)*randn);
            phi_2_tilde = wrapTo2Pi(4*pi/lambda_2*d_true(ite) + sigma(ite)*randn);
            phi_3_tilde = wrapTo2Pi(4*pi/lambda_3*d_true(ite) + sigma(ite)*randn);

            d_estimate = real(exp(-1j*4*pi/lambda_1*d).*exp(1j*phi_1_tilde) + exp(-1j*4*pi/lambda_2*d).*exp(1j*phi_2_tilde) + exp(-1j*4*pi/lambda_3*d).*exp(1j*phi_3_tilde));
        elseif N == 2
            phi_1_tilde = wrapTo2Pi(4*pi/lambda_1*d_true(ite) + sigma(ite)*randn);
            phi_2_tilde = wrapTo2Pi(4*pi/lambda_2*d_true(ite) + sigma(ite)*randn);

            d_estimate = real(exp(-1j*4*pi/lambda_1*d).*exp(1j*phi_1_tilde) + exp(-1j*4*pi/lambda_2*d).*exp(1j*phi_2_tilde));
        else
            phi_1_tilde = wrapTo2Pi(4*pi/lambda_1*d_true(ite) + sigma(ite)*randn);
            d_estimate = real(exp(-1j*4*pi/lambda_1*d).*exp(1j*phi_1_tilde));
        end
            
        [~,index] = max(d_estimate);
        distance_estimate = d(index);
        error = distance_estimate - d_true(ite);
        estimation_standard_deviation(ite) = estimation_standard_deviation(ite) + error^2;
    end
    estimation_standard_deviation(ite) = sqrt(estimation_standard_deviation(ite)/repeat);
end

%% Figure for all distance range

figure;
plot(d_true,estimation_standard_deviation,'LineWidth',4,'Color',[0 0 0.5])
xlabel('Target Distance (m)','fontweight','bold','FontSize',48)
ylabel({'Standard Deviation of'; 'Distance Estimate (m)'},'fontweight','bold','FontSize',48)
title('11 MHz - 13 MHz - 20 MHz')

set(gca,'FontSize',48); 
set(gca,'linew',4)

%% Figure zoomed into linear region 

figure;
plot(d_true,estimation_standard_deviation*100,'LineWidth',4,'Color',[0 0 0.5])
xlabel('Target Distance (m)','fontweight','bold','FontSize',48)
ylabel({'Standard Deviation of'; 'Distance Estimate (cm)'},'fontweight','bold','FontSize',48)
title('11 MHz - 13 MHz - 20 MHz')

set(gca,'FontSize',48); 
set(gca,'linew',4)

xlim([0 32])
ylim([0 19]) 
set(gca,'linew',4)
set(gca,'Color',[203 226 230]/255)

