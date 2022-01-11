%%
n = 0:0.001:1;
figure
plot(n,besselj(1,n),'Color',[0 0 0.5])
hold on
plot(n,-0.5*ones(size(n)),'Color',[0 0 0])
n_temp = 0.1436:0.001:0.2245;
x2 = [n_temp, fliplr(n_temp)];
inBetween = [besselj(1,n_temp) - 0.005, fliplr((-0.5 + 0.005)*ones(size(n_temp)))];
fill(x2, inBetween, [203 226 230]/255,'LineStyle','none');
xlabel('\phi_D (radians)','FontSize',60)
ylabel('J_1(\phi_D)','FontSize',60)
set(gca,'FontSize',60)
set(findall(gca,'Type','Line'),'LineWidth',4);
set(gca,'linewidth',4);
ylim([-0.5 0.7])
yticks([-0.5 -0.25 0 0.25 0.5 1])


%Plotting Static Phase Variation, Dynamic Phase Variation with Angle 
lambda = 532e-9;
L = 500e-6*3.7696e6/3.73386e6;
S_zy = 2e-5*3;
p_14 = 1/20;

n_o = 2.3 + (147-500)*0.5/50000;
n_e = 2.21 + (1000-500)*0.5/50000;

theta_vector = pi/2:-0.001:pi/2-pi/3;
psi_vector = pi/2:-0.001:pi/2-pi/3;
c_o = zeros(length(pi/2:-0.001:pi/2-pi/3),length(pi/2:-0.001:pi/2-pi/3));
c_e = zeros(length(pi/2:-0.001:pi/2-pi/3),length(pi/2:-0.001:pi/2-pi/3));
delta_phi_s = zeros(length(pi/2:-0.001:pi/2-pi/3),length(pi/2:-0.001:pi/2-pi/3));
delta_phi_d = zeros(length(pi/2:-0.001:pi/2-pi/3),length(pi/2:-0.001:pi/2-pi/3));

for u=1:size(c_o,1)
    for w=1:size(c_o,2)
        theta = theta_vector(u);
        psi = psi_vector(w);
        phi = acos(sin(theta).*sin(psi));
        
        theta_o_tilde = asin(sin(phi)/n_o);
        theta_e_tilde = asin(n_o*sin(phi)*sqrt((tan(theta))^2*(cos(psi))^2 + 1)/sqrt(n_o^2*n_e^2*((tan(theta))^2*(cos(psi))^2 + 1) - n_e^2*(sin(phi))^2 + n_o^2*(sin(phi))^2));
        theta_e = acos(n_o*sin(phi)/sqrt(n_o^2*n_e^2*((tan(theta))^2*(cos(psi))^2 + 1) - n_e^2*(sin(phi))^2 + n_o^2*(sin(phi))^2));
        
        n_eff = ((cos(theta_e))^2/n_o^2 + (sin(theta_e))^2/n_e^2)^-0.5;
        
        delta_phi_s(u,w) = 2*pi*L/lambda*(n_o*cos(theta_o_tilde) - n_eff*cos(theta_e_tilde));
        
        n_o_dynamic = -n_o*p_14*S_zy*n_o^2*(n_o^2*(cos(asin(sin(phi)/n_o)))^2*((tan(theta))^2*(cos(psi))^2 + 1) - (sin(phi))^2*(tan(theta))^2*(cos(psi))^2)/(n_o^2*((tan(theta))^2*(cos(psi))^2 + 1)*(cos(asin(sin(phi)/n_o)))^2 + (sin(phi))^2*(tan(theta))^2*(cos(psi))^2);
        n_e_dynamic = n_o*n_e/sqrt(n_o^2*(sin(theta_e))^2 + n_e^2*(cos(theta_e))^2)*-1*p_14*S_zy*n_o^2*(n_e^2*(cos(theta_e))^4*(tan(theta))^2*(cos(psi))^2 - (cos(theta_e_tilde))^2*(cos(theta_e))^2*n_e^2)/((n_o^2*(sin(theta_e))^2 + n_e^2*(cos(theta_e))^2)*((cos(theta_e))^4*(tan(theta))^2*(cos(psi))^2 + (cos(theta_e_tilde))^2*(cos(theta_e))^2 + ((cos(theta_e_tilde))^2 + ((cos(theta_e))^2*(tan(theta))^2*(cos(psi))^2)^2)));
        
        delta_phi_d(u,w) = 2*pi*L/lambda*(n_o_dynamic*cos(theta_o_tilde) - n_e_dynamic*cos(theta_e_tilde));
        
        denom = n_o.^2.*n_e.^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1) - n_e.^2.*(sin(phi)).^2 + n_o.^2.*(sin(phi)).^2;
        tan_term = (tan(theta)).^2.*(cos(psi)).^2 + 1;
    
        %% Expressing a_1, a_2, a_3 in terms of theta,psi
        %Solving ax^2 + bx + c = 0
        a = (tan(theta).*cos(psi)).^2 + ((cos(theta) + sin(theta).*tan(theta).*(cos(psi)).^2)./(sin(theta).*sin(psi))).^2 + 1;
        b = 2.*tan(theta).*cos(psi).*cos(asin(sin(phi)/n_o))./sqrt(((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)))) + 2.*cot(psi).*(cos(asin(sin(phi)/n_o))).*(cos(theta) + sin(theta).*tan(theta).*(cos(psi)).^2)./(sin(theta).*sin(psi).*sqrt(((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)))));
        c = (cos(asin(sin(phi)/n_o))).^2./((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1))) + (cot(psi)).^2.*(cos(asin(sin(phi)/n_o))).^2./((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1))) - 1;
        a_3_1 = (-b + sqrt(b.^2 - 4.*a.*c))./(2.*a);
        a_3_2 = (-b - sqrt(b.^2 - 4.*a.*c))./(2.*a);

        a_1_1 = a_3_1.*tan(theta).*cos(psi) + cos(asin(sin(phi)/n_o))./sqrt((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)));
        a_2_1 = -1*a_3_1.*(cos(theta) + sin(theta).*tan(theta).*(cos(psi)).^2)./(sin(theta).*sin(psi)) - cot(psi).*cos(asin(sin(phi)/n_o))./sqrt((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)));

        a_1_2 = a_3_2.*tan(theta).*cos(psi) + cos(asin(sin(phi)/n_o))./sqrt((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)));
        a_2_2 = -1*a_3_2.*(cos(theta) + sin(theta).*tan(theta).*(cos(psi)).^2)./(sin(theta).*sin(psi)) - cot(psi).*cos(asin(sin(phi)/n_o))./sqrt((cos(asin(sin(phi)/n_o))).^2 + (sin(phi)).^2.*(tan(theta)).^2.*(cos(psi)).^2./(n_o^2.*((tan(theta)).^2.*(cos(psi)).^2 + 1)));

        %Finding which solution is the correct one
        dot_product_1 = a_1_1.*cos(asin(sin(phi)/n_o)) - a_2_1.*sin(phi).*tan(theta).*cos(psi)./(n_o.*sqrt(tan(theta).^2.*(cos(psi)).^2 + 1));
        dot_product_2 = a_1_2.*cos(asin(sin(phi)/n_o)) - a_2_2.*sin(phi).*tan(theta).*cos(psi)./(n_o.*sqrt(tan(theta).^2.*(cos(psi)).^2 + 1));

        if dot_product_1 > dot_product_2 
            a_1 = a_1_1;
            a_2 = a_2_1;
            a_3 = a_3_1;
        else
            a_1 = a_1_2;
            a_2 = a_2_2;
            a_3 = a_3_2;
        end

        %% Expressing b_1, b_2, b_3 in terms of a_1, a_2, a_3 
        %% Expressing p_ei components in terms of p_oi components
        denominator = sqrt((a_2.*cos(theta) - a_3.*sin(theta).*sin(psi)).^2 +(a_3.*sin(theta).*cos(psi) - a_1.*cos(theta)).^2 + (a_1.*sin(theta).*sin(psi) - a_2.*sin(theta).*cos(psi)).^2);
        b_1 = (a_2.*cos(theta) - a_3.*sin(theta).*sin(psi))./denominator;
        b_2 = (-1*a_1.*cos(theta) + a_3.*sin(theta).*cos(psi))./denominator;
        b_3 = (a_1.*sin(theta).*sin(psi) - a_2.*sin(theta).*cos(psi))./denominator;

        %% Expressions for c_o and c_e
        denom_2 = sqrt((sin(theta).*cos(theta).*cos(psi) - (cos(theta)).^2 - (sin(theta)).^2.*(sin(psi)).^2).^2 + ((sin(theta)).^2.*cos(psi).*sin(psi) + sin(theta).*cos(theta).*sin(psi)).^2 + (sin(theta).*cos(theta).*cos(psi) - (sin(theta)).^2).^2);
        c_o(u,w) = (a_1.*(sin(theta).*cos(theta).*cos(psi) - (cos(theta)).^2 - (sin(theta)).^2.*(sin(psi)).^2) + a_2.*((sin(theta)).^2.*cos(psi).*sin(psi) + sin(theta).*cos(theta).*sin(psi)) + a_3.*(sin(theta).*cos(theta).*cos(psi) - (sin(theta)).^2))./denom_2;
        c_e(u,w) = (b_1.*(sin(theta).*cos(theta).*cos(psi) - (cos(theta)).^2 - (sin(theta)).^2.*(sin(psi)).^2) + b_2.*((sin(theta)).^2.*cos(psi).*sin(psi) + sin(theta).*cos(theta).*sin(psi)) + b_3.*(sin(theta).*cos(theta).*cos(psi) - (sin(theta)).^2))./denom_2;
    
    end
end

figure;
imagesc(theta_vector*180/pi,psi_vector*180/pi,abs(sin(delta_phi_s)))
hold on
xlabel('\theta({\circ})','FontSize',72)
ylabel('\psi({\circ})','FontSize',72)
title('|sin(\phi_{s})|','fontweight','bold','FontSize',48)
set(gca,'FontSize',72)
set(gca,'linewidth',3);
xticks([50 70 90])
yticks([50 70 90])
shading interp 

figure;
imagesc(theta_vector*180/pi,psi_vector*180/pi,-1*delta_phi_d)
xlabel('\theta({\circ})','FontSize',72)
ylabel('\psi({\circ})','FontSize',72)
set(gca,'FontSize',72)
set(gca,'linewidth',3);
xticks([50 70 90])
yticks([50 70 90])
shading interp 

figure;
imagesc(theta_vector*180/pi,psi_vector*180/pi,c_o.^2.*c_e.^2)
xlabel('\theta({\circ})','FontSize',72)
ylabel('\psi({\circ})','FontSize',72)
set(gca,'FontSize',72)
set(gca,'linewidth',3);
xticks([50 70 90])
yticks([50 70 90])
shading interp 

figure;
imagesc(theta_vector*180/pi,psi_vector*180/pi,100*abs(4*c_o.^2.*c_e.^2.*sin(delta_phi_s).*besselj(1,delta_phi_d))./(c_o.^4 + c_e.^4 + 2*c_o.^2.*c_e.^2.*cos(delta_phi_s).*besselj(0,delta_phi_d)))
xlabel('\theta({\circ})','FontSize',72)
ylabel('\psi({\circ})','FontSize',72)
set(gca,'FontSize',72)
set(gca,'linewidth',3);
xticks([50 70 90])
yticks([50 70 90])
%shading interp 

figure;
imagesc(theta_vector*180/pi,psi_vector*180/pi,100*abs(4*c_o.^2.*c_e.^2.*besselj(1,delta_phi_d))./(c_o.^4 + c_e.^4 + 2*c_o.^2.*c_e.^2.*0.*besselj(0,delta_phi_d)))
xlabel('\theta({\circ})','FontSize',72)
ylabel('\psi({\circ})','FontSize',72)
title({'Depth of Modulation', 'with Metasurface'},'FontSize',72)
set(gca,'FontSize',72)
set(gca,'linewidth',3);
xticks([50 70 90])
yticks([50 70 90])
shading interp 


