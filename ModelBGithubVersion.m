%%%%%%% MODEL B %%%%%%%

% This script is for model b, the compartmental model model comprised of 
% 23 ODEs describing the concentration of NP in each of seven organ
% compartments: lung heart kindey liver spleen gut other.
% This version also has an updated model architecture (contains pulmonary
% circuit and spleen and gut lead into the liver. All parameter values 
% correspond to values representative of a murine model.

clc;clear;close all;
clearvars;
hold all


organs = {'lung'; 'heart'; 'kidney'; 'liver'; 'spleen'; 'gut'; 'others'};

% Define the inital concentration that enters in the vein
y0 = [0.00000000039 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %inital conditions

global  A_i K_EC l_NC_b V_i_BL V_i_T V_vein V_hep;
global  K_i_on K_i_off K_i_deg K_i_up K_NS Q_i Q_vein Q_hep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Physiological Parameter Values from literature   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_i     = [.0019956689 .0009621096822 .0033612343 .0045212562 .0004750356 .003333333 .003333333]*1000; %dm^2, crossectional area of cell free layer above EC 

l_NC_b  = (100*10^(-9)); % m, length of nanocarrier bond

Q_i     =  [1200 138 763.9168803 284.6716878 350.48828292 333 1000]/(60*10^6); %L/sec, flow rate through each organ compartment
      
V_i_BL  = [43.65525742 20.04395171  1018.55584 127.7673516 2.226207538 90  700]/10^6; %L, volume of blood in each organ
      
V_i_T   = [33.26114851 16.03516137 56.02057122 75.3542703 7.91726 20 1000]/10^6; %L, tissue volume of each organ compartment
  
V_vein  = 466.9/10^6; % L, volume of blood in the veins (and also the arteries)
      
Q_hep = sum(Q_i(4:6)); %L/sec, flow out of the liver

V_hep = sum(V_i_BL(4:6)); %L/sec, total volume of blood passed through liver

Q_vein = Q_i(3)+Q_hep+Q_i(7); %L/sec, total flow through vein

%% %%%%%%%%%%%%%%%%%%%%%%%%          Rate Parameters          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_i_on  = [1534.3 1534.3 1534.3 1534.3 1534.3 1534.3 1534.3]; % Binding rate of NP to endothelial cell surface receptors

K_EC    = [1.23*10^42 10.6565662 41793.3924 2.28*10^12 7.74*10^21 3.33E15 3.33E15]; %From Ramakrishnan 2016

K_i_off = [104711.5975 104711.5975 104711.5975 104711.5975 104711.5975 104711.5975 104711.5975]; %kon/kec, the rate of NP unbinding from surface receptors

%%%%%%%%%%%%%%%%%    Rate parameters determined by sensitivity analysis    %%%%%%%%%%%%%%%%%

K_i_deg = [1/50000000 1/500000 1/500 1/200000 1/1000000 1/300000 1/3000 1/600000]; % Rate of NP degradation in each organ compartment
               
      
K_i_up  = [100 1/40 1/200 100 3 1/1000 1/1000]; %Rate of NP uptake via receptor mediated endocytosis

    
K_NS    = [100 1/40 1/700 100 2.25 1/50000 1/50000]; %Rate of nonspecific uptake via transyctosis

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Experimental Validation data from Dong 2019  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%defining the discrete time points
time_dong = [0 5 30 60 120]/60; % time in hours

%extracted biodistribution data from Dong for 4 nm NP
kidney_dong_4 = [0 1.01062037 1.58456911 1.47288535 1.91962037];
liver_dong_4 = [0 15.8596502 18.430215 21.17352269 22.466739];
spleen_dong_4 = [0 1.87607819 2.45123554 2.41433904 3.17581354];
%error associated with the biodistribution data for 4 nm NP
kidney_error_dong_4 = [0 0 0 0 0];
liver_error_dong_4 = [0 6.24280019 6.24280017 10.8439633 10.2825591];
spleen_error_dong_4 = [0 .95798404 1.22818467 1.17905729 1.40013052];

%extracted biodistribution data from Dong for 15 nm NP
liver_dong_15 = [0 12.1874148 9.98407359 9.98407359 14.3907561];
spleen_dong_15 = [0 1.16373108 1.32150949 1.28440914 1.45625307];
%error associated with the biodistribution data for 15 nm NP
liver_error_dong_15 = [0 0 0 0 0];
spleen_error_dong_15 = [0 0 0 0 0];

%extracted biodistribution data from Dong for 50 nm NP
kidney_dong_50 = [0 1.02615034 1.24951784 .057941532 0];
liver_dong_50 = [0 23.9385681 23.204121 33.8536037 53.3164512];
spleen_dong_50 =[0 1.65520879 3.58136929 4.62507145 6.17248221];
%error associated with the biodistribution data for 50 nm NP
kidney_error_dong_50 = [0 0 0 0 0];
liver_error_dong_50 = [0 14.6889416 8.07891788 11.3839297 11.751153];
spleen_error_dong_50 = [0 1.1299299 .073670695 1.00721336 1.2751591];

%extracted biodistribution data from Dong for 79 nm NP
kidney_dong_79 = [0 1.25882482 0 0 0];
liver_dong_79 = [0 45.0143095 51.8406162 78.7151613 75.3773773];
spleen_dong_79 = [0 3.66937842 3.77762034 5.23292392 4.99335442];
%error associated with the biodistribution data for 79 nm NP
kidney_error_dong_79 = [0 0 0 0 0];
liver_error_dong_79 = [0 4.44825857 9.26881861 8.89658649 6.67396066];
spleen_error_dong_79 = [0 .29455663 .2948624 .46670408 .39311589];

%extracted biodistribution data from Dong for 100 nm NP
kidney_dong_100 = [0 .1326803 0 0 0];
liver_dong_100 = [0 64.3331574 80.1237696 83.0615579 84.1632286];
spleen_dong_100 = [0 3.4820602 4.21158546 4.45035765 4.34021763];
%error associated with the biodistribution data for 100 nm NP
kidney_error_dong_100 = [0 1.12145142 0 0 0];
liver_error_dong_100 = [0 .293778832 .293778832 .403945893 .477390602];
spleen_error_dong_100 = [0 .39476755 .80183288 .69773748 .80787751];

%% %%%%%%%%%%%%%%%%%%%%%%%%%       Running the ODE simulation       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_step = .001; %time, seconds

[t,y] = ode15s(@odefun,(0:t_step:10000),y0); %un ode simulation for 10,000 seconds

conc=0.00000000039; %define inital NP concentration

%%%%%%%%%%%%%%%%%%  Sum total number of NP in the system to ensure mass conservation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

%calculating the total NP in the vascular compartment of each organ
M_total_bl = 0;
M_total_bl_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_bl = M_total_bl+ ...
                 y(i,1)*V_vein*t_step*K_i_deg(8)    +... % vein 
                 y(i,2)*V_vein*t_step*K_i_deg(8)    +... % artery
                 y(i,3)*V_i_BL(1)*t_step*K_i_deg(1) +...
                 y(i,4)*V_i_BL(2)*t_step*K_i_deg(2) +...
                 y(i,5)*V_i_BL(3)*t_step*K_i_deg(3) +...
                 y(i,6)*V_i_BL(4)*t_step*K_i_deg(4) +...
                 y(i,7)*V_i_BL(5)*t_step*K_i_deg(5) +...
             	 y(i,18)*V_i_BL(6)*t_step*K_i_deg(6) +... % gut
                 y(i,21)*V_i_BL(7)*t_step*K_i_deg(7); % others   
                 
    M_total_bl_vec(i) = M_total_bl;
end

%calculating the total NP in the endothelial cell compartment of each organ
M_total_star = 0;
M_total_star_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_star = M_total_star+ ...
                   y(i,8)*A_i(1)*l_NC_b*t_step*K_i_deg(1) +...
                   y(i,9)*A_i(2)*l_NC_b*t_step*K_i_deg(2) +...
                   y(i,10)*A_i(3)*l_NC_b*t_step*K_i_deg(3) +...
                   y(i,11)*A_i(4)*l_NC_b*t_step*K_i_deg(4) +...
                   y(i,12)*A_i(5)*l_NC_b*t_step*K_i_deg(5) +...
                   y(i,19)*A_i(6)*l_NC_b*t_step*K_i_deg(6) +... % gut
                   y(i,22)*A_i(7)*l_NC_b*t_step*K_i_deg(7); % others
    M_total_star_vec(i) = M_total_star;
end

% Calculating the total NP in the tissue compartment of each organ
M_total_T = 0;
M_total_T_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_T = M_total_T +...
                y(i,13)*V_i_T(1)*t_step*K_i_deg(1) +...
                y(i,14)*V_i_T(2)*t_step*K_i_deg(2) +...
                y(i,15)*V_i_T(3)*t_step*K_i_deg(3) +...
                y(i,16)*V_i_T(4)*t_step*K_i_deg(4) +...
                y(i,17)*V_i_T(5)*t_step*K_i_deg(5) +...
                y(i,20)*V_i_T(6)*t_step*K_i_deg(6) +...  % gut
                y(i,23)*V_i_T(7)*t_step*K_i_deg(7);      % others
    M_total_T_vec(i) = M_total_T;
end

%conservation plot
hold on
subplot(5,2,[9,10])
plot(t/3600,((y(:,1)*V_vein)+(y(:,2)*V_vein)+...
       (y(:,3)*V_i_BL(1) + y(:,4)*V_i_BL(2) + y(:,5)*V_i_BL(3) + y(:,6)*V_i_BL(4) + y(:,7)*V_i_BL(5) + y(:,18)*V_i_BL(6) + y(:,21)*V_i_BL(7))+...
       (y(:,8)*A_i(1)*l_NC_b + y(:,9)*A_i(2)*l_NC_b + y(:,10)*A_i(3)*l_NC_b + y(:,11)*A_i(4)*l_NC_b + y(:,12)*A_i(5)*l_NC_b + y(:,19)*A_i(6)*l_NC_b + y(:,22)*A_i(7)*l_NC_b)+...
       (y(:,13)*V_i_T(1) + y(:,14)*V_i_T(2) + y(:,15)*V_i_T(3) + y(:,16)*V_i_T(4) + y(:,17)*V_i_T(5) + y(:,20)*V_i_T(6)+ y(:,23)*V_i_T(7))+...
       (M_total_T_vec)+...
       (M_total_star_vec)+...
       (M_total_bl_vec)))
ylim([0 5*10^-13])
title('Total Mass Conservation')
xlabel('Time (Hours)')
ylabel('Moles NP')

%% %%%%%%%%%%%%%%%%%%%%%    Plot biodistribuiton of NP in each organ compartment    %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dong Plot
hold on

%plot NP concentration in artery/vein
subplot(5,2,1)
semilogx((t/3600),y(:,1)/conc)%vein
hold on
semilogx((t/3600),y(:,2)/conc)%artery
lgd=legend('C_{vein}','C_{artery}','Location',"northeastoutside");
title(lgd, 'Arteries/Veins')
title('Normalized Concentration of Nanoparticles')
xlim([10^-2 3])

%Plot NP concentration in Lung
subplot(5,2,2)
semilogx(t/3600,y(:,13)/conc)
xlim([10^-2 3])
hold on
plot(2, 7.5230,'^')%This is the singular data point provided by Dong for the lung and 100nm NP
xlim([10^-2 3])
lgd=legend('C_{Lungs}','Location',"northwestoutside");
title(lgd, 'Lungs')

%Plot NP concentration in Heart
subplot(5,2,3)
semilogx(t/3600,y(:,14)/conc)
xlim([10^-2 3])
hold on
plot(2, .3389,'^')%This is the singular data point provided by Dong for the heart and 100nm NP
xlim([10^-2 3])
lgd=legend('C_{Heart}','Location',"northeastoutside");
title(lgd, 'Heart')

%Plot NP concentration in Kidney
subplot(5,2,4)
semilogx(t/3600,y(:,15)/conc)
xlim([10^-2 3])
hold on 
errorbar(time_dong, kidney_dong_100, kidney_error_dong_100, '^')
xlim([10^-2 3])
lgd=legend('C_{Kidneys}','Location',"northwestoutside");
title(lgd, 'Kidneys')

%Plot NP concentration in Liver
subplot(5,2,5)
semilogx(t/3600,((y(:,16) +(M_total_T_vec+M_total_star_vec+M_total_bl_vec)/V_i_T(4))/conc)) % added moles degraded
xlim([10^-2 3])
hold on
errorbar(time_dong, liver_dong_100, liver_error_dong_100, 'o')
xlim([10^-2 3])
ylabel('C_{NP}/C_{Tot NP}')
lgd=legend('C_{Liver}','Location',"northeastoutside");
title(lgd, 'Liver')

%Plot NP concentration in Spleen
subplot(5,2,6)
semilogx(t/3600,y(:,17)/conc)
xlim([10^-2 3])
hold on
errorbar(time_dong, spleen_dong_100, spleen_error_dong_100,'^')
xlim([10^-2 3])
lgd=legend('C_{Spleen}','Location',"northwestoutside");
title(lgd, 'Spleen')

%Plot NP concentration in Gut
subplot(5,2,7)
semilogx(t/3600,y(:,20)/conc)
xlim([10^-2 3])
hold on
lgd=legend('C_{Gut}','Location',"northeastoutside");
title(lgd, 'Gut')

%Plot NP concentration in Other
subplot(5,2,8)
semilogx(t/3600,y(:,23)/conc)
xlim([10^-2 3])
hold on
xlabel('Time (Hours)')
lgd=legend('C_{Other}','Location',"northwestoutside");
title(lgd, 'Other')

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define the system of 23 linear ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = odefun(t,y);

global  A_i K_EC l_NC_b V_i_BL V_i_T V_vein V_hep;
global  K_i_on K_i_off K_i_deg K_i_up K_NS Q_i Q_vein Q_hep;

Q_hep = sum(Q_i(4:6)); 
V_hep = sum(V_i_BL(4:6));
Q_vein = Q_i(3)+Q_hep+Q_i(7); 

dydt = zeros(23,1);

% Vein 
dydt(1) = (Q_i(3)*y(5) + Q_hep*y(6) + Q_i(7)*y(21) - Q_vein*y(1) - y(1)*V_vein*K_i_deg(8)) /V_vein; 

% Artery 
dydt(2) = (Q_vein*y(4)-(Q_i(3)*y(2) + Q_i(4)*y(2)+ Q_i(5)*y(2) + Q_i(6)*y(2) + Q_i(7)*y(2) - y(2)*V_vein*K_i_deg(8)))/ V_vein;  

%%%%%%%%%%%%%%%%%%%%% Vascular compartment (organs i=1:5) %%%%%%%%%%%%%%%%%%%%%%%%%

% Lung 
dydt(3) = (Q_i(1)*y(4) - Q_i(1)*y(3)-K_i_on(1)*y(3)*V_i_BL(1)+K_i_off(1)*y(8)*A_i(1)*l_NC_b-V_i_BL(1)*K_NS(1)*y(3)-y(3)*V_i_BL(1)*K_i_deg(1))/(V_i_BL(1));

% Heart
dydt(4) = (Q_vein*y(1)+Q_i(1)*y(3)- Q_i(1)*y(4) - Q_vein*y(4)-K_i_on(2)*y(4)*V_i_BL(2)+K_i_off(2)*y(9)*A_i(2)*l_NC_b-V_i_BL(2)*K_NS(2)*y(4)-y(4)*V_i_BL(2)*K_i_deg(2))/(V_i_BL(2));

% Kidney
dydt(5) = (Q_i(3)*y(2)-Q_i(3)*y(5)-K_i_on(3)*y(5)*V_i_BL(3)+K_i_off(3)*y(10)*A_i(3)*l_NC_b-V_i_BL(3)*K_NS(3)*y(5)-y(5)*V_i_BL(3)*K_i_deg(3))/(V_i_BL(3));

% Liver 
dydt(6) = ((Q_i(5)*y(7) + Q_i(6)*y(18)+Q_i(4)*y(2)-Q_hep*y(6))-K_i_on(4)*y(6)*V_i_BL(4)+K_i_off(4)*y(11)*A_i(4)*l_NC_b-V_i_BL(4)*K_NS(4)*y(6)-y(6)*V_i_BL(4)*K_i_deg(4))/(V_i_BL(4));

% Spleen
dydt(7) = (Q_i(5)*y(2)-Q_i(5)*y(7)-K_i_on(5)*y(7)*V_i_BL(5)+K_i_off(5)*y(12)*A_i(5)*l_NC_b-V_i_BL(5)*K_NS(5)*y(7)-y(7)*V_i_BL(5)*K_i_deg(5))/(V_i_BL(5));

%%%%%%%%%%%%%%%%%%%%% Endothelial Cell compartment (organs i=1:5) %%%%%%%%%%%%%%%%%%%%%%%%%

% Lung
dydt(8) = (-A_i(1)*l_NC_b*K_i_up(1)*y(8)+K_i_on(1)*V_i_BL(1)*y(3)-K_i_off(1)*y(8)*A_i(1)*l_NC_b-y(8)*A_i(1)*l_NC_b*K_i_deg(1))/(A_i(1)*l_NC_b);

% Heart
dydt(9) = (-A_i(2)*l_NC_b*K_i_up(2)*y(9)+K_i_on(2)*V_i_BL(2)*y(4)-K_i_off(2)*y(9)*A_i(2)*l_NC_b-y(9)*A_i(2)*l_NC_b*K_i_deg(2))/(A_i(2)*l_NC_b);

% Kidney
dydt(10) = (-A_i(3)*l_NC_b*K_i_up(3)*y(10)+K_i_on(3)*V_i_BL(3)*y(5)-K_i_off(3)*y(10)*A_i(3)*l_NC_b-y(10)*A_i(3)*l_NC_b*K_i_deg(3))/(A_i(3)*l_NC_b);

% Liver 
dydt(11) = (-A_i(4)*l_NC_b*K_i_up(4)*y(11)+K_i_on(4)*V_i_BL(4)*y(6)-K_i_off(4)*y(11)*A_i(4)*l_NC_b-y(11)*A_i(4)*l_NC_b*K_i_deg(4))/(A_i(4)*l_NC_b);

% Spleen
dydt(12) = (-A_i(5)*l_NC_b*K_i_up(5)*y(12)+K_i_on(5)*V_i_BL(5)*y(7)-K_i_off(5)*y(12)*A_i(5)*l_NC_b-y(12)*A_i(5)*l_NC_b*K_i_deg(5))/(A_i(5)*l_NC_b);

%%%%%%%%%%%%%%%%%%%% Tissue Compartment (organs i=1:5)  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lung
dydt(13) = (A_i(1)*l_NC_b*K_i_up(1)*y(8)+V_i_BL(1)*K_NS(1)*y(3)-V_i_T(1)*K_i_deg(1)*y(13))/(V_i_T(1));

% Heart
dydt(14) = (A_i(2)*l_NC_b*K_i_up(2)*y(9)+V_i_BL(2)*K_NS(2)*y(4)-V_i_T(2)*K_i_deg(2)*y(14))/(V_i_T(2));

% Kidney
dydt(15) = (A_i(3)*l_NC_b*K_i_up(3)*y(10)+V_i_BL(3)*K_NS(3)*y(5)-V_i_T(3)*K_i_deg(3)*y(15))/(V_i_T(3));

% Liver
dydt(16) = (A_i(4)*l_NC_b*K_i_up(4)*y(11)+V_i_BL(4)*K_NS(4)*y(6)-V_i_T(4)*K_i_deg(4)*y(16))/(V_i_T(4));

% Spleen
dydt(17) = (A_i(5)*l_NC_b*K_i_up(5)*y(12)+V_i_BL(5)*K_NS(5)*y(7)-V_i_T(5)*K_i_deg(5)*y(17))/(V_i_T(5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gut Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gut Vasculature
dydt(18) = (Q_i(6)*y(2)-Q_i(6)*y(18)-K_i_on(6)*y(18)*V_i_BL(6)+K_i_off(6)*y(19)*A_i(6)*l_NC_b-V_i_BL(6)*K_NS(7)*y(18)-y(18)*V_i_BL(6)*K_i_deg(6))/(V_i_BL(6));

% Gut Endothelial Cell
dydt(19) = (-A_i(6)*l_NC_b*K_i_up(6)*y(19)+K_i_on(6)*V_i_BL(6)*y(18)-K_i_off(6)*y(19)*A_i(6)*l_NC_b-y(19)*A_i(6)*l_NC_b*K_i_deg(6))/(A_i(6)*l_NC_b);

% Gut Tissue
dydt(20) = (A_i(6)*l_NC_b*K_i_up(6)*y(19)+V_i_BL(6)*K_NS(6)*y(18)-V_i_T(6)*K_i_deg(6)*y(20))/(V_i_T(6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Other Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Other Vasculature
dydt(21) = (Q_i(7)*y(2)-Q_i(7)*y(21)-K_i_on(7)*y(21)*V_i_BL(7)+K_i_off(7)*y(22)*A_i(7)*l_NC_b-V_i_BL(7)*K_NS(7)*y(21)-y(21)*V_i_BL(7)*K_i_deg(7))/(V_i_BL(7));

% Other Endothelial Cell
dydt(22) = (-A_i(7)*l_NC_b*K_i_up(7)*y(22)+K_i_on(7)*V_i_BL(7)*y(21)-K_i_off(7)*y(22)*A_i(7)*l_NC_b-y(22)*A_i(7)*l_NC_b*K_i_deg(7))/(A_i(7)*l_NC_b);

% Other Tissue
dydt(23) = (A_i(7)*l_NC_b*K_i_up(7)*y(22)+V_i_BL(7)*K_NS(7)*y(21)-V_i_T(7)*K_i_deg(7)*y(23))/(V_i_T(7));

end
