%%%%%%% MODEL A %%%%%%%

% This script is for model a, the compartmental model model comprised of 
% 17 ODEs describing the concentration of NP in each of five organ
% compartments: lung heart kindey liver spleen. All parameter values 
% correspond to values representative of a murine model.

clc;clear;close all;
clearvars;


hold all
organs = {'lung';'heart';'kidney';'liver';'spleen'};

% Define the intial concentration that enters in the vein
y0 = [0.00000000039 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %initial conditions

global  A_i K_EC l_NC_b V_i_BL V_i_T V_vein;
global  K_i_on K_i_off K_i_deg K_i_up K_NS Q_i Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Variable Values from literature   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_i     = [.0019956689 .0009621096822 .0033612343 .0045212562 .0004750356]*1000; % dm^2, 

l_NC_b  = (10^(-6)); %dm, length of the nanocarrier bond
          
K_EC    = [1.23*10^42 10.6565662 41793.3924 2.28*10^12 7.74*10^21]; %from ramakrishnan 2016

K_i_on  = [1534.3 1534.3 1534.3 1534.3 1534.3]; %Binding rate of NP to endothelial cell surface receptors

K_i_off = [(K_i_on(1)/(log(K_EC(1)))) (K_i_on(2)/(log(K_EC(2)))) (K_i_on(3)/(log(K_EC(3)))) (K_i_on(4)/(log(K_EC(4)))) (K_i_on(5)/(log(K_EC(5))))]; % k_on/k_ec, the rate of NP unbinding from surface receptors
         
Q_i     = [1100.83758415 138.673326 763.9168803 284.6716878 224.48828292]/(60*10^6); %L/sec, blood flow rate in each organ
      
V_i_BL  = [43.65525742 20.04395171  1018.55584 127.7673516 2.226207538]/10^6; %L, volume of blood in each organ
          
V_i_T   = [33.26114851 16.03516137 56.02057122 75.3542703 7.91726]/10^6; %L, tissue volume of each organ compartment
          
V_vein  = 466.9/10^6; % L, volume of blood in the veins (and also the arteries
          
Q       =  Q_i(1)+Q_i(2)+Q_i(3)+Q_i(4)+Q_i(5); %L/sec total flow rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Variable Values From Sensitivity Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

K_i_deg = [1/500000     1/500000   1/500     1/20000   1/200000 ]; %rate of NP degredation in each organ compartment
          
        
K_i_up  = [1000            1/170      1/1000      100      1/18 ]; %rate of NP update via receptor mediated endocytosis
          
    
K_NS    = [1000            1/100       1/1000      100       1/20    ]; %rate of nonspecific uptake via transyctosis
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Experimental Validation data from Dong 2019   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%defining the discrete time points 
time_dong = [0 5 30 60 120]/60; %time in hours

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Running ODE simulation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_step = 0.001; %time, seconds

[t,y] = ode15s(@odefun,(0:t_step:10000),y0); %run ode simulation for 10,000 seconds

conc=0.00000000039; %define initial NP concentration

%%%%%%%%%%%%%%%%%%  Sum total number of NP in the system to ensure mass conservation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on

%calculating the total NP in the vascular compartment of each organ
M_total_bl = 0;
M_total_bl_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_bl = M_total_bl+ ...
                 y(i,3)*V_i_BL(1)*t_step*K_i_deg(1) +...
                 y(i,4)*V_i_BL(2)*t_step*K_i_deg(2) +...
                 y(i,5)*V_i_BL(3)*t_step*K_i_deg(3) +...
                 y(i,6)*V_i_BL(4)*t_step*K_i_deg(4) +...
                 y(i,7)*V_i_BL(5)*t_step*K_i_deg(5);
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
                   y(i,12)*A_i(5)*l_NC_b*t_step*K_i_deg(5);
    M_total_star_vec(i) = M_total_star;
end


%%calculating the total NP in the tissue compartment of each organ
M_total_T = 0;
M_total_T_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_T = M_total_T +...
                y(i,13)*V_i_T(1)*t_step*K_i_deg(1) +...
                y(i,14)*V_i_T(2)*t_step*K_i_deg(2) +...
                y(i,15)*V_i_T(3)*t_step*K_i_deg(3) +...
                y(i,16)*V_i_T(4)*t_step*K_i_deg(4) +...
                y(i,17)*V_i_T(5)*t_step*K_i_deg(5);
    M_total_T_vec(i) = M_total_T;
end

%conservation plot
subplot(4,2,[7,8])
plot(t/3600,((y(:,1)*V_vein)+...
       (y(:,2)*V_vein)+...
       (y(:,3)*V_i_BL(1) + y(:,4)*V_i_BL(2) + y(:,5)*V_i_BL(3) + y(:,6)*V_i_BL(4) + y(:,7)*V_i_BL(5))+...
       (y(:,8)*A_i(1)*l_NC_b + y(:,9)*A_i(2)*l_NC_b + y(:,10)*A_i(3)*l_NC_b + y(:,11)*A_i(4)*l_NC_b + y(:,12)*A_i(5)*l_NC_b)+...
       (y(:,13)*V_i_T(1) + y(:,14)*V_i_T(2) + y(:,15)*V_i_T(3) + y(:,16)*V_i_T(4) + y(:,17)*V_i_T(5))+...
       (M_total_T_vec)+...
       (M_total_star_vec)+...
       (M_total_bl_vec)))
   
ylim([1*10e-14 3*10e-14])
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot biodistribution of NP in each organ compartment %%%%%%%%%%%%%%%%%%%%%%%%%%


%Plot NP concentration in artery/vein
subplot(4,2,1)
semilogx(t/3600,y(:,1)/conc,'^') %vein
hold on
semilogx(t/3600,y(:,2)/conc,'.') %artery
lgd=legend('C_{Vein}^{T}', 'C_{Artery}^{T}','Location',"northeast");
title(lgd, 'Arteries/Veins')

%Plot NP concentration in Lung
subplot(4,2,2)
semilogx(t/3600,y(:,13)/conc)
xlim([10^-2 3])
hold on
plot(2, 7.5230,'o') %This is the singular data point provided by Dong for the lung and 100nm NP
xlim([10^-2 3])
legend('C_{Lung}^{T}')

%Plot NP concentration in heart
subplot(4,2,3)
semilogx(t/3600,y(:,14)/conc)
xlim([10^-2 3])
hold on
plot(2, .3389,'o')  %This is the singular data point provided by Dong for the heart and 100nm NP
xlim([10^-2 3])
legend('C_{Heart}^{T}')

%Plot NP concentration in the kidneys
subplot(4,2,4)
semilogx(t/3600,y(:,15)/conc)
xlim([10^-2 3])
hold on 
errorbar(time_dong, kidney_dong_100, kidney_error_dong_100, '^')
xlim([10^-2 3])
legend('C_{Kidneys}^{T}')

%PLot the NP concentration in the liver
subplot(4,2,5)
semilogx(t/3600,((y(:,16) +(M_total_T_vec+M_total_star_vec+M_total_bl_vec)/V_i_T(4))/conc))
xlim([10^-2 3])
hold on
errorbar(time_dong, liver_dong_100, liver_error_dong_100, '^')
xlim([10^-2 3])
legend('C_{Liver}^{T}')

%Plot the NP concentration in the spleen
subplot(4,2,6)
semilogx(t/3600,y(:,17)/conc)
xlim([10^-2 3])
hold on
errorbar(time_dong, spleen_dong_100, spleen_error_dong_100,'^')
xlim([10^-2 3])
legend('C_{Spleen}^{T}')

legend('C_{Spleen}^{T}')
xlabel('Time (Hours)')
ylabel('C_{NP}/C_{NP Tot}')

ylabel('Moles NP')
xlabel('Time Hours')
   

%%%%%%%%%%%%%%%%%%% Define the system of 17 linear ODEs %%%%%%%%%%%%%%%%%   
function dydt = odefun(t,y)

global  A_i K_EC l_NC_b V_i_BL V_i_T V_vein;
global  K_i_on K_i_off K_i_deg K_i_up K_NS Q_i Q;         

dydt= zeros(17,1);

%artery
dydt(1) = ((Q_i(1)*y(3)+Q_i(2)*y(4)+Q_i(3)*y(5)+Q_i(4)*y(6)+Q_i(5)*y(7))-Q*y(1))/(V_vein);

%vein
dydt(2) = ((Q*y(1)-(Q_i(1)*y(2)+Q_i(2)*y(2)+Q_i(3)*y(2)+Q_i(4)*y(2)+Q_i(5)*y(2)))/(V_vein));

%%%%%%%%%%%%%%%%%%%%% Vascular compartment %%%%%%%%%%%%%%%%%%%%%%%%%

%lung 
dydt(3) = (Q_i(1)*y(2)-Q_i(1)*y(3)-K_i_on(1)*y(3)*V_i_BL(1)+K_i_off(1)*y(8)*A_i(1)*l_NC_b-V_i_BL(1)*K_NS(1)*y(3)-y(3)*V_i_BL(1)*K_i_deg(1))/(V_i_BL(1));

%heart 
dydt(4) = (Q_i(2)*y(2)-Q_i(2)*y(4)-K_i_on(2)*y(4)*V_i_BL(2)+K_i_off(2)*y(9)*A_i(2)*l_NC_b-V_i_BL(2)*K_NS(2)*y(4)-y(4)*V_i_BL(2)*K_i_deg(2))/(V_i_BL(2));

%kidney
dydt(5) = (Q_i(3)*y(2)-Q_i(3)*y(5)-K_i_on(3)*y(5)*V_i_BL(3)+K_i_off(3)*y(10)*A_i(3)*l_NC_b-V_i_BL(3)*K_NS(3)*y(5)-y(5)*V_i_BL(3)*K_i_deg(3))/(V_i_BL(3));

%liver
dydt(6) = (Q_i(4)*y(2)-Q_i(4)*y(6)-K_i_on(4)*y(6)*V_i_BL(4)+K_i_off(4)*y(11)*A_i(4)*l_NC_b-V_i_BL(4)*K_NS(4)*y(6)-y(6)*V_i_BL(4)*K_i_deg(4))/(V_i_BL(4));

%spleen
dydt(7) = (Q_i(5)*y(2)-Q_i(5)*y(7)-K_i_on(5)*y(7)*V_i_BL(5)+K_i_off(5)*y(12)*A_i(5)*l_NC_b-V_i_BL(5)*K_NS(5)*y(7)-y(7)*V_i_BL(5)*K_i_deg(5))/(V_i_BL(5));

%%%%%%%%%%%%%%%%%%%%% Endothelial Cell compartment %%%%%%%%%%%%%%%%%%%%%%%%%

%lung
dydt(8) = (-A_i(1)*l_NC_b*K_i_up(1)*y(8)+K_i_on(1)*V_i_BL(1)*y(3)-K_i_off(1)*y(8)*A_i(1)*l_NC_b-y(8)*A_i(1)*l_NC_b*K_i_deg(1))/(A_i(1)*l_NC_b);

%heart
dydt(9) = (-A_i(2)*l_NC_b*K_i_up(2)*y(9)+K_i_on(2)*V_i_BL(2)*y(4)-K_i_off(2)*y(9)*A_i(2)*l_NC_b-y(9)*A_i(2)*l_NC_b*K_i_deg(2))/(A_i(2)*l_NC_b);

%kidney
dydt(10) = (-A_i(3)*l_NC_b*K_i_up(3)*y(10)+K_i_on(3)*V_i_BL(3)*y(5)-K_i_off(3)*y(10)*A_i(3)*l_NC_b-y(10)*A_i(3)*l_NC_b*K_i_deg(3))/(A_i(3)*l_NC_b);

%liver
dydt(11) = (-A_i(4)*l_NC_b*K_i_up(4)*y(11)+K_i_on(4)*V_i_BL(4)*y(6)-K_i_off(4)*y(11)*A_i(4)*l_NC_b-y(11)*A_i(4)*l_NC_b*K_i_deg(4))/(A_i(4)*l_NC_b);

%spleen
dydt(12) = (-A_i(5)*l_NC_b*K_i_up(5)*y(12)+K_i_on(5)*V_i_BL(5)*y(7)-K_i_off(5)*y(12)*A_i(5)*l_NC_b-y(12)*A_i(5)*l_NC_b*K_i_deg(5))/(A_i(5)*l_NC_b);

%%%%%%%%%%%%%%%%%%%% Tissue Compartment %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lung
dydt(13) = (A_i(1)*l_NC_b*K_i_up(1)*y(8)+V_i_BL(1)*K_NS(1)*y(3)-V_i_T(1)*K_i_deg(1)*y(13))/(V_i_T(1));

%heart
dydt(14) = (A_i(2)*l_NC_b*K_i_up(2)*y(9)+V_i_BL(2)*K_NS(2)*y(4)-V_i_T(2)*K_i_deg(2)*y(14))/(V_i_T(2));

%kidney
dydt(15) = (A_i(3)*l_NC_b*K_i_up(3)*y(10)+V_i_BL(3)*K_NS(3)*y(5)-V_i_T(3)*K_i_deg(3)*y(15))/(V_i_T(3));

%liver
dydt(16) = (A_i(4)*l_NC_b*K_i_up(4)*y(11)+V_i_BL(4)*K_NS(4)*y(6)-V_i_T(4)*K_i_deg(4)*y(16))/(V_i_T(4));

%spleen
dydt(17) = (A_i(5)*l_NC_b*K_i_up(5)*y(12)+V_i_BL(5)*K_NS(5)*y(7)-V_i_T(5)*K_i_deg(5)*y(17))/(V_i_T(5));

end
