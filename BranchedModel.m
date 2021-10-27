%%%%%%% BRANCHED MODEL %%%%%%%

% This script is for the branched model, the model comprised of 
% 457 ODEs describing the concentration of NP in each of seven organ
% compartments: lung heart kindey liver spleen gut other, but using a
% branched architecture in in each of the organ compartments. All parameter 
% values correspond to values representative of a murine model.


clc;clear; close all;
clearvars;

hold all

% Define the inital concentration that enters the vein 
y0 = zeros( 1, 457);
y0(1) = 0.00000000039;

global num_gen num_branches_ls all_D all_SA all_V all_Q  V_T V_vein
global K_on K_off K_EC K_NS K_up K_deg K_deg_vein l_NC_b Q Kup_term Kns_term Q_vein V_i_BL

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Branching characterisitics for each organ used in all organs %%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 3; %vessel length to diameter ratio, constant for all generations

%---------------------------------------calculate diameter of each branch----------------------------------------------

D = 0.0006; %initial diameter of artery entering element (meters)
all_D = [];
num_gen = 0;
if D==0.0002
    all_D = [0.0002];
    num_gen = 1;
else 
    while D >= .000015            % meters
        D = ((D^3)/2)^(1/3);      %calculate diameter of branch in next generation given branch diameter in previous generation
        all_D = [all_D D];      
        num_gen = num_gen+1;
    end
end


if num_gen>=2                      
    for i = length(all_D):-1:1
        all_D = [all_D all_D(i)];  %backwards iteration adds values to existing matrix to account for other half of branching architecture
    end
else
    all_D = [0.0002];  % if we specify there to be only one generation we don't do backwards iteration
end


%------------------------------------calculate length of each branch---------------------------------------------------

L = x*all_D; %vector containing all the lengths of the blood vessels in each generation

%-------------------------Calculate the total surface area of all branches combined----------------------------------

C = all_D * 3.14; %circumfrence of one branch at each generation
SA = C .* L;      %surface area of one branch at every generation
all_SA = [];
num_branches_ls = [];

%calculating total surface area across the branching element (forward)
if num_gen >=2
    for i = 1:num_gen
        num_branches_ls = [num_branches_ls 2^i] ; 
        num_branches = 2^i;             %number of branches in each generation
        SA_tot = SA(i)*num_branches;    %total surface area of all branches in a generation
        all_SA = [all_SA SA_tot];       %add that surface area to the vector
    end
else
    all_SA = [SA];
    num_branches = [1];
    num_branches_ls = [1];
end
%calculating total surface area across the branching element (back)
if num_gen >=2
    for i = length(all_SA):-1:1          %total SA for each generation entire network
        all_SA = [all_SA all_SA(i)];
        num_branches_ls = [num_branches_ls num_branches_ls(i)];%total num branches each generation entire network
    end
else
    all_SA = all_SA;
    num_branches_ls = [1];
end

network_SA = sum(all_SA);                %surface area of vasculature branches (SA of EC also)

%-----------------------calculate total volume of all branches combined---------------------------------------------

A = ((all_D/2).^2)*3.14;
V = A.*L;
all_V = [];
for i = 1:num_gen
    num_branches = 2^i;
    V_tot = V(i)*num_branches*1000;        %total volume of all branches in a generation L
    all_V = [all_V V_tot];                 %add that volume to the vector
end

if num_gen >=2
    for i = length(all_V):-1:1             %total vol for each generation entire network
        all_V = [all_V all_V(i)];
    end
else
    all_V = all_V;
end
network_V = sum(all_V);                    %Whole volume of vasculature in branching network

%------------------------------------- calculate flow in each generation ------------------------------------------------
all_Q = zeros(7,16);
for j = 1:7                                %for each organ in the system
    Q = [74.83758415 138.673326 763.9168803 842.1600 224.48828292 333 1000]/(60*10^6); %L/sec, flow in each organ
    if num_gen>1
        for i = 1:num_gen
            Q(j) = (4*Q(j))/6.3504;         %calculates flow in a branch of new generation
            all_Q_a(j,i) = Q(j);            %appends matrix with new flow
        end
    else
        Q(j) = Q(j);                        %if only one generation Q stays the same
    end
end


if num_gen >=2                              %if more than one generation mirror the matrix to account for all 32 branches
   all_Q_b = fliplr(all_Q_a);
   all_Q = [all_Q_a all_Q_b];
else
    Q = [74.83758415 138.673326 763.9168803 842.1600 224.48828292 333 1000]/(60*10^6);
    for i=1:7
        all_Q = [Q(j)];                     %otherwise Q stays the same
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Physiological Parameter Values from literature   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_NC_b  = (100*10^(-9)); % m, length of np cell receptor bond
      
V_T   = [33.26114851 16.03516137 56.02057122 75.3542703 7.91726 20 1000]/10^6; %L, total tissue volume
         
V_vein  =  466.9/10^6; % L, volume of vein and also artery

Q = [74.83758415 138.673326 763.9168803 284.6716878 224.48828292 333 1000]/(60*10^6); %L/sec, flow in each organ

Q_vein = sum(Q); %L/sec, flow in vein and artery

%% %%%%%%%%%%%%%%%%%%%%%%%%          Rate parameters from literature          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_on  = [120.9 120.9 120.9 120.9 120.9 120.9 120.9]; % Binding rate of NP to endothelial cell surface receptors

K_EC    = [exp(96.96) exp(2.37) exp(10.64) exp(28.46) exp(50.4) exp(35.74) exp(35.74)]; %From Ramakrishnan 2016

K_off = [1.25 51.1 11.36 4.25 2.4 3.38 3.38]; %kon/ln(kec), the rate of NP unbinding from surface receptors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Rate Parameters from sensitivity analysis      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organs = [lung heart kidney liver spleen gut others vein/art]

K_deg = [1/500000 1/500000 1/2000 1/200000 1/20000 1/30000 1/30000 1/6000]; %rate of NP degradation in each organ compartment
         
K_up  = [1 1/10 1/300 1/200 1/10 1/100 1/100]; %rate of NP uptake via receptor mediated endocytosis
    
K_NS  = [1/500 1/500 1/500 1/50 1/70 1/50 1/500]];  % rate of nonspecific uptake via transyctosis

%%%%%%%%%%%%%%%%%%%%%%%% Determining the number of branching elements per organ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%replicate all_v and all_sa
for i = 2:7
    all_V(i,:) = all_V(i-1,:);
    all_SA(i,:) = all_SA(i-1,:);
end 

V_i_BL =[43.65525742 20.04395171  1018.55584 127.7673516 2.226207538 90  700]/10^6;  %L, volume of blood in each organ
phi =([43.65525742 20.04395171  1018.55584 127.7673516 2.226207538 90  700]/10^6)/network_V; %L number of branching elements in the each compartment compartment (element volume/volume of blood in lung)

%recalculate the volume, surface, area, and rate parameters scaled with phi value
for i = 1:length(phi)
    all_V(i,:) = phi(i).*all_V(i,:);
    all_SA(i,:) = phi(i).*all_SA(i,:);

    K_on(i) = phi(i)*K_on(i);
    K_EC(i) = phi(i)*K_EC(i);
    K_off(i) = phi(i)*(K_on(i)./(log(K_EC(i))));
    K_deg(i) = phi(i)*K_deg(i);
    K_up(i) = phi(i)*K_up(i);
    K_NS(i) = phi(i)*K_NS(i);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%       Running ODE simulation       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_step = 10^-9; %time seconds

[t,y] = ode15s(@odefun1,(0:t_step:(10^-3)),y0); %run ode simulation for 1 ms

conc=0.00000000039; %define inital np concentration

%%%%%%%%%%%%%%%%%%%%%%% Plot biodistribuiton of NP in each organ compartment %%%%%%%%%%%%%%%%%%%%%%%%

num_gen = 32;

%Plot NP concentration in artery and vein
figure(1)
subplot(8,1,1)
semilogx((t/3600),y(:,1)/conc, 'o')
hold on
semilogx((t/3600),y(:,2)/conc, 'o')
legend('C_{vein}','C_{artery}')
title('Conc. in Vein & Artery vs time')

%Plot NP concentration in Lung
subplot(8,1,2)
semilogx(t/3600,y(:,451)/conc)
hold on
legend('C_{Lung}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Heart
subplot(8,1,3)
semilogx(t/3600,y(:,452)/conc)
hold on
legend('C_{Heart}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Kidneys
subplot(8,1,4)
semilogx(t/3600,y(:,453)/conc)
hold on
legend('C_{Kidneys}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Liver
subplot(8,1,5)
semilogx(t/3600,y(:,454)/conc)
hold on
legend('C_{Liver}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Spleen
subplot(8,1,6)
semilogx(t/3600,y(:,455)/conc)
hold on
legend('C_{Spleen}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Gut
subplot(8,1,7)
semilogx(t/3600,y(:,456)/conc)
hold on
legend('C_{Gut}^{T}')
title('Normalized Conc of NP')

%Plot NP concentration in Other
subplot(8,1,8)
semilogx(t/3600,y(:,457)/conc)
hold on
legend('C_{Other}^{T}')
title('Normalized Conc of NP')


%%%%%%%%%%%%%%%%%%%% Sum total number of NP degraded in the system to esure mass conservation %%%%%%%%%%%%%%%%%%%%%%%%%


%Calculating total NP degraded in the arteries and veins
M_total_art_ven = 0;
M_total_art_ven_vec = zeros(length(t),1);
for i = 1:length(t)
    M_total_gen_art_ven = 0;
    for j = 1:2
        M_total_art_ven_j = y(i,j)*V_vein*t_step*K_deg(8);
        M_total_gen_art_ven = M_total_gen_art_ven + M_total_art_ven_j;
    end
    M_total_art_ven = M_total_art_ven + M_total_gen_art_ven;
    M_total_art_ven_vec(i) = M_total_art_ven;
end


%Calculating total NP degraded in the blood (vascular branches)
M_total_bl = 0;
M_total_bl_vec = zeros(length(t),1);
for i=1:length(t) %for time
    M_total_gen = 0;
    for j = 3:34 %for branches lung
            M_total_vas_j = y(i,j)*all_V(1,j-2)*t_step*K_deg(1);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 35:66 %for branches heart
            M_total_vas_j = y(i,j)*all_V(2,j-34)*t_step*K_deg(2);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 67:98 %for branches kidneys
            M_total_vas_j = y(i,j)*all_V(3,j-66)*t_step*K_deg(3);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 99:130 %for branches liver
            M_total_vas_j = y(i,j)*all_V(4,j-98)*t_step*K_deg(4);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 131:162 %for branches spleen
            M_total_vas_j = y(i,j)*all_V(5,j-130)*t_step*K_deg(5);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 163:194 %for branches gut 
            M_total_vas_j = y(i,j)*all_V(6,j-162)*t_step*K_deg(6);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    for j = 195:226 %for branches other
            M_total_vas_j = y(i,j)*all_V(7,j-194)*t_step*K_deg(7);
            M_total_gen = M_total_gen + M_total_vas_j;
    end
    M_total_gen = M_total_gen + M_total_vas_j;
    M_total_bl_vec(i) = M_total_gen;
end 

%Calculating total NP degraded in endothelial cells of all branches
M_total_star = 0;
M_total_star_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_gen_EC = 0;
    for j = 227:258
        M_total_j_star = y(i,j)*all_SA(1,j-226)*l_NC_b*t_step*K_deg(1); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 259:290
        M_total_j_star = y(i,j)*all_SA(2,j-258)*l_NC_b*t_step*K_deg(2); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 291:322
        M_total_j_star = y(i,j)*all_SA(3,j-290)*l_NC_b*t_step*K_deg(3); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 323:354
        M_total_j_star = y(i,j)*all_SA(4,j-322)*l_NC_b*t_step*K_deg(4); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 355:386
        M_total_j_star = y(i,j)*all_SA(5,j-354)*l_NC_b*t_step*K_deg(5); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 387:418
        M_total_j_star = y(i,j)*all_SA(6,j-386)*l_NC_b*t_step*K_deg(6); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    for j = 419:450
        M_total_j_star = y(i,j)*all_SA(7,j-418)*l_NC_b*t_step*K_deg(7); 
        M_total_gen_EC= M_total_gen_EC + M_total_j_star;
    end
    M_total_star= M_total_star + M_total_gen_EC;
    M_total_star_vec(i) = M_total_star;
end

%Calculating the total NP degraded in the tissue compartments
M_total_T = 0;
M_total_T_vec = zeros(length(t),1);
for i=1:length(t)
    M_total_gen_T = 0;
    for j= 451
        M_total_j_T = y(i,j)*V_T(1)*t_step*K_deg(1);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 452
        M_total_j_T = y(i,j)*V_T(2)*t_step*K_deg(2);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 453
        M_total_j_T = y(i,j)*V_T(3)*t_step*K_deg(3);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 454
        M_total_j_T = y(i,j)*V_T(4)*t_step*K_deg(4);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 455
        M_total_j_T = y(i,j)*V_T(5)*t_step*K_deg(5);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 456
        M_total_j_T = y(i,j)*V_T(6)*t_step*K_deg(6);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    for j= 457
        M_total_j_T = y(i,j)*V_T(7)*t_step*K_deg(7);
        M_total_gen_T = M_total_gen_T + M_total_j_T;
    end
    M_total_T = M_total_T + M_total_gen_T;      % others
    M_total_T_vec(i) = M_total_T;
end

%%%%%%%%%%%% Calculating total NP in the vasculature of each branch %%%%%%%%%%%%%%%%%%
Vas_total = 0;
for i=3:34
    Vas_total = Vas_total + y(:,i)*all_V(1,i-2);
end
for i=35:66
    Vas_total = Vas_total + y(:,i)*all_V(2,i-34);
end
for i=67:98
    Vas_total = Vas_total + y(:,i)*all_V(3,i-66);
end
for i=99:130
    Vas_total = Vas_total + y(:,i)*all_V(4,i-98);
end
for i=131:162
    Vas_total = Vas_total + y(:,i)*all_V(5,i-130);
end
for i=163:194
    Vas_total = Vas_total + y(:,i)*all_V(6,i-162);
end
for i=195:226
    Vas_total = Vas_total + y(:,i)*all_V(7,i-194);
end


%%%%%%%%%%%% Calculating total NP in the endothelial cells of each branch %%%%%%%%%%%%%%%%%%
EC_total = 0;
for i=227:258
    EC_total = EC_total + y(:,i)*all_SA(1,i-226)*l_NC_b;
end
for i=259:290
    EC_total = EC_total + y(:,i)*all_SA(2,i-258)*l_NC_b;
end
for i=291:322
    EC_total = EC_total + y(:,i)*all_SA(3,i-290)*l_NC_b;
end
for i=323:354
    EC_total = EC_total + y(:,i)*all_SA(4,i-322)*l_NC_b;
end
for i=355:386
    EC_total = EC_total + y(:,i)*all_SA(5,i-354)*l_NC_b;
end
for i=387:418
    EC_total = EC_total + y(:,i)*all_SA(6,i-386)*l_NC_b;
end
for i=419:450
    EC_total = EC_total + y(:,i)*all_SA(7,i-418)*l_NC_b;
end


%%%%%%%%%%%% Calculating total NP in the tissue of each branch %%%%%%%%%%%%%%%%%%
T_total = 0;
for i=451
    T_total = T_total + y(:,i)*V_T(1);
end
for i=452
    T_total = T_total + y(:,i)*V_T(2);
end
for i=453
    T_total = T_total + y(:,i)*V_T(3);
end
for i=454
    T_total = T_total + y(:,i)*V_T(4);
end
for i=455
    T_total = T_total + y(:,i)*V_T(5);
end
for i=456
    T_total = T_total + y(:,i)*V_T(6);
end
for i=457
    T_total = T_total + y(:,i)*V_T(7);
end


%%%%%%%%%%%%%%%%%% This figure shows mass conservation %%%%%%%%%%%%%%%%%
figure(2)
xlabel('Time(s)')
ylabel('Moles NP')
semilogx(t/3600,(((y(:,1)*V_vein)+(y(:,2)*V_vein)+... %conc in art/ven
       (Vas_total)+...                                %conc in vasc
       (EC_total)+...                                 %conc in EC
       (T_total)+...                                  %conc in Tissue  
       (M_total_art_ven_vec)+...                      %degraded in art/vein
       (M_total_T_vec)+...                            %degraded in tissue
       (M_total_star_vec)+...                         %degraded in *
       (M_total_bl_vec))), 'o')                       %degraded in blood
ylim([1.5*10^-13 2.2*10^-13])
%xlim([0 10])
title('Total Mass Conservation')

%%%%%%%%%%%%%%%%% This figure shows branch Diameter vs generation %%%%%%%%%%%%%%%%%
figure(3)
xlabel('Generation')
ylabel('Diameter (m)')
hold on
bar(all_D)

%%%%%%%%%%%%%%%%%%%%%% Define the system of 457 ODEs %%%%%%%%%%%%%%%%%%%%%%

function dydt = odefun1(t,y);
global num_gen num_branches_ls all_D all_SA all_V all_Q V_T V_vein
global K_on K_off K_EC K_NS K_up K_deg K_deg_vein l_NC_b Q Kup_term Kns_term Q_vein V_i_BL

dydt = zeros(457,1);
num_equns = 0;
num_gen = 16*2;
Q_hep = sum(Q(4:6));
V_hep = sum(V_i_BL(4:6));

%vein
dydt(1) = (Q(3)*y(num_gen*3+2)+Q_hep*y(num_gen*4+2)+Q(5)*y(num_gen*5+2)+Q(6)*y(num_gen*6+2)+Q(7)*y(num_gen*7+2) - Q_vein*y(1) - y(1)*V_vein*K_deg(8))/V_vein;
num_equns = num_equns + 1;

%arteries
dydt(2) = ((Q_vein*y(num_gen*2+2)-(Q(3)+Q(4)+Q(5)+Q(6)+Q(7))*y(2)-y(2)*V_vein*K_deg(8)))/V_vein;
num_equns = num_equns + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%   Generating Vasculature ODEs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------LUNGS-----------------------------------------
for j = 1 %for lung
    for i=1:num_gen  
        if i == 1
            dydt(num_equns+1) = ((Q(2)*y(num_gen*2+2)) - ...
                                (all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (2 <= i) && (i<= (num_gen/2))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i-1)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (((num_gen/2)+1) <= i) && (i<= (num_gen-1))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if i == num_gen && num_gen~=1
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(Q(j)*y(num_equns+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
    end
end

%----------------------------------------HEART-----------------------------------------
for j = 2 %for heart
    for i=1:num_gen  
        if i == 1
            dydt(num_equns+1) = ((Q_vein*y(1)) +(Q(1)*y(num_gen+2)) - ...
                                (all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (2 <= i) && (i<= (num_gen/2))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i-1)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (((num_gen/2)+1) <= i) && (i<= (num_gen-1))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if i == num_gen && num_gen~=1
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(Q(j)*y(num_equns+1)+Q_vein*y(2)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
    end
end

%----------------------------------------KIDNEYS-----------------------------------------
for j = 3 %for kidney
    for i=1:num_gen  
        if i == 1
            dydt(num_equns+1) = ((Q_vein*y(1)) - ...
                                (all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (2 <= i) && (i<= (num_gen/2))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i-1)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (((num_gen/2)+1) <= i) && (i<= (num_gen-1))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if i == num_gen && num_gen~=1
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(Q(j)*y(num_equns+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
    end
end
%----------------------------------------LIVER-----------------------------------------
for j = 4 %for liver
    for i=1:num_gen  
        if i == 1
            dydt(num_equns+1) = ((Q(5)*y(num_gen*5+2)+Q(6)*y(num_gen*6+2)+Q(7)*y(num_gen*7+2)) - ...
                                (all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                V_hep;
            num_equns = num_equns + 1;
        end
        if (2 <= i) && (i<= (num_gen/2))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i-1)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                V_hep;
            num_equns = num_equns + 1;
        end
        if (((num_gen/2)+1) <= i) && (i<= (num_gen-1))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                V_hep;
            num_equns = num_equns + 1;
        end
        if i == num_gen && num_gen~=1
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(Q_hep*y(num_equns+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                V_hep;
            num_equns = num_equns + 1;
        end
    end
end
%----------------------------------------SPLEEN GUT OTHER-----------------------------------------
for j = 5:7 %for spleen, gut, other
    for i=1:num_gen  
        if i == 1
            dydt(num_equns+1) = ((Q_vein*y(1)) - ...
                                (all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (2 <= i) && (i<= (num_gen/2))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i-1)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if (((num_gen/2)+1) <= i) && (i<= (num_gen-1))
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(all_Q(j,i)*y(num_equns+1)*num_branches_ls(i+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
        if i == num_gen && num_gen~=1
            dydt(num_equns+1) = ((all_Q(j,i-1)*y(num_equns)*num_branches_ls(i)) ...
                                -(Q(j)*y(num_equns+1)) ...
                                -(K_on(j)*y(num_equns+1)*all_V(i)) ...
                                +(K_off(j)*y((j*num_gen)+2+i)*all_SA(i)*l_NC_b) ...
                                -(K_NS(j)*y(num_equns+1)*all_V(i)) ...
                                -(K_deg(j)*y(num_equns+1)*all_V(i))) / ...
                                all_V(i);
            num_equns = num_equns + 1;
        end
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%  Generating EC ODEs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:7 %for lung heart kidney liver spleen gut other
    for i = 1:num_gen
            dydt(num_equns +1) = ((-K_up(j)*y(num_equns+1)*all_SA(i)*l_NC_b) ...
                                 +(K_on(j)*y(2+i)*all_V(i))  ...
                                 -(K_off(j)*y(num_equns+1)*all_SA(i)*l_NC_b)  ...
                                 -(K_deg(j)*y(num_equns+1)*all_SA(i)*l_NC_b)) / ...
                                 (all_SA(i)*l_NC_b);
            num_equns = num_equns + 1;
    end
end


    %%%%%%%%%%%%%%%%%%%%%   Generating Tissue ODEs  %%%%%%%%%%%%%%%%%%%%%%
for j = 1:7 %for each organ
    Kup_term = [];
    Kns_term = [];

    Kup = [];
    Kns = [];

    for i = 1:num_gen
        Kup_gen = K_up(j) * y((j*num_gen)+2+i) * all_SA(i) * l_NC_b;
        Kup = [Kup Kup_gen];
        Kns_gen = K_NS(j) * y(2+i) * all_V(i);
        Kns = [Kns Kns_gen];
    end
    Kup_term = [Kup_term sum(Kup)];
    Kns_term = [Kns_term sum(Kns)];
    
    dydt(num_equns+1) = ((Kup_term) + (Kns_term) - K_deg(j)*y(num_equns+1)*V_T) / ...
                         (V_T);
                     num_equns = num_equns + 1;
end
end
