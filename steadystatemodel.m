%%%%%%% STEADY STATE MODEL %%%%%%%

% This script is for the steady state model, the model adapted from 
% Ramakrishnan 2016 to include nonspecific uptake via transcytosis. 
% kappa p values were otpimized. 

clc
clear
close all

organs = {'lung';'heart';'kidney';'liver';'spleen'};

%-------------------------------------------------------------------
%Constants throughout organs
lscript_ecb = 10^(-7); %m = distance from an EC or macrophage surface receptor 
d_ec = 5*10^(-6); %m = Diameter of endothelial cell
l_cap = 25.8*10^(-6); %m = size of endothelial cell 
c_out = 10^(-15); %M = injected concentration of NP
lscript_mb = 10^(-7); %m = distance from macrophage surface receptor
d_m = 5*10^(-6); %m = diameter of macrophage cell
l_mb = 10*10^(-6); %m = size of macrophage cell

%-------------------------------------------------------------------
% Kappa p (nonspecific binding rate of NP) 
% values taken from the following references and used in fig
references = {'Murciano','Papede','Muro/Hsu','Kohlar', 'Carme'};
kappa_p_lu = [0.79 1.9 22.35 5.0 10.24];
kappa_p_h = [0.913 0.4 3.01 1.16 1.07];
kappa_p_k = [2.10 0.5 4.5 1.11 1.95];
kappa_p_li = [1.48 11 50.95 15.78 56.83];
kappa_p_s = [1.22 15.0 67.96 22.22 63.73];


kappa_p = [nthroot((prod(kappa_p_lu)),5) nthroot((prod(kappa_p_h)), 5)/1000000 nthroot((prod(kappa_p_k)), 5)/1000000 ...
          nthroot((prod(kappa_p_li)), 5)*1000000000000 nthroot((prod(kappa_p_s)), 5)*10000000000000];

%-------------------------------------------------------------------
%Phi EC values
phi_ec = [.3 .3 .3 .3 .3]; % Concentration of EC in target tissue
phi_m = [.03 .03 .03 .03 .03]; % Concentration of macrophage cells in target tissue
phi_ma = [.03 .03 .03 .03 .03]; % 

%-------------------------------------------------------------------
%association constants for different antibody levels in different organs
k_ec_flat_low = [4.65*10^18 3.285060633 7100.073303 22595869600 1.24*10^13];
k_ec_flat_med = [4.39*10^25 2.832680049 3083.202198 2423458474 3.67*10^12];
k_ec_flat_high = [2.89*10^23 3.526502771 10680.54584 69467039755 2.72*10^13];


k_ec_low = [2.42*10^16 1.88457365 646.227279 45277209.7 1.8877*10^11]; % The association constant for binding to endothelial cells
k_ec_med = [1.4*10^24 4.13474953 7112.57997 2.2692*10^10 2.42*10^15];
k_ec_high = [1.23*10^42 10.6565662 41793.3924 2.28*10^12 7.74*10^21];


k_m_low = [ 47718.5533 47718.5533 47718.5533 47718.5533 47718.5533]; % The association constant for unactivated macrophages
k_m_med = [2460953.89 2460953.89 2460953.89 2460953.89 2460953.89];
k_m_high = [2.55*10^13 2.55*10^13 2.55*10^13 2.55*10^13 2.55*10^13];


k_ma_low = [1216469103 1216469103 1216469103 1216469103 1216469103]; % The association constant for activated macrophages
k_ma_med = [6.84*10^15 6.84*10^15 6.84*10^15 6.84*10^15 6.84*10^15];
k_ma_high = [3.66*10^22 3.66*10^22 3.66*10^22 3.66*10^22 3.66*10^22];

%-------------------------------------------------------------------
%Idg calculations

%flat
for i = 1:length(kappa_p)
    k_ec_low_fl = k_ec_flat_low(i);
    kappa_p_fl = kappa_p(i);
    phi_ec_fl = phi_ec(i);
    percent_idg_fl_low(i) = (kappa_p_fl + (((phi_ec_fl * k_ec_low_fl * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb); 
end
for j = 1:length(kappa_p)
    k_ec_med_fl = k_ec_flat_med(j);
    kappa_p_fl = kappa_p(j);
    phi_ec_fl = phi_ec(j);
    percent_idg_fl_med(j) = (kappa_p_fl + (((phi_ec_fl * k_ec_med_fl * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb);    
end
for k = 1:length(kappa_p)
    k_ec_high_fl = k_ec_flat_high(k);
    kappa_p_fl = kappa_p(k);
    phi_ec_fl = phi_ec(k);
    percent_idg_fl_high(k) = (kappa_p_fl + (((phi_ec_fl * k_ec_high_fl * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb);  
end

%membrane only
for i = 1:length(kappa_p)
    k_ec_low_m = k_ec_low(i);
    kappa_p_m = kappa_p(i);
    phi_ec_m = phi_ec(i);
    percent_idg_m_low(i) = (kappa_p_m + (((phi_ec_m * k_ec_low_m * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb); 
end
for j = 1:length(kappa_p)
    k_ec_med_m = k_ec_med(j);
    kappa_p_m = kappa_p(j);
    phi_ec_m = phi_ec(j);
    percent_idg_m_med(j) = (kappa_p_m + (((phi_ec_m * k_ec_med_m * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb);    
end
for k = 1:length(kappa_p)
    k_ec_high_m = k_ec_high(k);
    kappa_p_m = kappa_p(k);
    phi_ec_m = phi_ec(k);
    percent_idg_m_high(k) = (kappa_p_m + (((phi_ec_m * k_ec_high_m * lscript_ecb)/(d_ec)) * c_out)) * (l_cap/lscript_ecb);  
end
    
%membrane and other cells
for i = 1:length(kappa_p)
    k_ec_low_mc = k_ec_low(i);
    k_m_low_mc = k_m_low(i);
    kappa_p_mc = kappa_p(i);
    phi_ec_mc = phi_ec(i);
    phi_m_mc = phi_m(i);
   percent_idg_mc_low(i) = (kappa_p_mc + (((phi_ec_mc * k_ec_low_mc * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_m_mc * k_m_low_mc * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));
end
for j = 1:length(kappa_p)
    k_ec_med_mc = k_ec_med(j);
    k_m_med_mc = k_m_med(j);
    kappa_p_mc = kappa_p(j);
    phi_ec_mc = phi_ec(j);
    phi_m_mc = phi_m(j);
   percent_idg_mc_med(j) = (kappa_p_mc + (((phi_ec_mc * k_ec_med_mc * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_m_mc * k_m_med_mc * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));
end
for k = 1:length(kappa_p)
    k_ec_high_mc = k_ec_high(k);
    k_m_high_mc = k_m_high(k);
    kappa_p_mc = kappa_p(k);
    phi_ec_mc = phi_ec(k);
    phi_m_mc = phi_m(k);
   percent_idg_mc_high(k) = (kappa_p_mc + (((phi_ec_mc * k_ec_high_mc * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_m_mc * k_m_high_mc * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));
end

%membrane and other cells*
for i = 1:length(kappa_p)
    k_ec_low_mca = k_ec_low(i);
    k_ma_low_mca = k_ma_low(i);
    kappa_p_mca = kappa_p(i);
    phi_ec_mca = phi_ec(i);
    phi_ma_mca = phi_ma(i);
   percent_idg_mca_low(i) = (kappa_p_mca+ (((phi_ec_mca * k_ec_low_mca * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_ma_mca * k_ma_low_mca * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));     
end
for j = 1:length(kappa_p)
    k_ec_med_mca = k_ec_med(j);
    k_ma_med_mca = k_ma_med(j);
    kappa_p_mca = kappa_p(j);
    phi_ec_mca = phi_ec(j);
    phi_ma_mca = phi_ma(j);
   percent_idg_mca_med(j) = (kappa_p_mca + (((phi_ec_mca * k_ec_med_mca * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_ma_mca * k_ma_med_mca * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));      
end
for k = 1:length(kappa_p)
    k_ec_high_mca = k_ec_high(k);
    k_ma_high_mca = k_ma_high(k);
    kappa_p_mca = kappa_p(k);
    phi_ec_mca = phi_ec(k);
    phi_ma_mca = phi_ma(k);
   percent_idg_mca_high(k) = (kappa_p_mca + (((phi_ec_mca * k_ec_high_mca * lscript_ecb)/(d_ec)) * c_out))...
          * (l_cap/lscript_ecb) + ((((phi_ma_mca * k_ma_high_mca * lscript_mb)/(d_m))...
          * c_out) * (l_cap / lscript_mb));     
end

%-------------------------------------------------------------------
% nexp for experimental values taken from fig-4-5-6.py:
nexp_high = [125 2 NaN 72 75]/(125);
nexp_med = [91 2 NaN 60 60]/(125);
nexp_low = [23 1 NaN 85 115]/(125);

%nsim for flat surface
nsim_fl_low = log(percent_idg_fl_low)/log(percent_idg_fl_high(1));
nsim_fl_med = log(percent_idg_fl_med)/log(percent_idg_fl_high(1));
nsim_fl_high = log(percent_idg_fl_high)/log(percent_idg_fl_high(1));

% nsim for membrane
nsim_m_low = log(percent_idg_m_low)/log(percent_idg_m_high(1));
nsim_m_med = log(percent_idg_m_med)/log(percent_idg_m_high(1));
nsim_m_high = log(percent_idg_m_high)/log(percent_idg_m_high(1));

%nsim for membrane and other cells
nsim_mc_low = log(percent_idg_mc_low)/log(percent_idg_mc_high(1));
nsim_mc_med = log(percent_idg_mc_med)/log(percent_idg_mc_high(1));
nsim_mc_high = log(percent_idg_mc_high)/log(percent_idg_mc_high(1));

%nsim for membrane and other cells*
nsim_mca_low = log(percent_idg_mca_low)/log(percent_idg_mca_high(1));
nsim_mca_med = log(percent_idg_mca_med)/log(percent_idg_mca_high(1));
nsim_mca_high = log(percent_idg_mca_high)/log(percent_idg_mca_high(1));

%-------------------------------------------------------------------
%correlation coefficent calculations

nexp_lung = [nexp_low(1) nexp_med(1) nexp_high(1)];
nexp_heart =[nexp_low(2) nexp_med(2) nexp_high(2)];
nexp_liver = [nexp_low(4) nexp_med(4) nexp_high(4)];
nexp_spleen = [nexp_low(5) nexp_med(5) nexp_high(5)];
nexp_total = [nexp_lung;...
              nexp_heart;...
              nexp_liver;...
              nexp_spleen];
          
nsim_lung_fl = [nsim_fl_low(1) nsim_fl_med(1) nsim_fl_high(1)];
nsim_heart_fl = [nsim_fl_low(2) nsim_fl_med(2) nsim_fl_high(2)];
nsim_liver_fl = [nsim_fl_low(4) nsim_fl_med(4) nsim_fl_high(4)];
nsim_spleen_fl = [nsim_fl_low(5) nsim_fl_med(5) nsim_fl_high(5)];
nsim_total_fl = [nsim_lung_fl;...
                 nsim_heart_fl;...
                 nsim_liver_fl;...
                 nsim_spleen_fl];

nsim_lung_m = [nsim_m_low(1) nsim_m_med(1) nsim_m_high(1)];
nsim_heart_m = [nsim_m_low(2) nsim_m_med(2) nsim_m_high(2)];
nsim_liver_m = [nsim_m_low(4) nsim_m_med(4) nsim_m_high(4)];
nsim_spleen_m = [nsim_m_low(5) nsim_m_med(5) nsim_m_high(5)];
nsim_total_m = [nsim_lung_m;...
                nsim_heart_m;...
                nsim_liver_m;...
                nsim_spleen_m];

nsim_lung_mc = [nsim_mc_low(1) nsim_mc_med(1) nsim_mc_high(1)];
nsim_heart_mc = [nsim_mc_low(2) nsim_mc_med(2) nsim_mc_high(2)];
nsim_liver_mc = [nsim_mc_low(4) nsim_mc_med(4) nsim_mc_high(4)];
nsim_spleen_mc = [nsim_mc_low(5) nsim_mc_med(5) nsim_mc_high(5)];
nsim_total_mc = [nsim_lung_mc;...
                 nsim_heart_mc;...
                 nsim_liver_mc;...
                 nsim_spleen_mc];

nsim_lung_mca = [nsim_mca_low(1) nsim_mca_med(1) nsim_mca_high(1)];
nsim_heart_mca = [nsim_mca_low(2) nsim_mca_med(2) nsim_mca_high(2)];
nsim_liver_mca = [nsim_mca_low(4) nsim_mca_med(4) nsim_mca_high(4)];
nsim_spleen_mca = [nsim_mca_low(5) nsim_mca_med(5) nsim_mca_high(5)];
nsim_total_mca = [nsim_lung_mca;...
                  nsim_heart_mca;...
                  nsim_liver_mca;...
                  nsim_spleen_mca];
              
R_lung_fl = corrcoef(nexp_lung, nsim_lung_fl);             
R_lung_m = corrcoef(nexp_lung, nsim_lung_m);
R_lung_mc = corrcoef(nexp_lung, nsim_lung_mc);
R_lung_mca = corrcoef(nexp_lung, nsim_lung_mca);

R_fl = corrcoef(nexp_total, nsim_total_fl);
R_m = corrcoef(nexp_total, nsim_total_m);
R_mc = corrcoef(nexp_total, nsim_total_mc);
R_mca = corrcoef(nexp_total, nsim_total_mca);



%-------------------------------------------------------------------
%Plots created

figure
subplot(5,1,1)
               exp = [nexp_low(1) nexp_med(1) nexp_high(1);...
                      nexp_low(2) nexp_med(2) nexp_high(2);...
                      nexp_low(3) nexp_med(3) nexp_high(3);...
                      nexp_low(4) nexp_med(4) nexp_high(4);...
                      nexp_low(5) nexp_med(5) nexp_high(5)];
bar(exp)
title('\it in vivo')
ylabel('\eta_{exp}')
set(gca,'xtick',[1:5],'xticklabel',organs)
axis([.5 5.5 0 1.1])
legend({'41','100', '162'},'Location','best')

subplot(5,1,2)
flat = [nsim_fl_low(1) nsim_fl_med(1) nsim_fl_high(1);...
                      nsim_fl_low(2) nsim_fl_med(2) nsim_fl_high(2);...
                      nsim_fl_low(3) nsim_fl_med(3) nsim_fl_high(3);...
                      nsim_fl_low(4) nsim_fl_med(4) nsim_fl_high(4);...
                      nsim_fl_low(5) nsim_fl_med(5) nsim_fl_high(5)];
bar(flat)
title('Flat')
ylabel('\eta_{sim}')
set(gca,'xtick',[1:5],'xticklabel',organs)
axis([.5 5.5 0 2.2])
dim = [.2 .5 .3 .3];
R_lung_fl_txt=['  R Lung = ',num2str((R_lung_fl(1,2)))];
R_fl_txt = ['  R Total = ',num2str((R_fl(1,2)))];
R_lung_fl_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_lung_fl_txt); 
R_fl_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_fl_txt); 
set(R_lung_fl_T, 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
set(R_fl_T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(5,1,3)
membrane = [nsim_m_low(1) nsim_m_med(1) nsim_m_high(1);...
            nsim_m_low(2) nsim_m_med(2) nsim_m_high(2);...
            nsim_m_low(3) nsim_m_med(3) nsim_m_high(3);...
            nsim_m_low(4) nsim_m_med(4) nsim_m_high(4);...
            nsim_m_low(5) nsim_m_med(5) nsim_m_high(5)];
bar(membrane)
title('Membrane')
ylabel('\eta_{sim}')
set(gca,'xtick',[1:5],'xticklabel',organs)
axis([.5 5.5 0 1.2])
R_lung_m_txt=['  R Lung = ',num2str((R_lung_m(1,2)))];
R_m_txt = ['  R Total = ',num2str((R_m(1,2)))];
R_lung_m_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_lung_m_txt); 
R_m_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_m_txt); 
set(R_lung_m_T, 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
set(R_m_T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(5,1,4)
membrane_cells = [nsim_mc_low(1) nsim_mc_med(1) nsim_mc_high(1);...
                  nsim_mc_low(2) nsim_mc_med(2) nsim_mc_high(2);...
                  nsim_mc_low(3) nsim_mc_med(3) nsim_mc_high(3);...
                  nsim_mc_low(4) nsim_mc_med(4) nsim_mc_high(4);...
                  nsim_mc_low(5) nsim_mc_med(5) nsim_mc_high(5)];
bar(membrane_cells)
title('Membrane + Other Cells')
ylabel('\eta_{sim}')
set(gca,'xtick',[1:5],'xticklabel',organs)
axis([.5 5.5 0 1.7])
R_lung_mc_txt=['  R Lung = ',num2str((R_lung_mc(1,2)))];
R_mc_txt = ['  R Total = ',num2str((R_mc(1,2)))];
R_lung_mc_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_lung_mc_txt); 
R_mc_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_mc_txt); 
set(R_lung_mc_T, 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
set(R_mc_T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');


subplot(5,1,5)
membrane_cells_act = [nsim_mca_low(1) nsim_mca_med(1) nsim_mca_high(1);...
                      nsim_mca_low(2) nsim_mca_med(2) nsim_mca_high(2);...
                      nsim_mca_low(3) nsim_mca_med(3) nsim_mca_high(3);...
                      nsim_mca_low(4) nsim_mca_med(4) nsim_mca_high(4);...
                      nsim_mca_low(5) nsim_mca_med(5) nsim_mca_high(5)];
bar(membrane_cells_act)
title('Membrane + Other Cells*')
ylabel('\eta_{sim}')
xlabel('organ')
set(gca,'xtick',[1:5],'xticklabel',organs)
axis([.5 5.5 0 1.7])
dim = [.2 .5 .3 .3];
R_lung_mca_txt=['  R Lung = ',num2str((R_lung_mca(1,2)))];
R_mca_txt = ['  R Total = ',num2str((R_mca(1,2)))];
R_lung_mca_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_lung_mca_txt); 
R_mca_T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), R_mca_txt); 
set(R_lung_mca_T, 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'left');
set(R_mca_T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

