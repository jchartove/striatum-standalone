function [T,soma_V,soma_somaGolombK_a,soma_somaGolombK_b,soma_somaGolombKdr_n,soma_somaGolombNa_h,dend_V,dend_dendGolombK_a,dend_dendGolombK_b,dend_dendGolombKdr_n,dend_dendGolombNa_h,D1_V,D1_naCurrentMSN_m,D1_naCurrentMSN_h,D1_kCurrentMSN_m,D1_mCurrentMSN_m,D2_V,D2_naCurrentMSN_m,D2_naCurrentMSN_h,D2_kCurrentMSN_m,D2_mCurrentMSN_m,soma_soma_somaSomaiSYN_s,D1_soma_somaMSNiSYN_s,D1_D1_gabaRecInputMSN_s,D2_soma_somaMSNiSYN_s,D2_D2_gabaRecInputMSN_s,soma_somaGolombK_gd2,soma_somaLeak_gl2,dend_dendGolombK_gd_dend,dend_dendGolombKdr_gkdr_dend,dend_dendGolombNa_gna_dend,dend_dendInput_tonic2,dend_dendLeak_gl_dend,dend_dendiMultiPoissonExp_Ge,dend_dendiMultiPoissonExp_Gi,D1_naCurrentMSN_V_IC,D1_kCurrentMSN_V_IC,D1_mCurrentMSN_Qs,D1_mCurrentMSN_V_IC,D1_injectedCurrentD1_freq,D1_AMPAMSN_t1,D1_AMPAMSN_AMPA_onset2,D1_AMPAMSN_psp,D1_AMPAMSN_psp2,D1_AMPAMSN_halfpop,D1_AMPAMSN_psptime,D1_AMPAMSN_psptime2,D1_AMPAMSN_cellmask1,D1_AMPAMSN_cellmask2,D2_naCurrentMSN_V_IC,D2_kCurrentMSN_V_IC,D2_mCurrentMSN_Qs,D2_mCurrentMSN_V_IC,D2_injectedCurrentD2_freq,D2_AMPAMSN_t1,D2_AMPAMSN_AMPA_onset2,D2_AMPAMSN_psp,D2_AMPAMSN_psp2,D2_AMPAMSN_halfpop,D2_AMPAMSN_psptime,D2_AMPAMSN_psptime2,D2_AMPAMSN_cellmask1,D2_AMPAMSN_cellmask2,soma_soma_somaSomaiSYN_gsyn2,soma_soma_somaSomaiSYN_mask,dend_dend_dendDendiGAP_g_GAP2,dend_dend_dendDendiGAP_mask,D1_soma_somaMSNiSYN_indegree,D1_soma_somaMSNiSYN_mask,D1_D1_gabaRecInputMSN_netcon,D2_soma_somaMSNiSYN_indegree,D2_soma_somaMSNiSYN_mask,D2_D2_gabaRecInputMSN_netcon]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);


% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
% seed the random number generator
rng_wrapper(p.random_seed);
soma_somaGolombK_gd2 =  p.soma_somaGolombK_gd + p.soma_somaGolombK_gd.*p.soma_somaGolombK_gd_het.*2.*(rand(1,p.soma_Npop) - 0.5);
soma_somaLeak_gl2 =  p.soma_somaLeak_gl + p.soma_somaLeak_gl.*p.soma_somaLeak_gl_het.*2.*(rand(1,p.soma_Npop) - 0.5);
dend_dendGolombK_gd_dend =  p.dend_dendGolombK_gd/10 + p.dend_dendGolombK_gd.*p.dend_dendGolombK_gd_het.*2.*(0.1.*rand(1,p.dend_Npop) - 0.05);
dend_dendGolombKdr_gkdr_dend =  p.dend_dendGolombKdr_gkdr/10;
dend_dendGolombNa_gna_dend =  p.dend_dendGolombNa_gna/10;
dend_dendInput_tonic2 =  (p.dend_dendInput_tonic + 7.*p.dend_dendInput_DA).*(1 + p.dend_dendInput_tonic_het.*2.*(rand(1,p.dend_Npop) - 0.5));
dend_dendLeak_gl_dend =  p.dend_dendLeak_gl/10 + p.dend_dendLeak_gl.*p.dend_dendLeak_gl_het.*2.*(0.1.*rand(1,p.dend_Npop) - 0.05);
dend_dendiMultiPoissonExp_Ge =  corrPoisson(p.dend_Npop, p.dend_dendiMultiPoissonExp_N_einputs, p.dend_dendiMultiPoissonExp_rate, p.dend_dendiMultiPoissonExp_tau_i, p.dend_dendiMultiPoissonExp_tau_1, 2, .5, p.dend_dendiMultiPoissonExp_Tfinal, dt, p.dend_dendiMultiPoissonExp_fraction_shared);
dend_dendiMultiPoissonExp_Gi =  corrPoisson(p.dend_Npop, p.dend_dendiMultiPoissonExp_N_iinputs, p.dend_dendiMultiPoissonExp_rate, p.dend_dendiMultiPoissonExp_tau_i, p.dend_dendiMultiPoissonExp_tau_1, 5, .5, p.dend_dendiMultiPoissonExp_Tfinal, dt, p.dend_dendiMultiPoissonExp_fraction_shared);
D1_naCurrentMSN_V_IC = -6.300000000000000e+01;
D1_kCurrentMSN_V_IC = -6.300000000000000e+01;
D1_mCurrentMSN_Qs =  p.D1_mCurrentMSN_Q10^(.1*(37-23));
D1_mCurrentMSN_V_IC = -6.300000000000000e+01;
D1_injectedCurrentD1_freq =  1/(2*pi);
D1_AMPAMSN_t1 =  0:p.D1_AMPAMSN_Tfinal;
D1_AMPAMSN_AMPA_onset2 =  p.D1_AMPAMSN_AMPA_onset + p.D1_AMPAMSN_onset2_delay;
D1_AMPAMSN_psp =  p.D1_AMPAMSN_tau_i*(exp(-max(D1_AMPAMSN_t1 - p.D1_AMPAMSN_tau_1,0)/p.D1_AMPAMSN_tau_d) - exp(-max(D1_AMPAMSN_t1 - p.D1_AMPAMSN_tau_1,0)/p.D1_AMPAMSN_tau_r))/(p.D1_AMPAMSN_tau_d - p.D1_AMPAMSN_tau_r);
D1_AMPAMSN_psp2 =  D1_AMPAMSN_psp(D1_AMPAMSN_psp > eps);
D1_AMPAMSN_halfpop =  p.D1_Npop/2;
D1_AMPAMSN_psptime =  horzcat(zeros(1,p.D1_AMPAMSN_AMPA_onset),D1_AMPAMSN_psp2,zeros(1,p.D1_AMPAMSN_Tfinal-(p.D1_AMPAMSN_AMPA_onset+length(D1_AMPAMSN_psp2))+1));
D1_AMPAMSN_psptime2 =  horzcat(zeros(1,D1_AMPAMSN_AMPA_onset2),D1_AMPAMSN_psp2,zeros(1,p.D1_AMPAMSN_Tfinal-(D1_AMPAMSN_AMPA_onset2+length(D1_AMPAMSN_psp2))+1));
D1_AMPAMSN_cellmask1 =  horzcat(zeros(1,D1_AMPAMSN_halfpop),ones(1,D1_AMPAMSN_halfpop));
D1_AMPAMSN_cellmask2 =  horzcat(zeros(1,D1_AMPAMSN_halfpop/2),ones(1,D1_AMPAMSN_halfpop),zeros(1,D1_AMPAMSN_halfpop/2));
D2_naCurrentMSN_V_IC = -6.300000000000000e+01;
D2_kCurrentMSN_V_IC = -6.300000000000000e+01;
D2_mCurrentMSN_Qs =  p.D2_mCurrentMSN_Q10^(.1*(37-23));
D2_mCurrentMSN_V_IC = -6.300000000000000e+01;
D2_injectedCurrentD2_freq =  1/(2*pi);
D2_AMPAMSN_t1 =  0:p.D2_AMPAMSN_Tfinal;
D2_AMPAMSN_AMPA_onset2 =  p.D2_AMPAMSN_AMPA_onset + p.D2_AMPAMSN_onset2_delay;
D2_AMPAMSN_psp =  p.D2_AMPAMSN_tau_i*(exp(-max(D2_AMPAMSN_t1 - p.D2_AMPAMSN_tau_1,0)/p.D2_AMPAMSN_tau_d) - exp(-max(D2_AMPAMSN_t1 - p.D2_AMPAMSN_tau_1,0)/p.D2_AMPAMSN_tau_r))/(p.D2_AMPAMSN_tau_d - p.D2_AMPAMSN_tau_r);
D2_AMPAMSN_psp2 =  D2_AMPAMSN_psp(D2_AMPAMSN_psp > eps);
D2_AMPAMSN_halfpop =  p.D2_Npop/2;
D2_AMPAMSN_psptime =  horzcat(zeros(1,p.D2_AMPAMSN_AMPA_onset),D2_AMPAMSN_psp2,zeros(1,p.D2_AMPAMSN_Tfinal-(p.D2_AMPAMSN_AMPA_onset+length(D2_AMPAMSN_psp2))+1));
D2_AMPAMSN_psptime2 =  horzcat(zeros(1,D2_AMPAMSN_AMPA_onset2),D2_AMPAMSN_psp2,zeros(1,p.D2_AMPAMSN_Tfinal-(D2_AMPAMSN_AMPA_onset2+length(D2_AMPAMSN_psp2))+1));
D2_AMPAMSN_cellmask1 =  horzcat(zeros(1,D2_AMPAMSN_halfpop),ones(1,D2_AMPAMSN_halfpop));
D2_AMPAMSN_cellmask2 =  horzcat(zeros(1,D2_AMPAMSN_halfpop/2),ones(1,D2_AMPAMSN_halfpop),zeros(1,D2_AMPAMSN_halfpop/2));
soma_soma_somaSomaiSYN_gsyn2 =  p.soma_soma_somaSomaiSYN_gsyn - p.soma_soma_somaSomaiSYN_DA.*(p.soma_soma_somaSomaiSYN_gsyn - 0.005);
soma_soma_somaSomaiSYN_mask =  genmask(p.soma_Npop,p.soma_Npop,p.soma_soma_somaSomaiSYN_i_con,soma_soma_somaSomaiSYN_gsyn2,1,1,0);
dend_dend_dendDendiGAP_g_GAP2 =  p.dend_dend_dendDendiGAP_g_GAP + p.dend_dend_dendDendiGAP_DA.*(0.3-p.dend_dend_dendDendiGAP_g_GAP);
dend_dend_dendDendiGAP_mask =  genmask(p.dend_Npop,p.dend_Npop,p.dend_dend_dendDendiGAP_gcon,dend_dend_dendDendiGAP_g_GAP2,0,1,0);
D1_soma_somaMSNiSYN_indegree =  p.D1_soma_somaMSNiSYN_i_con*(p.soma_Npop-p.D1_soma_somaMSNiSYN_ko);
D1_soma_somaMSNiSYN_mask =  genmask(p.soma_Npop,p.D1_Npop,p.D1_soma_somaMSNiSYN_i_con,p.D1_soma_somaMSNiSYN_m_gsyn,1,0,p.D1_soma_somaMSNiSYN_ko);
D1_D1_gabaRecInputMSN_netcon =  ones(p.D1_Npop)-eye(p.D1_Npop);
D2_soma_somaMSNiSYN_indegree =  p.D2_soma_somaMSNiSYN_i_con*(p.soma_Npop-p.D2_soma_somaMSNiSYN_ko);
D2_soma_somaMSNiSYN_mask =  genmask(p.soma_Npop,p.D2_Npop,p.D2_soma_somaMSNiSYN_i_con,p.D2_soma_somaMSNiSYN_m_gsyn,1,0,p.D2_soma_somaMSNiSYN_ko);
D2_D2_gabaRecInputMSN_netcon =  ones(p.D2_Npop)-eye(p.D2_Npop);


% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng_wrapper(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
soma_V_last = -90 * ones(1,p.soma_Npop);
soma_V = zeros(nsamp,p.soma_Npop);
soma_V(1,:) = soma_V_last;
soma_somaGolombK_a_last =  0.15 + 0.60*rand(1,p.soma_Npop);
soma_somaGolombK_a = zeros(nsamp,p.soma_Npop);
soma_somaGolombK_a(1,:) = soma_somaGolombK_a_last;
soma_somaGolombK_b_last =  0.15 + 0.45*rand(1,p.soma_Npop);
soma_somaGolombK_b = zeros(nsamp,p.soma_Npop);
soma_somaGolombK_b(1,:) = soma_somaGolombK_b_last;
soma_somaGolombKdr_n_last =  0.05 + 0.45*rand(1,p.soma_Npop);
soma_somaGolombKdr_n = zeros(nsamp,p.soma_Npop);
soma_somaGolombKdr_n(1,:) = soma_somaGolombKdr_n_last;
soma_somaGolombNa_h_last =  0.05 + 0.85*rand(1,p.soma_Npop);
soma_somaGolombNa_h = zeros(nsamp,p.soma_Npop);
soma_somaGolombNa_h(1,:) = soma_somaGolombNa_h_last;
dend_V_last = -90 * ones(1,p.dend_Npop);
dend_V = zeros(nsamp,p.dend_Npop);
dend_V(1,:) = dend_V_last;
dend_dendGolombK_a_last =  0.15 + 0.60*rand(1,p.dend_Npop);
dend_dendGolombK_a = zeros(nsamp,p.dend_Npop);
dend_dendGolombK_a(1,:) = dend_dendGolombK_a_last;
dend_dendGolombK_b_last =  0.15 + 0.45*rand(1,p.dend_Npop);
dend_dendGolombK_b = zeros(nsamp,p.dend_Npop);
dend_dendGolombK_b(1,:) = dend_dendGolombK_b_last;
dend_dendGolombKdr_n_last =  0.05 + 0.45*rand(1,p.dend_Npop);
dend_dendGolombKdr_n = zeros(nsamp,p.dend_Npop);
dend_dendGolombKdr_n(1,:) = dend_dendGolombKdr_n_last;
dend_dendGolombNa_h_last =  0.05 + 0.85*rand(1,p.dend_Npop);
dend_dendGolombNa_h = zeros(nsamp,p.dend_Npop);
dend_dendGolombNa_h(1,:) = dend_dendGolombNa_h_last;
D1_V_last = -63 * ones(1,p.D1_Npop);
D1_V = zeros(nsamp,p.D1_Npop);
D1_V(1,:) = D1_V_last;
D1_naCurrentMSN_m_last =  0.32*(D1_naCurrentMSN_V_IC+54)/(1-exp(-(D1_naCurrentMSN_V_IC+54)/4))/(0.32*(D1_naCurrentMSN_V_IC+54)/(1-exp(-(D1_naCurrentMSN_V_IC+54)/4))+0.28*(D1_naCurrentMSN_V_IC+27)/(exp((D1_naCurrentMSN_V_IC+27)/5)-1))*randn(1,p.D1_Npop);
D1_naCurrentMSN_m = zeros(nsamp,p.D1_Npop);
D1_naCurrentMSN_m(1,:) = D1_naCurrentMSN_m_last;
D1_naCurrentMSN_h_last =  0.128*exp(-(D1_naCurrentMSN_V_IC+50)/18)/(0.128*exp(-(D1_naCurrentMSN_V_IC+50)/18)+4/(1+exp(-(D1_naCurrentMSN_V_IC+27)/5)))*randn(1,p.D1_Npop);
D1_naCurrentMSN_h = zeros(nsamp,p.D1_Npop);
D1_naCurrentMSN_h(1,:) = D1_naCurrentMSN_h_last;
D1_kCurrentMSN_m_last =  0.032*(D1_kCurrentMSN_V_IC+52)/(1-exp(-(D1_kCurrentMSN_V_IC+52)/5))/(0.032*(D1_kCurrentMSN_V_IC+52)/(1-exp(-(D1_kCurrentMSN_V_IC+52)/5))+0.5*exp(-(D1_kCurrentMSN_V_IC+57)/40))*randn(1,p.D1_Npop);
D1_kCurrentMSN_m = zeros(nsamp,p.D1_Npop);
D1_kCurrentMSN_m(1,:) = D1_kCurrentMSN_m_last;
D1_mCurrentMSN_m_last =  p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp(-(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9))/(p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp(-(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9))-p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp((D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9)))*randn(1,p.D1_Npop);
D1_mCurrentMSN_m = zeros(nsamp,p.D1_Npop);
D1_mCurrentMSN_m(1,:) = D1_mCurrentMSN_m_last;
D2_V_last = -63 * ones(1,p.D2_Npop);
D2_V = zeros(nsamp,p.D2_Npop);
D2_V(1,:) = D2_V_last;
D2_naCurrentMSN_m_last =  0.32*(D2_naCurrentMSN_V_IC+54)/(1-exp(-(D2_naCurrentMSN_V_IC+54)/4))/(0.32*(D2_naCurrentMSN_V_IC+54)/(1-exp(-(D2_naCurrentMSN_V_IC+54)/4))+0.28*(D2_naCurrentMSN_V_IC+27)/(exp((D2_naCurrentMSN_V_IC+27)/5)-1))*randn(1,p.D2_Npop);
D2_naCurrentMSN_m = zeros(nsamp,p.D2_Npop);
D2_naCurrentMSN_m(1,:) = D2_naCurrentMSN_m_last;
D2_naCurrentMSN_h_last =  0.128*exp(-(D2_naCurrentMSN_V_IC+50)/18)/(0.128*exp(-(D2_naCurrentMSN_V_IC+50)/18)+4/(1+exp(-(D2_naCurrentMSN_V_IC+27)/5)))*randn(1,p.D2_Npop);
D2_naCurrentMSN_h = zeros(nsamp,p.D2_Npop);
D2_naCurrentMSN_h(1,:) = D2_naCurrentMSN_h_last;
D2_kCurrentMSN_m_last =  0.032*(D2_kCurrentMSN_V_IC+52)/(1-exp(-(D2_kCurrentMSN_V_IC+52)/5))/(0.032*(D2_kCurrentMSN_V_IC+52)/(1-exp(-(D2_kCurrentMSN_V_IC+52)/5))+0.5*exp(-(D2_kCurrentMSN_V_IC+57)/40))*randn(1,p.D2_Npop);
D2_kCurrentMSN_m = zeros(nsamp,p.D2_Npop);
D2_kCurrentMSN_m(1,:) = D2_kCurrentMSN_m_last;
D2_mCurrentMSN_m_last =  p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp(-(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9))/(p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp(-(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9))-p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp((D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9)))*randn(1,p.D2_Npop);
D2_mCurrentMSN_m = zeros(nsamp,p.D2_Npop);
D2_mCurrentMSN_m(1,:) = D2_mCurrentMSN_m_last;
soma_soma_somaSomaiSYN_s_last =  p.soma_soma_somaSomaiSYN_IC+p.soma_soma_somaSomaiSYN_IC_noise.*rand(1,p.soma_Npop);
soma_soma_somaSomaiSYN_s = zeros(nsamp,p.soma_Npop);
soma_soma_somaSomaiSYN_s(1,:) = soma_soma_somaSomaiSYN_s_last;
D1_soma_somaMSNiSYN_s_last =  p.D1_soma_somaMSNiSYN_IC+p.D1_soma_somaMSNiSYN_IC_noise.*rand(1,p.soma_Npop);
D1_soma_somaMSNiSYN_s = zeros(nsamp,p.soma_Npop);
D1_soma_somaMSNiSYN_s(1,:) = D1_soma_somaMSNiSYN_s_last;
D1_D1_gabaRecInputMSN_s_last = zeros(1,p.D1_Npop);
D1_D1_gabaRecInputMSN_s = zeros(nsamp,p.D1_Npop);
D1_D1_gabaRecInputMSN_s(1,:) = D1_D1_gabaRecInputMSN_s_last;
D2_soma_somaMSNiSYN_s_last =  p.D2_soma_somaMSNiSYN_IC+p.D2_soma_somaMSNiSYN_IC_noise.*rand(1,p.soma_Npop);
D2_soma_somaMSNiSYN_s = zeros(nsamp,p.soma_Npop);
D2_soma_somaMSNiSYN_s(1,:) = D2_soma_somaMSNiSYN_s_last;
D2_D2_gabaRecInputMSN_s_last = zeros(1,p.D2_Npop);
D2_D2_gabaRecInputMSN_s = zeros(nsamp,p.D2_Npop);
D2_D2_gabaRecInputMSN_s(1,:) = D2_D2_gabaRecInputMSN_s_last;


% ###########################################################
% Memory check:
% ###########################################################
try 
  memoryUsed = memoryUsageCallerGB(); 
  fprintf('Total Memory Used <= %i GB \n', ceil(memoryUsed)); 
end 


% ###########################################################
% Numerical integration:
% ###########################################################
% seed the random number generator
rng_wrapper(p.random_seed);
n=2;
for k=2:ntime
  t=T(k-1);
  soma_V_k1 = (p.soma_Iapp + ((-(( soma_somaGolombK_gd2.*(soma_somaGolombK_a_last.^3).*soma_somaGolombK_b_last.*(soma_V_last-p.soma_somaGolombK_vk))))+((-(( p.soma_somaGolombKdr_gkdr.*(soma_somaGolombKdr_n_last.^p.soma_somaGolombKdr_nexp).*(soma_V_last-p.soma_somaGolombKdr_vk))))+(((( p.soma_somaInput_soma_tonic)))+((-(( p.soma_somaGolombNa_gna.*soma_somaGolombNa_h_last.*((( 1./(1+exp(-(soma_V_last-p.soma_somaGolombNa_thetam)./p.soma_somaGolombNa_sigmam)))).^3).*(soma_V_last-p.soma_somaGolombNa_vna))))+((-(( soma_somaLeak_gl2.*(soma_V_last-p.soma_somaLeak_vl))))+((-(( (soma_soma_somaSomaiSYN_gsyn2.*(soma_soma_somaSomaiSYN_s_last*soma_soma_somaSomaiSYN_mask).*(soma_V_last-p.soma_soma_somaSomaiSYN_Esyn)))))+(((( p.soma_dend_dendSomaiCOM_gCOM.*(dend_V_last-soma_V_last))))))))))) )/p.soma_Cm;
  soma_somaGolombK_a_k1 = ((( 1./(1+exp(-(soma_V_last-p.soma_somaGolombK_theta_a)./p.soma_somaGolombK_sigma_a))))-soma_somaGolombK_a_last)./p.soma_somaGolombK_tau_a;
  soma_somaGolombK_b_k1 = ((( 1./(1+exp(-(soma_V_last-p.soma_somaGolombK_theta_b)./p.soma_somaGolombK_sigma_b))))-soma_somaGolombK_b_last)./p.soma_somaGolombK_taub;
  soma_somaGolombKdr_n_k1 = ((( 1./(1+exp(-(soma_V_last-p.soma_somaGolombKdr_thetan)./p.soma_somaGolombKdr_sigman))))-soma_somaGolombKdr_n_last)./(( p.soma_somaGolombKdr_tau_mult*(p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp((soma_V_last-p.soma_somaGolombKdr_phin1)./p.soma_somaGolombKdr_sigman1))) .* (p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp(-(soma_V_last-p.soma_somaGolombKdr_phin2)./p.soma_somaGolombKdr_sigman2)))));
  soma_somaGolombNa_h_k1 = ((( 1./(1+exp(-(soma_V_last-p.soma_somaGolombNa_thetah)./p.soma_somaGolombNa_sigmah))))-soma_somaGolombNa_h_last)./(( p.soma_somaGolombNa_tauminh + p.soma_somaGolombNa_taumaxh./(1+exp(-(soma_V_last-p.soma_somaGolombNa_phih)./p.soma_somaGolombNa_sigmath))));
  dend_V_k1 = (p.dend_Iapp + ((-(( dend_dendGolombK_gd_dend.*(dend_dendGolombK_a_last.^3).*dend_dendGolombK_b_last.*(dend_V_last-p.dend_dendGolombK_vk))))+((-(( dend_dendGolombKdr_gkdr_dend.*(dend_dendGolombKdr_n_last.^p.dend_dendGolombKdr_nexp).*(dend_V_last-p.dend_dendGolombKdr_vk))))+((-(( dend_dendGolombNa_gna_dend.*dend_dendGolombNa_h_last.*((( 1./(1+exp(-(dend_V_last-p.dend_dendGolombNa_thetam)./p.dend_dendGolombNa_sigmam)))).^3).*(dend_V_last-p.dend_dendGolombNa_vna))))+(((( dend_dendInput_tonic2 + p.dend_dendInput_fs_noise.*randn(1,p.dend_Npop))))+((-(( dend_dendLeak_gl_dend.*(dend_V_last-p.dend_dendLeak_vl))))+((-(( (( (( p.dend_dendiMultiPoissonExp_g_esyn.*dend_dendiMultiPoissonExp_Ge(:, max(1,round(t/dt)))')).*(dend_V_last - p.dend_dendiMultiPoissonExp_E_esyn))) + (( (( p.dend_dendiMultiPoissonExp_g_isyn.*dend_dendiMultiPoissonExp_Gi(:, max(1,round(t/dt)))')).*(dend_V_last - p.dend_dendiMultiPoissonExp_E_isyn))))))+(((( p.dend_soma_somaDendiCOM_gCOM.*(soma_V_last-dend_V_last))))+(((( dend_dend_dendDendiGAP_g_GAP2.*sum((( ((dend_V_last'*ones(1,size(dend_V_last,2)))'-(dend_V_last'*ones(1,size(dend_V_last,2))))')).*dend_dend_dendDendiGAP_mask,1)))))))))))) )/p.dend_Cm;
  dend_dendGolombK_a_k1 = ((( 1./(1+exp(-(dend_V_last-p.dend_dendGolombK_theta_a)./p.dend_dendGolombK_sigma_a))))-dend_dendGolombK_a_last)./p.dend_dendGolombK_tau_a;
  dend_dendGolombK_b_k1 = ((( 1./(1+exp(-(dend_V_last-p.dend_dendGolombK_theta_b)./p.dend_dendGolombK_sigma_b))))-dend_dendGolombK_b_last)./p.dend_dendGolombK_taub;
  dend_dendGolombKdr_n_k1 = ((( 1./(1+exp(-(dend_V_last-p.dend_dendGolombKdr_thetan)./p.dend_dendGolombKdr_sigman))))-dend_dendGolombKdr_n_last)./(( p.dend_dendGolombKdr_tau_mult*(p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp((dend_V_last-p.dend_dendGolombKdr_phin1)./p.dend_dendGolombKdr_sigman1))) .* (p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp(-(dend_V_last-p.dend_dendGolombKdr_phin2)./p.dend_dendGolombKdr_sigman2)))));
  dend_dendGolombNa_h_k1 = ((( 1./(1+exp(-(dend_V_last-p.dend_dendGolombNa_thetah)./p.dend_dendGolombNa_sigmah))))-dend_dendGolombNa_h_last)./(( p.dend_dendGolombNa_tauminh + p.dend_dendGolombNa_taumaxh./(1+exp(-(dend_V_last-p.dend_dendGolombNa_phih)./p.dend_dendGolombNa_sigmath))));
  D1_V_k1 = (p.D1_Iapp + ((-(( p.D1_naCurrentMSN_g_na*D1_naCurrentMSN_m_last.^3.*D1_naCurrentMSN_h_last.*(D1_V_last-p.D1_naCurrentMSN_E_na))))+((-((p.D1_kCurrentMSN_g_k*D1_kCurrentMSN_m_last.^4.*(D1_V_last-p.D1_kCurrentMSN_E_k))))+((-(( p.D1_mCurrentMSN_g_m*D1_mCurrentMSN_m_last.*(D1_V_last-p.D1_mCurrentMSN_E_m))))+((-((p.D1_leakCurrentMSN_g_l*(D1_V_last-p.D1_leakCurrentMSN_E_l))))+(((( p.D1_injectedCurrentD1_injectedCurrent + p.D1_injectedCurrentD1_sinmult*sin(2*pi*D1_injectedCurrentD1_freq*t/1000) + p.D1_injectedCurrentD1_DAmult*p.D1_injectedCurrentD1_DA)))+(((( p.D1_noisyInputMSN_sigma_noise.*randn(1,p.D1_Npop).*sqrt(dt))))+(((( D1_AMPAMSN_cellmask1.*D1_AMPAMSN_psptime(round(t)+1) + D1_AMPAMSN_cellmask2.*D1_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D1_soma_somaMSNiSYN_m_gsyn.*(D1_soma_somaMSNiSYN_s_last*D1_soma_somaMSNiSYN_mask).*(D1_V_last-p.D1_soma_somaMSNiSYN_Esyn)))))+((-(( p.D1_D1_gabaRecInputMSN_g_gaba.*(D1_D1_gabaRecInputMSN_s_last*D1_D1_gabaRecInputMSN_netcon).*(D1_V_last-p.D1_D1_gabaRecInputMSN_E_gaba))))))))))))) )/p.D1_Cm;
  D1_naCurrentMSN_m_k1 = (( 0.32*(D1_V_last+54)./(1-exp(-(D1_V_last+54)/4)))).*(1-D1_naCurrentMSN_m_last)-(( 0.28*(D1_V_last+27)./(exp((D1_V_last+27)/5)-1))).*D1_naCurrentMSN_m_last;
  D1_naCurrentMSN_h_k1 = (( 0.128*exp(-(D1_V_last+50)/18))).*(1-D1_naCurrentMSN_h_last)-(( 4./(1+exp(-(D1_V_last+27)/5)))).*D1_naCurrentMSN_h_last;
  D1_kCurrentMSN_m_k1 = (( 0.032*(D1_V_last+52)./(1-exp(-(D1_V_last+52)/5)))).*(1-D1_kCurrentMSN_m_last)-(( 0.5*exp(-(D1_V_last+57)/40))).*D1_kCurrentMSN_m_last;
  D1_mCurrentMSN_m_k1 = (( D1_mCurrentMSN_Qs*1e-4*(D1_V_last-p.D1_mCurrentMSN_vhalf)./(1-exp(-(D1_V_last-p.D1_mCurrentMSN_vhalf)/9)))).*(1-D1_mCurrentMSN_m_last) - (( -D1_mCurrentMSN_Qs*1e-4*(D1_V_last-p.D1_mCurrentMSN_vhalf)./(1-exp((D1_V_last-p.D1_mCurrentMSN_vhalf)/9)))).*D1_mCurrentMSN_m_last;
  D2_V_k1 = (p.D2_Iapp + ((-(( p.D2_naCurrentMSN_g_na*D2_naCurrentMSN_m_last.^3.*D2_naCurrentMSN_h_last.*(D2_V_last-p.D2_naCurrentMSN_E_na))))+((-((p.D2_kCurrentMSN_g_k*D2_kCurrentMSN_m_last.^4.*(D2_V_last-p.D2_kCurrentMSN_E_k))))+((-(( p.D2_mCurrentMSN_g_m*D2_mCurrentMSN_m_last.*(D2_V_last-p.D2_mCurrentMSN_E_m))))+((-((p.D2_leakCurrentMSN_g_l*(D2_V_last-p.D2_leakCurrentMSN_E_l))))+(((( p.D2_injectedCurrentD2_injectedCurrent + p.D2_injectedCurrentD2_sinmult*sin(2*pi*D2_injectedCurrentD2_freq*t/1000) - p.D2_injectedCurrentD2_DAmult*p.D2_injectedCurrentD2_DA)))+(((( p.D2_noisyInputMSN_sigma_noise.*randn(1,p.D2_Npop).*sqrt(dt))))+(((( D2_AMPAMSN_cellmask1.*D2_AMPAMSN_psptime(round(t)+1) + D2_AMPAMSN_cellmask2.*D2_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D2_soma_somaMSNiSYN_m_gsyn.*(D2_soma_somaMSNiSYN_s_last*D2_soma_somaMSNiSYN_mask).*(D2_V_last-p.D2_soma_somaMSNiSYN_Esyn)))))+((-(( p.D2_D2_gabaRecInputMSN_g_gaba.*(D2_D2_gabaRecInputMSN_s_last*D2_D2_gabaRecInputMSN_netcon).*(D2_V_last-p.D2_D2_gabaRecInputMSN_E_gaba))))))))))))) )/p.D2_Cm;
  D2_naCurrentMSN_m_k1 = (( 0.32*(D2_V_last+54)./(1-exp(-(D2_V_last+54)/4)))).*(1-D2_naCurrentMSN_m_last)-(( 0.28*(D2_V_last+27)./(exp((D2_V_last+27)/5)-1))).*D2_naCurrentMSN_m_last;
  D2_naCurrentMSN_h_k1 = (( 0.128*exp(-(D2_V_last+50)/18))).*(1-D2_naCurrentMSN_h_last)-(( 4./(1+exp(-(D2_V_last+27)/5)))).*D2_naCurrentMSN_h_last;
  D2_kCurrentMSN_m_k1 = (( 0.032*(D2_V_last+52)./(1-exp(-(D2_V_last+52)/5)))).*(1-D2_kCurrentMSN_m_last)-(( 0.5*exp(-(D2_V_last+57)/40))).*D2_kCurrentMSN_m_last;
  D2_mCurrentMSN_m_k1 = (( D2_mCurrentMSN_Qs*1e-4*(D2_V_last-p.D2_mCurrentMSN_vhalf)./(1-exp(-(D2_V_last-p.D2_mCurrentMSN_vhalf)/9)))).*(1-D2_mCurrentMSN_m_last) - (( -D2_mCurrentMSN_Qs*1e-4*(D2_V_last-p.D2_mCurrentMSN_vhalf)./(1-exp((D2_V_last-p.D2_mCurrentMSN_vhalf)/9)))).*D2_mCurrentMSN_m_last;
  soma_soma_somaSomaiSYN_s_k1 = -soma_soma_somaSomaiSYN_s_last./p.soma_soma_somaSomaiSYN_tauD + ((1-soma_soma_somaSomaiSYN_s_last)/p.soma_soma_somaSomaiSYN_tauR).*(1+tanh(soma_V_last/10));
  D1_soma_somaMSNiSYN_s_k1 = -D1_soma_somaMSNiSYN_s_last./p.D1_soma_somaMSNiSYN_tauD + ((1-D1_soma_somaMSNiSYN_s_last)/p.D1_soma_somaMSNiSYN_tauR).*(1+tanh(soma_V_last/10));
  D1_D1_gabaRecInputMSN_s_k1 = -D1_D1_gabaRecInputMSN_s_last./p.D1_D1_gabaRecInputMSN_tau_gaba + 2*(1+tanh(D1_V_last/4)).*(1-D1_D1_gabaRecInputMSN_s_last);
  D2_soma_somaMSNiSYN_s_k1 = -D2_soma_somaMSNiSYN_s_last./p.D2_soma_somaMSNiSYN_tauD + ((1-D2_soma_somaMSNiSYN_s_last)/p.D2_soma_somaMSNiSYN_tauR).*(1+tanh(soma_V_last/10));
  D2_D2_gabaRecInputMSN_s_k1 = -D2_D2_gabaRecInputMSN_s_last./p.D2_D2_gabaRecInputMSN_tau_gaba + 2*(1+tanh(D2_V_last/4)).*(1-D2_D2_gabaRecInputMSN_s_last);

  t = t + .5*dt;
  soma_V_k2 = (p.soma_Iapp + ((-(( soma_somaGolombK_gd2.*(((soma_somaGolombK_a_last + .5*dt*soma_somaGolombK_a_k1)).^3).*((soma_somaGolombK_b_last + .5*dt*soma_somaGolombK_b_k1)).*(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombK_vk))))+((-(( p.soma_somaGolombKdr_gkdr.*(((soma_somaGolombKdr_n_last + .5*dt*soma_somaGolombKdr_n_k1)).^p.soma_somaGolombKdr_nexp).*(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombKdr_vk))))+(((( p.soma_somaInput_soma_tonic)))+((-(( p.soma_somaGolombNa_gna.*((soma_somaGolombNa_h_last + .5*dt*soma_somaGolombNa_h_k1)).*((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombNa_thetam)./p.soma_somaGolombNa_sigmam)))).^3).*(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombNa_vna))))+((-(( soma_somaLeak_gl2.*(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaLeak_vl))))+((-(( (soma_soma_somaSomaiSYN_gsyn2.*(((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k1))*soma_soma_somaSomaiSYN_mask).*(((soma_V_last + .5*dt*soma_V_k1))-p.soma_soma_somaSomaiSYN_Esyn)))))+(((( p.soma_dend_dendSomaiCOM_gCOM.*(((dend_V_last + .5*dt*dend_V_k1))-((soma_V_last + .5*dt*soma_V_k1))))))))))))) )/p.soma_Cm;
  soma_somaGolombK_a_k2 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombK_theta_a)./p.soma_somaGolombK_sigma_a))))-((soma_somaGolombK_a_last + .5*dt*soma_somaGolombK_a_k1)))./p.soma_somaGolombK_tau_a;
  soma_somaGolombK_b_k2 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombK_theta_b)./p.soma_somaGolombK_sigma_b))))-((soma_somaGolombK_b_last + .5*dt*soma_somaGolombK_b_k1)))./p.soma_somaGolombK_taub;
  soma_somaGolombKdr_n_k2 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombKdr_thetan)./p.soma_somaGolombKdr_sigman))))-((soma_somaGolombKdr_n_last + .5*dt*soma_somaGolombKdr_n_k1)))./(( p.soma_somaGolombKdr_tau_mult*(p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp((((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombKdr_phin1)./p.soma_somaGolombKdr_sigman1))) .* (p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombKdr_phin2)./p.soma_somaGolombKdr_sigman2)))));
  soma_somaGolombNa_h_k2 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombNa_thetah)./p.soma_somaGolombNa_sigmah))))-((soma_somaGolombNa_h_last + .5*dt*soma_somaGolombNa_h_k1)))./(( p.soma_somaGolombNa_tauminh + p.soma_somaGolombNa_taumaxh./(1+exp(-(((soma_V_last + .5*dt*soma_V_k1))-p.soma_somaGolombNa_phih)./p.soma_somaGolombNa_sigmath))));
  dend_V_k2 = (p.dend_Iapp + ((-(( dend_dendGolombK_gd_dend.*(((dend_dendGolombK_a_last + .5*dt*dend_dendGolombK_a_k1)).^3).*((dend_dendGolombK_b_last + .5*dt*dend_dendGolombK_b_k1)).*(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombK_vk))))+((-(( dend_dendGolombKdr_gkdr_dend.*(((dend_dendGolombKdr_n_last + .5*dt*dend_dendGolombKdr_n_k1)).^p.dend_dendGolombKdr_nexp).*(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombKdr_vk))))+((-(( dend_dendGolombNa_gna_dend.*((dend_dendGolombNa_h_last + .5*dt*dend_dendGolombNa_h_k1)).*((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombNa_thetam)./p.dend_dendGolombNa_sigmam)))).^3).*(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombNa_vna))))+(((( dend_dendInput_tonic2 + p.dend_dendInput_fs_noise.*randn(1,p.dend_Npop))))+((-(( dend_dendLeak_gl_dend.*(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendLeak_vl))))+((-(( (( (( p.dend_dendiMultiPoissonExp_g_esyn.*dend_dendiMultiPoissonExp_Ge(:, max(1,round(t/dt)))')).*(((dend_V_last + .5*dt*dend_V_k1)) - p.dend_dendiMultiPoissonExp_E_esyn))) + (( (( p.dend_dendiMultiPoissonExp_g_isyn.*dend_dendiMultiPoissonExp_Gi(:, max(1,round(t/dt)))')).*(((dend_V_last + .5*dt*dend_V_k1)) - p.dend_dendiMultiPoissonExp_E_isyn))))))+(((( p.dend_soma_somaDendiCOM_gCOM.*(((soma_V_last + .5*dt*soma_V_k1))-((dend_V_last + .5*dt*dend_V_k1))))))+(((( dend_dend_dendDendiGAP_g_GAP2.*sum((( ((((dend_V_last + .5*dt*dend_V_k1))'*ones(1,size(((dend_V_last + .5*dt*dend_V_k1)),2)))'-(((dend_V_last + .5*dt*dend_V_k1))'*ones(1,size(((dend_V_last + .5*dt*dend_V_k1)),2))))')).*dend_dend_dendDendiGAP_mask,1)))))))))))) )/p.dend_Cm;
  dend_dendGolombK_a_k2 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombK_theta_a)./p.dend_dendGolombK_sigma_a))))-((dend_dendGolombK_a_last + .5*dt*dend_dendGolombK_a_k1)))./p.dend_dendGolombK_tau_a;
  dend_dendGolombK_b_k2 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombK_theta_b)./p.dend_dendGolombK_sigma_b))))-((dend_dendGolombK_b_last + .5*dt*dend_dendGolombK_b_k1)))./p.dend_dendGolombK_taub;
  dend_dendGolombKdr_n_k2 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombKdr_thetan)./p.dend_dendGolombKdr_sigman))))-((dend_dendGolombKdr_n_last + .5*dt*dend_dendGolombKdr_n_k1)))./(( p.dend_dendGolombKdr_tau_mult*(p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp((((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombKdr_phin1)./p.dend_dendGolombKdr_sigman1))) .* (p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombKdr_phin2)./p.dend_dendGolombKdr_sigman2)))));
  dend_dendGolombNa_h_k2 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombNa_thetah)./p.dend_dendGolombNa_sigmah))))-((dend_dendGolombNa_h_last + .5*dt*dend_dendGolombNa_h_k1)))./(( p.dend_dendGolombNa_tauminh + p.dend_dendGolombNa_taumaxh./(1+exp(-(((dend_V_last + .5*dt*dend_V_k1))-p.dend_dendGolombNa_phih)./p.dend_dendGolombNa_sigmath))));
  D1_V_k2 = (p.D1_Iapp + ((-(( p.D1_naCurrentMSN_g_na*((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k1)).^3.*((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k1)).*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_naCurrentMSN_E_na))))+((-((p.D1_kCurrentMSN_g_k*((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k1)).^4.*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_kCurrentMSN_E_k))))+((-(( p.D1_mCurrentMSN_g_m*((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k1)).*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_mCurrentMSN_E_m))))+((-((p.D1_leakCurrentMSN_g_l*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_leakCurrentMSN_E_l))))+(((( p.D1_injectedCurrentD1_injectedCurrent + p.D1_injectedCurrentD1_sinmult*sin(2*pi*D1_injectedCurrentD1_freq*t/1000) + p.D1_injectedCurrentD1_DAmult*p.D1_injectedCurrentD1_DA)))+(((( p.D1_noisyInputMSN_sigma_noise.*randn(1,p.D1_Npop).*sqrt(dt))))+(((( D1_AMPAMSN_cellmask1.*D1_AMPAMSN_psptime(round(t)+1) + D1_AMPAMSN_cellmask2.*D1_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D1_soma_somaMSNiSYN_m_gsyn.*(((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k1))*D1_soma_somaMSNiSYN_mask).*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_soma_somaMSNiSYN_Esyn)))))+((-(( p.D1_D1_gabaRecInputMSN_g_gaba.*(((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k1))*D1_D1_gabaRecInputMSN_netcon).*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_D1_gabaRecInputMSN_E_gaba))))))))))))) )/p.D1_Cm;
  D1_naCurrentMSN_m_k2 = (( 0.32*(((D1_V_last + .5*dt*D1_V_k1))+54)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k1))+54)/4)))).*(1-((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k1)))-(( 0.28*(((D1_V_last + .5*dt*D1_V_k1))+27)./(exp((((D1_V_last + .5*dt*D1_V_k1))+27)/5)-1))).*((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k1));
  D1_naCurrentMSN_h_k2 = (( 0.128*exp(-(((D1_V_last + .5*dt*D1_V_k1))+50)/18))).*(1-((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k1)))-(( 4./(1+exp(-(((D1_V_last + .5*dt*D1_V_k1))+27)/5)))).*((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k1));
  D1_kCurrentMSN_m_k2 = (( 0.032*(((D1_V_last + .5*dt*D1_V_k1))+52)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k1))+52)/5)))).*(1-((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k1)))-(( 0.5*exp(-(((D1_V_last + .5*dt*D1_V_k1))+57)/40))).*((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k1));
  D1_mCurrentMSN_m_k2 = (( D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_mCurrentMSN_vhalf)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k1))-p.D1_mCurrentMSN_vhalf)/9)))).*(1-((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k1))) - (( -D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + .5*dt*D1_V_k1))-p.D1_mCurrentMSN_vhalf)./(1-exp((((D1_V_last + .5*dt*D1_V_k1))-p.D1_mCurrentMSN_vhalf)/9)))).*((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k1));
  D2_V_k2 = (p.D2_Iapp + ((-(( p.D2_naCurrentMSN_g_na*((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k1)).^3.*((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k1)).*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_naCurrentMSN_E_na))))+((-((p.D2_kCurrentMSN_g_k*((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k1)).^4.*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_kCurrentMSN_E_k))))+((-(( p.D2_mCurrentMSN_g_m*((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k1)).*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_mCurrentMSN_E_m))))+((-((p.D2_leakCurrentMSN_g_l*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_leakCurrentMSN_E_l))))+(((( p.D2_injectedCurrentD2_injectedCurrent + p.D2_injectedCurrentD2_sinmult*sin(2*pi*D2_injectedCurrentD2_freq*t/1000) - p.D2_injectedCurrentD2_DAmult*p.D2_injectedCurrentD2_DA)))+(((( p.D2_noisyInputMSN_sigma_noise.*randn(1,p.D2_Npop).*sqrt(dt))))+(((( D2_AMPAMSN_cellmask1.*D2_AMPAMSN_psptime(round(t)+1) + D2_AMPAMSN_cellmask2.*D2_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D2_soma_somaMSNiSYN_m_gsyn.*(((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k1))*D2_soma_somaMSNiSYN_mask).*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_soma_somaMSNiSYN_Esyn)))))+((-(( p.D2_D2_gabaRecInputMSN_g_gaba.*(((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k1))*D2_D2_gabaRecInputMSN_netcon).*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_D2_gabaRecInputMSN_E_gaba))))))))))))) )/p.D2_Cm;
  D2_naCurrentMSN_m_k2 = (( 0.32*(((D2_V_last + .5*dt*D2_V_k1))+54)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k1))+54)/4)))).*(1-((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k1)))-(( 0.28*(((D2_V_last + .5*dt*D2_V_k1))+27)./(exp((((D2_V_last + .5*dt*D2_V_k1))+27)/5)-1))).*((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k1));
  D2_naCurrentMSN_h_k2 = (( 0.128*exp(-(((D2_V_last + .5*dt*D2_V_k1))+50)/18))).*(1-((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k1)))-(( 4./(1+exp(-(((D2_V_last + .5*dt*D2_V_k1))+27)/5)))).*((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k1));
  D2_kCurrentMSN_m_k2 = (( 0.032*(((D2_V_last + .5*dt*D2_V_k1))+52)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k1))+52)/5)))).*(1-((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k1)))-(( 0.5*exp(-(((D2_V_last + .5*dt*D2_V_k1))+57)/40))).*((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k1));
  D2_mCurrentMSN_m_k2 = (( D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_mCurrentMSN_vhalf)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k1))-p.D2_mCurrentMSN_vhalf)/9)))).*(1-((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k1))) - (( -D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + .5*dt*D2_V_k1))-p.D2_mCurrentMSN_vhalf)./(1-exp((((D2_V_last + .5*dt*D2_V_k1))-p.D2_mCurrentMSN_vhalf)/9)))).*((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k1));
  soma_soma_somaSomaiSYN_s_k2 = -((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k1))./p.soma_soma_somaSomaiSYN_tauD + ((1-((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k1)))/p.soma_soma_somaSomaiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k1))/10));
  D1_soma_somaMSNiSYN_s_k2 = -((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k1))./p.D1_soma_somaMSNiSYN_tauD + ((1-((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k1)))/p.D1_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k1))/10));
  D1_D1_gabaRecInputMSN_s_k2 = -((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k1))./p.D1_D1_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D1_V_last + .5*dt*D1_V_k1))/4)).*(1-((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k1)));
  D2_soma_somaMSNiSYN_s_k2 = -((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k1))./p.D2_soma_somaMSNiSYN_tauD + ((1-((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k1)))/p.D2_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k1))/10));
  D2_D2_gabaRecInputMSN_s_k2 = -((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k1))./p.D2_D2_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D2_V_last + .5*dt*D2_V_k1))/4)).*(1-((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k1)));

  soma_V_k3 = (p.soma_Iapp + ((-(( soma_somaGolombK_gd2.*(((soma_somaGolombK_a_last + .5*dt*soma_somaGolombK_a_k2)).^3).*((soma_somaGolombK_b_last + .5*dt*soma_somaGolombK_b_k2)).*(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombK_vk))))+((-(( p.soma_somaGolombKdr_gkdr.*(((soma_somaGolombKdr_n_last + .5*dt*soma_somaGolombKdr_n_k2)).^p.soma_somaGolombKdr_nexp).*(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombKdr_vk))))+(((( p.soma_somaInput_soma_tonic)))+((-(( p.soma_somaGolombNa_gna.*((soma_somaGolombNa_h_last + .5*dt*soma_somaGolombNa_h_k2)).*((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombNa_thetam)./p.soma_somaGolombNa_sigmam)))).^3).*(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombNa_vna))))+((-(( soma_somaLeak_gl2.*(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaLeak_vl))))+((-(( (soma_soma_somaSomaiSYN_gsyn2.*(((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k2))*soma_soma_somaSomaiSYN_mask).*(((soma_V_last + .5*dt*soma_V_k2))-p.soma_soma_somaSomaiSYN_Esyn)))))+(((( p.soma_dend_dendSomaiCOM_gCOM.*(((dend_V_last + .5*dt*dend_V_k2))-((soma_V_last + .5*dt*soma_V_k2))))))))))))) )/p.soma_Cm;
  soma_somaGolombK_a_k3 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombK_theta_a)./p.soma_somaGolombK_sigma_a))))-((soma_somaGolombK_a_last + .5*dt*soma_somaGolombK_a_k2)))./p.soma_somaGolombK_tau_a;
  soma_somaGolombK_b_k3 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombK_theta_b)./p.soma_somaGolombK_sigma_b))))-((soma_somaGolombK_b_last + .5*dt*soma_somaGolombK_b_k2)))./p.soma_somaGolombK_taub;
  soma_somaGolombKdr_n_k3 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombKdr_thetan)./p.soma_somaGolombKdr_sigman))))-((soma_somaGolombKdr_n_last + .5*dt*soma_somaGolombKdr_n_k2)))./(( p.soma_somaGolombKdr_tau_mult*(p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp((((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombKdr_phin1)./p.soma_somaGolombKdr_sigman1))) .* (p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombKdr_phin2)./p.soma_somaGolombKdr_sigman2)))));
  soma_somaGolombNa_h_k3 = ((( 1./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombNa_thetah)./p.soma_somaGolombNa_sigmah))))-((soma_somaGolombNa_h_last + .5*dt*soma_somaGolombNa_h_k2)))./(( p.soma_somaGolombNa_tauminh + p.soma_somaGolombNa_taumaxh./(1+exp(-(((soma_V_last + .5*dt*soma_V_k2))-p.soma_somaGolombNa_phih)./p.soma_somaGolombNa_sigmath))));
  dend_V_k3 = (p.dend_Iapp + ((-(( dend_dendGolombK_gd_dend.*(((dend_dendGolombK_a_last + .5*dt*dend_dendGolombK_a_k2)).^3).*((dend_dendGolombK_b_last + .5*dt*dend_dendGolombK_b_k2)).*(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombK_vk))))+((-(( dend_dendGolombKdr_gkdr_dend.*(((dend_dendGolombKdr_n_last + .5*dt*dend_dendGolombKdr_n_k2)).^p.dend_dendGolombKdr_nexp).*(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombKdr_vk))))+((-(( dend_dendGolombNa_gna_dend.*((dend_dendGolombNa_h_last + .5*dt*dend_dendGolombNa_h_k2)).*((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombNa_thetam)./p.dend_dendGolombNa_sigmam)))).^3).*(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombNa_vna))))+(((( dend_dendInput_tonic2 + p.dend_dendInput_fs_noise.*randn(1,p.dend_Npop))))+((-(( dend_dendLeak_gl_dend.*(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendLeak_vl))))+((-(( (( (( p.dend_dendiMultiPoissonExp_g_esyn.*dend_dendiMultiPoissonExp_Ge(:, max(1,round(t/dt)))')).*(((dend_V_last + .5*dt*dend_V_k2)) - p.dend_dendiMultiPoissonExp_E_esyn))) + (( (( p.dend_dendiMultiPoissonExp_g_isyn.*dend_dendiMultiPoissonExp_Gi(:, max(1,round(t/dt)))')).*(((dend_V_last + .5*dt*dend_V_k2)) - p.dend_dendiMultiPoissonExp_E_isyn))))))+(((( p.dend_soma_somaDendiCOM_gCOM.*(((soma_V_last + .5*dt*soma_V_k2))-((dend_V_last + .5*dt*dend_V_k2))))))+(((( dend_dend_dendDendiGAP_g_GAP2.*sum((( ((((dend_V_last + .5*dt*dend_V_k2))'*ones(1,size(((dend_V_last + .5*dt*dend_V_k2)),2)))'-(((dend_V_last + .5*dt*dend_V_k2))'*ones(1,size(((dend_V_last + .5*dt*dend_V_k2)),2))))')).*dend_dend_dendDendiGAP_mask,1)))))))))))) )/p.dend_Cm;
  dend_dendGolombK_a_k3 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombK_theta_a)./p.dend_dendGolombK_sigma_a))))-((dend_dendGolombK_a_last + .5*dt*dend_dendGolombK_a_k2)))./p.dend_dendGolombK_tau_a;
  dend_dendGolombK_b_k3 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombK_theta_b)./p.dend_dendGolombK_sigma_b))))-((dend_dendGolombK_b_last + .5*dt*dend_dendGolombK_b_k2)))./p.dend_dendGolombK_taub;
  dend_dendGolombKdr_n_k3 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombKdr_thetan)./p.dend_dendGolombKdr_sigman))))-((dend_dendGolombKdr_n_last + .5*dt*dend_dendGolombKdr_n_k2)))./(( p.dend_dendGolombKdr_tau_mult*(p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp((((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombKdr_phin1)./p.dend_dendGolombKdr_sigman1))) .* (p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombKdr_phin2)./p.dend_dendGolombKdr_sigman2)))));
  dend_dendGolombNa_h_k3 = ((( 1./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombNa_thetah)./p.dend_dendGolombNa_sigmah))))-((dend_dendGolombNa_h_last + .5*dt*dend_dendGolombNa_h_k2)))./(( p.dend_dendGolombNa_tauminh + p.dend_dendGolombNa_taumaxh./(1+exp(-(((dend_V_last + .5*dt*dend_V_k2))-p.dend_dendGolombNa_phih)./p.dend_dendGolombNa_sigmath))));
  D1_V_k3 = (p.D1_Iapp + ((-(( p.D1_naCurrentMSN_g_na*((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k2)).^3.*((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k2)).*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_naCurrentMSN_E_na))))+((-((p.D1_kCurrentMSN_g_k*((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k2)).^4.*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_kCurrentMSN_E_k))))+((-(( p.D1_mCurrentMSN_g_m*((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k2)).*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_mCurrentMSN_E_m))))+((-((p.D1_leakCurrentMSN_g_l*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_leakCurrentMSN_E_l))))+(((( p.D1_injectedCurrentD1_injectedCurrent + p.D1_injectedCurrentD1_sinmult*sin(2*pi*D1_injectedCurrentD1_freq*t/1000) + p.D1_injectedCurrentD1_DAmult*p.D1_injectedCurrentD1_DA)))+(((( p.D1_noisyInputMSN_sigma_noise.*randn(1,p.D1_Npop).*sqrt(dt))))+(((( D1_AMPAMSN_cellmask1.*D1_AMPAMSN_psptime(round(t)+1) + D1_AMPAMSN_cellmask2.*D1_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D1_soma_somaMSNiSYN_m_gsyn.*(((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k2))*D1_soma_somaMSNiSYN_mask).*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_soma_somaMSNiSYN_Esyn)))))+((-(( p.D1_D1_gabaRecInputMSN_g_gaba.*(((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k2))*D1_D1_gabaRecInputMSN_netcon).*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_D1_gabaRecInputMSN_E_gaba))))))))))))) )/p.D1_Cm;
  D1_naCurrentMSN_m_k3 = (( 0.32*(((D1_V_last + .5*dt*D1_V_k2))+54)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k2))+54)/4)))).*(1-((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k2)))-(( 0.28*(((D1_V_last + .5*dt*D1_V_k2))+27)./(exp((((D1_V_last + .5*dt*D1_V_k2))+27)/5)-1))).*((D1_naCurrentMSN_m_last + .5*dt*D1_naCurrentMSN_m_k2));
  D1_naCurrentMSN_h_k3 = (( 0.128*exp(-(((D1_V_last + .5*dt*D1_V_k2))+50)/18))).*(1-((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k2)))-(( 4./(1+exp(-(((D1_V_last + .5*dt*D1_V_k2))+27)/5)))).*((D1_naCurrentMSN_h_last + .5*dt*D1_naCurrentMSN_h_k2));
  D1_kCurrentMSN_m_k3 = (( 0.032*(((D1_V_last + .5*dt*D1_V_k2))+52)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k2))+52)/5)))).*(1-((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k2)))-(( 0.5*exp(-(((D1_V_last + .5*dt*D1_V_k2))+57)/40))).*((D1_kCurrentMSN_m_last + .5*dt*D1_kCurrentMSN_m_k2));
  D1_mCurrentMSN_m_k3 = (( D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_mCurrentMSN_vhalf)./(1-exp(-(((D1_V_last + .5*dt*D1_V_k2))-p.D1_mCurrentMSN_vhalf)/9)))).*(1-((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k2))) - (( -D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + .5*dt*D1_V_k2))-p.D1_mCurrentMSN_vhalf)./(1-exp((((D1_V_last + .5*dt*D1_V_k2))-p.D1_mCurrentMSN_vhalf)/9)))).*((D1_mCurrentMSN_m_last + .5*dt*D1_mCurrentMSN_m_k2));
  D2_V_k3 = (p.D2_Iapp + ((-(( p.D2_naCurrentMSN_g_na*((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k2)).^3.*((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k2)).*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_naCurrentMSN_E_na))))+((-((p.D2_kCurrentMSN_g_k*((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k2)).^4.*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_kCurrentMSN_E_k))))+((-(( p.D2_mCurrentMSN_g_m*((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k2)).*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_mCurrentMSN_E_m))))+((-((p.D2_leakCurrentMSN_g_l*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_leakCurrentMSN_E_l))))+(((( p.D2_injectedCurrentD2_injectedCurrent + p.D2_injectedCurrentD2_sinmult*sin(2*pi*D2_injectedCurrentD2_freq*t/1000) - p.D2_injectedCurrentD2_DAmult*p.D2_injectedCurrentD2_DA)))+(((( p.D2_noisyInputMSN_sigma_noise.*randn(1,p.D2_Npop).*sqrt(dt))))+(((( D2_AMPAMSN_cellmask1.*D2_AMPAMSN_psptime(round(t)+1) + D2_AMPAMSN_cellmask2.*D2_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D2_soma_somaMSNiSYN_m_gsyn.*(((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k2))*D2_soma_somaMSNiSYN_mask).*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_soma_somaMSNiSYN_Esyn)))))+((-(( p.D2_D2_gabaRecInputMSN_g_gaba.*(((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k2))*D2_D2_gabaRecInputMSN_netcon).*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_D2_gabaRecInputMSN_E_gaba))))))))))))) )/p.D2_Cm;
  D2_naCurrentMSN_m_k3 = (( 0.32*(((D2_V_last + .5*dt*D2_V_k2))+54)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k2))+54)/4)))).*(1-((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k2)))-(( 0.28*(((D2_V_last + .5*dt*D2_V_k2))+27)./(exp((((D2_V_last + .5*dt*D2_V_k2))+27)/5)-1))).*((D2_naCurrentMSN_m_last + .5*dt*D2_naCurrentMSN_m_k2));
  D2_naCurrentMSN_h_k3 = (( 0.128*exp(-(((D2_V_last + .5*dt*D2_V_k2))+50)/18))).*(1-((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k2)))-(( 4./(1+exp(-(((D2_V_last + .5*dt*D2_V_k2))+27)/5)))).*((D2_naCurrentMSN_h_last + .5*dt*D2_naCurrentMSN_h_k2));
  D2_kCurrentMSN_m_k3 = (( 0.032*(((D2_V_last + .5*dt*D2_V_k2))+52)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k2))+52)/5)))).*(1-((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k2)))-(( 0.5*exp(-(((D2_V_last + .5*dt*D2_V_k2))+57)/40))).*((D2_kCurrentMSN_m_last + .5*dt*D2_kCurrentMSN_m_k2));
  D2_mCurrentMSN_m_k3 = (( D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_mCurrentMSN_vhalf)./(1-exp(-(((D2_V_last + .5*dt*D2_V_k2))-p.D2_mCurrentMSN_vhalf)/9)))).*(1-((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k2))) - (( -D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + .5*dt*D2_V_k2))-p.D2_mCurrentMSN_vhalf)./(1-exp((((D2_V_last + .5*dt*D2_V_k2))-p.D2_mCurrentMSN_vhalf)/9)))).*((D2_mCurrentMSN_m_last + .5*dt*D2_mCurrentMSN_m_k2));
  soma_soma_somaSomaiSYN_s_k3 = -((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k2))./p.soma_soma_somaSomaiSYN_tauD + ((1-((soma_soma_somaSomaiSYN_s_last + .5*dt*soma_soma_somaSomaiSYN_s_k2)))/p.soma_soma_somaSomaiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k2))/10));
  D1_soma_somaMSNiSYN_s_k3 = -((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k2))./p.D1_soma_somaMSNiSYN_tauD + ((1-((D1_soma_somaMSNiSYN_s_last + .5*dt*D1_soma_somaMSNiSYN_s_k2)))/p.D1_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k2))/10));
  D1_D1_gabaRecInputMSN_s_k3 = -((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k2))./p.D1_D1_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D1_V_last + .5*dt*D1_V_k2))/4)).*(1-((D1_D1_gabaRecInputMSN_s_last + .5*dt*D1_D1_gabaRecInputMSN_s_k2)));
  D2_soma_somaMSNiSYN_s_k3 = -((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k2))./p.D2_soma_somaMSNiSYN_tauD + ((1-((D2_soma_somaMSNiSYN_s_last + .5*dt*D2_soma_somaMSNiSYN_s_k2)))/p.D2_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + .5*dt*soma_V_k2))/10));
  D2_D2_gabaRecInputMSN_s_k3 = -((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k2))./p.D2_D2_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D2_V_last + .5*dt*D2_V_k2))/4)).*(1-((D2_D2_gabaRecInputMSN_s_last + .5*dt*D2_D2_gabaRecInputMSN_s_k2)));

  t = t + .5*dt;
  soma_V_k4 = (p.soma_Iapp + ((-(( soma_somaGolombK_gd2.*(((soma_somaGolombK_a_last + dt*soma_somaGolombK_a_k3)).^3).*((soma_somaGolombK_b_last + dt*soma_somaGolombK_b_k3)).*(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombK_vk))))+((-(( p.soma_somaGolombKdr_gkdr.*(((soma_somaGolombKdr_n_last + dt*soma_somaGolombKdr_n_k3)).^p.soma_somaGolombKdr_nexp).*(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombKdr_vk))))+(((( p.soma_somaInput_soma_tonic)))+((-(( p.soma_somaGolombNa_gna.*((soma_somaGolombNa_h_last + dt*soma_somaGolombNa_h_k3)).*((( 1./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombNa_thetam)./p.soma_somaGolombNa_sigmam)))).^3).*(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombNa_vna))))+((-(( soma_somaLeak_gl2.*(((soma_V_last + dt*soma_V_k3))-p.soma_somaLeak_vl))))+((-(( (soma_soma_somaSomaiSYN_gsyn2.*(((soma_soma_somaSomaiSYN_s_last + dt*soma_soma_somaSomaiSYN_s_k3))*soma_soma_somaSomaiSYN_mask).*(((soma_V_last + dt*soma_V_k3))-p.soma_soma_somaSomaiSYN_Esyn)))))+(((( p.soma_dend_dendSomaiCOM_gCOM.*(((dend_V_last + dt*dend_V_k3))-((soma_V_last + dt*soma_V_k3))))))))))))) )/p.soma_Cm;
  soma_somaGolombK_a_k4 = ((( 1./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombK_theta_a)./p.soma_somaGolombK_sigma_a))))-((soma_somaGolombK_a_last + dt*soma_somaGolombK_a_k3)))./p.soma_somaGolombK_tau_a;
  soma_somaGolombK_b_k4 = ((( 1./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombK_theta_b)./p.soma_somaGolombK_sigma_b))))-((soma_somaGolombK_b_last + dt*soma_somaGolombK_b_k3)))./p.soma_somaGolombK_taub;
  soma_somaGolombKdr_n_k4 = ((( 1./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombKdr_thetan)./p.soma_somaGolombKdr_sigman))))-((soma_somaGolombKdr_n_last + dt*soma_somaGolombKdr_n_k3)))./(( p.soma_somaGolombKdr_tau_mult*(p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp((((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombKdr_phin1)./p.soma_somaGolombKdr_sigman1))) .* (p.soma_somaGolombKdr_tauminn+p.soma_somaGolombKdr_taumaxn./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombKdr_phin2)./p.soma_somaGolombKdr_sigman2)))));
  soma_somaGolombNa_h_k4 = ((( 1./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombNa_thetah)./p.soma_somaGolombNa_sigmah))))-((soma_somaGolombNa_h_last + dt*soma_somaGolombNa_h_k3)))./(( p.soma_somaGolombNa_tauminh + p.soma_somaGolombNa_taumaxh./(1+exp(-(((soma_V_last + dt*soma_V_k3))-p.soma_somaGolombNa_phih)./p.soma_somaGolombNa_sigmath))));
  dend_V_k4 = (p.dend_Iapp + ((-(( dend_dendGolombK_gd_dend.*(((dend_dendGolombK_a_last + dt*dend_dendGolombK_a_k3)).^3).*((dend_dendGolombK_b_last + dt*dend_dendGolombK_b_k3)).*(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombK_vk))))+((-(( dend_dendGolombKdr_gkdr_dend.*(((dend_dendGolombKdr_n_last + dt*dend_dendGolombKdr_n_k3)).^p.dend_dendGolombKdr_nexp).*(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombKdr_vk))))+((-(( dend_dendGolombNa_gna_dend.*((dend_dendGolombNa_h_last + dt*dend_dendGolombNa_h_k3)).*((( 1./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombNa_thetam)./p.dend_dendGolombNa_sigmam)))).^3).*(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombNa_vna))))+(((( dend_dendInput_tonic2 + p.dend_dendInput_fs_noise.*randn(1,p.dend_Npop))))+((-(( dend_dendLeak_gl_dend.*(((dend_V_last + dt*dend_V_k3))-p.dend_dendLeak_vl))))+((-(( (( (( p.dend_dendiMultiPoissonExp_g_esyn.*dend_dendiMultiPoissonExp_Ge(:, max(1,round(t/dt)))')).*(((dend_V_last + dt*dend_V_k3)) - p.dend_dendiMultiPoissonExp_E_esyn))) + (( (( p.dend_dendiMultiPoissonExp_g_isyn.*dend_dendiMultiPoissonExp_Gi(:, max(1,round(t/dt)))')).*(((dend_V_last + dt*dend_V_k3)) - p.dend_dendiMultiPoissonExp_E_isyn))))))+(((( p.dend_soma_somaDendiCOM_gCOM.*(((soma_V_last + dt*soma_V_k3))-((dend_V_last + dt*dend_V_k3))))))+(((( dend_dend_dendDendiGAP_g_GAP2.*sum((( ((((dend_V_last + dt*dend_V_k3))'*ones(1,size(((dend_V_last + dt*dend_V_k3)),2)))'-(((dend_V_last + dt*dend_V_k3))'*ones(1,size(((dend_V_last + dt*dend_V_k3)),2))))')).*dend_dend_dendDendiGAP_mask,1)))))))))))) )/p.dend_Cm;
  dend_dendGolombK_a_k4 = ((( 1./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombK_theta_a)./p.dend_dendGolombK_sigma_a))))-((dend_dendGolombK_a_last + dt*dend_dendGolombK_a_k3)))./p.dend_dendGolombK_tau_a;
  dend_dendGolombK_b_k4 = ((( 1./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombK_theta_b)./p.dend_dendGolombK_sigma_b))))-((dend_dendGolombK_b_last + dt*dend_dendGolombK_b_k3)))./p.dend_dendGolombK_taub;
  dend_dendGolombKdr_n_k4 = ((( 1./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombKdr_thetan)./p.dend_dendGolombKdr_sigman))))-((dend_dendGolombKdr_n_last + dt*dend_dendGolombKdr_n_k3)))./(( p.dend_dendGolombKdr_tau_mult*(p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp((((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombKdr_phin1)./p.dend_dendGolombKdr_sigman1))) .* (p.dend_dendGolombKdr_tauminn+p.dend_dendGolombKdr_taumaxn./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombKdr_phin2)./p.dend_dendGolombKdr_sigman2)))));
  dend_dendGolombNa_h_k4 = ((( 1./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombNa_thetah)./p.dend_dendGolombNa_sigmah))))-((dend_dendGolombNa_h_last + dt*dend_dendGolombNa_h_k3)))./(( p.dend_dendGolombNa_tauminh + p.dend_dendGolombNa_taumaxh./(1+exp(-(((dend_V_last + dt*dend_V_k3))-p.dend_dendGolombNa_phih)./p.dend_dendGolombNa_sigmath))));
  D1_V_k4 = (p.D1_Iapp + ((-(( p.D1_naCurrentMSN_g_na*((D1_naCurrentMSN_m_last + dt*D1_naCurrentMSN_m_k3)).^3.*((D1_naCurrentMSN_h_last + dt*D1_naCurrentMSN_h_k3)).*(((D1_V_last + dt*D1_V_k3))-p.D1_naCurrentMSN_E_na))))+((-((p.D1_kCurrentMSN_g_k*((D1_kCurrentMSN_m_last + dt*D1_kCurrentMSN_m_k3)).^4.*(((D1_V_last + dt*D1_V_k3))-p.D1_kCurrentMSN_E_k))))+((-(( p.D1_mCurrentMSN_g_m*((D1_mCurrentMSN_m_last + dt*D1_mCurrentMSN_m_k3)).*(((D1_V_last + dt*D1_V_k3))-p.D1_mCurrentMSN_E_m))))+((-((p.D1_leakCurrentMSN_g_l*(((D1_V_last + dt*D1_V_k3))-p.D1_leakCurrentMSN_E_l))))+(((( p.D1_injectedCurrentD1_injectedCurrent + p.D1_injectedCurrentD1_sinmult*sin(2*pi*D1_injectedCurrentD1_freq*t/1000) + p.D1_injectedCurrentD1_DAmult*p.D1_injectedCurrentD1_DA)))+(((( p.D1_noisyInputMSN_sigma_noise.*randn(1,p.D1_Npop).*sqrt(dt))))+(((( D1_AMPAMSN_cellmask1.*D1_AMPAMSN_psptime(round(t)+1) + D1_AMPAMSN_cellmask2.*D1_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D1_soma_somaMSNiSYN_m_gsyn.*(((D1_soma_somaMSNiSYN_s_last + dt*D1_soma_somaMSNiSYN_s_k3))*D1_soma_somaMSNiSYN_mask).*(((D1_V_last + dt*D1_V_k3))-p.D1_soma_somaMSNiSYN_Esyn)))))+((-(( p.D1_D1_gabaRecInputMSN_g_gaba.*(((D1_D1_gabaRecInputMSN_s_last + dt*D1_D1_gabaRecInputMSN_s_k3))*D1_D1_gabaRecInputMSN_netcon).*(((D1_V_last + dt*D1_V_k3))-p.D1_D1_gabaRecInputMSN_E_gaba))))))))))))) )/p.D1_Cm;
  D1_naCurrentMSN_m_k4 = (( 0.32*(((D1_V_last + dt*D1_V_k3))+54)./(1-exp(-(((D1_V_last + dt*D1_V_k3))+54)/4)))).*(1-((D1_naCurrentMSN_m_last + dt*D1_naCurrentMSN_m_k3)))-(( 0.28*(((D1_V_last + dt*D1_V_k3))+27)./(exp((((D1_V_last + dt*D1_V_k3))+27)/5)-1))).*((D1_naCurrentMSN_m_last + dt*D1_naCurrentMSN_m_k3));
  D1_naCurrentMSN_h_k4 = (( 0.128*exp(-(((D1_V_last + dt*D1_V_k3))+50)/18))).*(1-((D1_naCurrentMSN_h_last + dt*D1_naCurrentMSN_h_k3)))-(( 4./(1+exp(-(((D1_V_last + dt*D1_V_k3))+27)/5)))).*((D1_naCurrentMSN_h_last + dt*D1_naCurrentMSN_h_k3));
  D1_kCurrentMSN_m_k4 = (( 0.032*(((D1_V_last + dt*D1_V_k3))+52)./(1-exp(-(((D1_V_last + dt*D1_V_k3))+52)/5)))).*(1-((D1_kCurrentMSN_m_last + dt*D1_kCurrentMSN_m_k3)))-(( 0.5*exp(-(((D1_V_last + dt*D1_V_k3))+57)/40))).*((D1_kCurrentMSN_m_last + dt*D1_kCurrentMSN_m_k3));
  D1_mCurrentMSN_m_k4 = (( D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + dt*D1_V_k3))-p.D1_mCurrentMSN_vhalf)./(1-exp(-(((D1_V_last + dt*D1_V_k3))-p.D1_mCurrentMSN_vhalf)/9)))).*(1-((D1_mCurrentMSN_m_last + dt*D1_mCurrentMSN_m_k3))) - (( -D1_mCurrentMSN_Qs*1e-4*(((D1_V_last + dt*D1_V_k3))-p.D1_mCurrentMSN_vhalf)./(1-exp((((D1_V_last + dt*D1_V_k3))-p.D1_mCurrentMSN_vhalf)/9)))).*((D1_mCurrentMSN_m_last + dt*D1_mCurrentMSN_m_k3));
  D2_V_k4 = (p.D2_Iapp + ((-(( p.D2_naCurrentMSN_g_na*((D2_naCurrentMSN_m_last + dt*D2_naCurrentMSN_m_k3)).^3.*((D2_naCurrentMSN_h_last + dt*D2_naCurrentMSN_h_k3)).*(((D2_V_last + dt*D2_V_k3))-p.D2_naCurrentMSN_E_na))))+((-((p.D2_kCurrentMSN_g_k*((D2_kCurrentMSN_m_last + dt*D2_kCurrentMSN_m_k3)).^4.*(((D2_V_last + dt*D2_V_k3))-p.D2_kCurrentMSN_E_k))))+((-(( p.D2_mCurrentMSN_g_m*((D2_mCurrentMSN_m_last + dt*D2_mCurrentMSN_m_k3)).*(((D2_V_last + dt*D2_V_k3))-p.D2_mCurrentMSN_E_m))))+((-((p.D2_leakCurrentMSN_g_l*(((D2_V_last + dt*D2_V_k3))-p.D2_leakCurrentMSN_E_l))))+(((( p.D2_injectedCurrentD2_injectedCurrent + p.D2_injectedCurrentD2_sinmult*sin(2*pi*D2_injectedCurrentD2_freq*t/1000) - p.D2_injectedCurrentD2_DAmult*p.D2_injectedCurrentD2_DA)))+(((( p.D2_noisyInputMSN_sigma_noise.*randn(1,p.D2_Npop).*sqrt(dt))))+(((( D2_AMPAMSN_cellmask1.*D2_AMPAMSN_psptime(round(t)+1) + D2_AMPAMSN_cellmask2.*D2_AMPAMSN_psptime2(round(t)+1))))+((-(( (p.D2_soma_somaMSNiSYN_m_gsyn.*(((D2_soma_somaMSNiSYN_s_last + dt*D2_soma_somaMSNiSYN_s_k3))*D2_soma_somaMSNiSYN_mask).*(((D2_V_last + dt*D2_V_k3))-p.D2_soma_somaMSNiSYN_Esyn)))))+((-(( p.D2_D2_gabaRecInputMSN_g_gaba.*(((D2_D2_gabaRecInputMSN_s_last + dt*D2_D2_gabaRecInputMSN_s_k3))*D2_D2_gabaRecInputMSN_netcon).*(((D2_V_last + dt*D2_V_k3))-p.D2_D2_gabaRecInputMSN_E_gaba))))))))))))) )/p.D2_Cm;
  D2_naCurrentMSN_m_k4 = (( 0.32*(((D2_V_last + dt*D2_V_k3))+54)./(1-exp(-(((D2_V_last + dt*D2_V_k3))+54)/4)))).*(1-((D2_naCurrentMSN_m_last + dt*D2_naCurrentMSN_m_k3)))-(( 0.28*(((D2_V_last + dt*D2_V_k3))+27)./(exp((((D2_V_last + dt*D2_V_k3))+27)/5)-1))).*((D2_naCurrentMSN_m_last + dt*D2_naCurrentMSN_m_k3));
  D2_naCurrentMSN_h_k4 = (( 0.128*exp(-(((D2_V_last + dt*D2_V_k3))+50)/18))).*(1-((D2_naCurrentMSN_h_last + dt*D2_naCurrentMSN_h_k3)))-(( 4./(1+exp(-(((D2_V_last + dt*D2_V_k3))+27)/5)))).*((D2_naCurrentMSN_h_last + dt*D2_naCurrentMSN_h_k3));
  D2_kCurrentMSN_m_k4 = (( 0.032*(((D2_V_last + dt*D2_V_k3))+52)./(1-exp(-(((D2_V_last + dt*D2_V_k3))+52)/5)))).*(1-((D2_kCurrentMSN_m_last + dt*D2_kCurrentMSN_m_k3)))-(( 0.5*exp(-(((D2_V_last + dt*D2_V_k3))+57)/40))).*((D2_kCurrentMSN_m_last + dt*D2_kCurrentMSN_m_k3));
  D2_mCurrentMSN_m_k4 = (( D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + dt*D2_V_k3))-p.D2_mCurrentMSN_vhalf)./(1-exp(-(((D2_V_last + dt*D2_V_k3))-p.D2_mCurrentMSN_vhalf)/9)))).*(1-((D2_mCurrentMSN_m_last + dt*D2_mCurrentMSN_m_k3))) - (( -D2_mCurrentMSN_Qs*1e-4*(((D2_V_last + dt*D2_V_k3))-p.D2_mCurrentMSN_vhalf)./(1-exp((((D2_V_last + dt*D2_V_k3))-p.D2_mCurrentMSN_vhalf)/9)))).*((D2_mCurrentMSN_m_last + dt*D2_mCurrentMSN_m_k3));
  soma_soma_somaSomaiSYN_s_k4 = -((soma_soma_somaSomaiSYN_s_last + dt*soma_soma_somaSomaiSYN_s_k3))./p.soma_soma_somaSomaiSYN_tauD + ((1-((soma_soma_somaSomaiSYN_s_last + dt*soma_soma_somaSomaiSYN_s_k3)))/p.soma_soma_somaSomaiSYN_tauR).*(1+tanh(((soma_V_last + dt*soma_V_k3))/10));
  D1_soma_somaMSNiSYN_s_k4 = -((D1_soma_somaMSNiSYN_s_last + dt*D1_soma_somaMSNiSYN_s_k3))./p.D1_soma_somaMSNiSYN_tauD + ((1-((D1_soma_somaMSNiSYN_s_last + dt*D1_soma_somaMSNiSYN_s_k3)))/p.D1_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + dt*soma_V_k3))/10));
  D1_D1_gabaRecInputMSN_s_k4 = -((D1_D1_gabaRecInputMSN_s_last + dt*D1_D1_gabaRecInputMSN_s_k3))./p.D1_D1_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D1_V_last + dt*D1_V_k3))/4)).*(1-((D1_D1_gabaRecInputMSN_s_last + dt*D1_D1_gabaRecInputMSN_s_k3)));
  D2_soma_somaMSNiSYN_s_k4 = -((D2_soma_somaMSNiSYN_s_last + dt*D2_soma_somaMSNiSYN_s_k3))./p.D2_soma_somaMSNiSYN_tauD + ((1-((D2_soma_somaMSNiSYN_s_last + dt*D2_soma_somaMSNiSYN_s_k3)))/p.D2_soma_somaMSNiSYN_tauR).*(1+tanh(((soma_V_last + dt*soma_V_k3))/10));
  D2_D2_gabaRecInputMSN_s_k4 = -((D2_D2_gabaRecInputMSN_s_last + dt*D2_D2_gabaRecInputMSN_s_k3))./p.D2_D2_gabaRecInputMSN_tau_gaba + 2*(1+tanh(((D2_V_last + dt*D2_V_k3))/4)).*(1-((D2_D2_gabaRecInputMSN_s_last + dt*D2_D2_gabaRecInputMSN_s_k3)));

  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  soma_V_last = soma_V_last+(dt/6)*(soma_V_k1 + 2*(soma_V_k2 + soma_V_k3) + soma_V_k4);
  soma_somaGolombK_a_last = soma_somaGolombK_a_last+(dt/6)*(soma_somaGolombK_a_k1 + 2*(soma_somaGolombK_a_k2 + soma_somaGolombK_a_k3) + soma_somaGolombK_a_k4);
  soma_somaGolombK_b_last = soma_somaGolombK_b_last+(dt/6)*(soma_somaGolombK_b_k1 + 2*(soma_somaGolombK_b_k2 + soma_somaGolombK_b_k3) + soma_somaGolombK_b_k4);
  soma_somaGolombKdr_n_last = soma_somaGolombKdr_n_last+(dt/6)*(soma_somaGolombKdr_n_k1 + 2*(soma_somaGolombKdr_n_k2 + soma_somaGolombKdr_n_k3) + soma_somaGolombKdr_n_k4);
  soma_somaGolombNa_h_last = soma_somaGolombNa_h_last+(dt/6)*(soma_somaGolombNa_h_k1 + 2*(soma_somaGolombNa_h_k2 + soma_somaGolombNa_h_k3) + soma_somaGolombNa_h_k4);
  dend_V_last = dend_V_last+(dt/6)*(dend_V_k1 + 2*(dend_V_k2 + dend_V_k3) + dend_V_k4);
  dend_dendGolombK_a_last = dend_dendGolombK_a_last+(dt/6)*(dend_dendGolombK_a_k1 + 2*(dend_dendGolombK_a_k2 + dend_dendGolombK_a_k3) + dend_dendGolombK_a_k4);
  dend_dendGolombK_b_last = dend_dendGolombK_b_last+(dt/6)*(dend_dendGolombK_b_k1 + 2*(dend_dendGolombK_b_k2 + dend_dendGolombK_b_k3) + dend_dendGolombK_b_k4);
  dend_dendGolombKdr_n_last = dend_dendGolombKdr_n_last+(dt/6)*(dend_dendGolombKdr_n_k1 + 2*(dend_dendGolombKdr_n_k2 + dend_dendGolombKdr_n_k3) + dend_dendGolombKdr_n_k4);
  dend_dendGolombNa_h_last = dend_dendGolombNa_h_last+(dt/6)*(dend_dendGolombNa_h_k1 + 2*(dend_dendGolombNa_h_k2 + dend_dendGolombNa_h_k3) + dend_dendGolombNa_h_k4);
  D1_V_last = D1_V_last+(dt/6)*(D1_V_k1 + 2*(D1_V_k2 + D1_V_k3) + D1_V_k4);
  D1_naCurrentMSN_m_last = D1_naCurrentMSN_m_last+(dt/6)*(D1_naCurrentMSN_m_k1 + 2*(D1_naCurrentMSN_m_k2 + D1_naCurrentMSN_m_k3) + D1_naCurrentMSN_m_k4);
  D1_naCurrentMSN_h_last = D1_naCurrentMSN_h_last+(dt/6)*(D1_naCurrentMSN_h_k1 + 2*(D1_naCurrentMSN_h_k2 + D1_naCurrentMSN_h_k3) + D1_naCurrentMSN_h_k4);
  D1_kCurrentMSN_m_last = D1_kCurrentMSN_m_last+(dt/6)*(D1_kCurrentMSN_m_k1 + 2*(D1_kCurrentMSN_m_k2 + D1_kCurrentMSN_m_k3) + D1_kCurrentMSN_m_k4);
  D1_mCurrentMSN_m_last = D1_mCurrentMSN_m_last+(dt/6)*(D1_mCurrentMSN_m_k1 + 2*(D1_mCurrentMSN_m_k2 + D1_mCurrentMSN_m_k3) + D1_mCurrentMSN_m_k4);
  D2_V_last = D2_V_last+(dt/6)*(D2_V_k1 + 2*(D2_V_k2 + D2_V_k3) + D2_V_k4);
  D2_naCurrentMSN_m_last = D2_naCurrentMSN_m_last+(dt/6)*(D2_naCurrentMSN_m_k1 + 2*(D2_naCurrentMSN_m_k2 + D2_naCurrentMSN_m_k3) + D2_naCurrentMSN_m_k4);
  D2_naCurrentMSN_h_last = D2_naCurrentMSN_h_last+(dt/6)*(D2_naCurrentMSN_h_k1 + 2*(D2_naCurrentMSN_h_k2 + D2_naCurrentMSN_h_k3) + D2_naCurrentMSN_h_k4);
  D2_kCurrentMSN_m_last = D2_kCurrentMSN_m_last+(dt/6)*(D2_kCurrentMSN_m_k1 + 2*(D2_kCurrentMSN_m_k2 + D2_kCurrentMSN_m_k3) + D2_kCurrentMSN_m_k4);
  D2_mCurrentMSN_m_last = D2_mCurrentMSN_m_last+(dt/6)*(D2_mCurrentMSN_m_k1 + 2*(D2_mCurrentMSN_m_k2 + D2_mCurrentMSN_m_k3) + D2_mCurrentMSN_m_k4);
  soma_soma_somaSomaiSYN_s_last = soma_soma_somaSomaiSYN_s_last+(dt/6)*(soma_soma_somaSomaiSYN_s_k1 + 2*(soma_soma_somaSomaiSYN_s_k2 + soma_soma_somaSomaiSYN_s_k3) + soma_soma_somaSomaiSYN_s_k4);
  D1_soma_somaMSNiSYN_s_last = D1_soma_somaMSNiSYN_s_last+(dt/6)*(D1_soma_somaMSNiSYN_s_k1 + 2*(D1_soma_somaMSNiSYN_s_k2 + D1_soma_somaMSNiSYN_s_k3) + D1_soma_somaMSNiSYN_s_k4);
  D1_D1_gabaRecInputMSN_s_last = D1_D1_gabaRecInputMSN_s_last+(dt/6)*(D1_D1_gabaRecInputMSN_s_k1 + 2*(D1_D1_gabaRecInputMSN_s_k2 + D1_D1_gabaRecInputMSN_s_k3) + D1_D1_gabaRecInputMSN_s_k4);
  D2_soma_somaMSNiSYN_s_last = D2_soma_somaMSNiSYN_s_last+(dt/6)*(D2_soma_somaMSNiSYN_s_k1 + 2*(D2_soma_somaMSNiSYN_s_k2 + D2_soma_somaMSNiSYN_s_k3) + D2_soma_somaMSNiSYN_s_k4);
  D2_D2_gabaRecInputMSN_s_last = D2_D2_gabaRecInputMSN_s_last+(dt/6)*(D2_D2_gabaRecInputMSN_s_k1 + 2*(D2_D2_gabaRecInputMSN_s_k2 + D2_D2_gabaRecInputMSN_s_k3) + D2_D2_gabaRecInputMSN_s_k4);

  if mod(k,downsample_factor)==0 % store this time point

  % ------------------------------------------------------------
  % Store state variables:
  % ------------------------------------------------------------
    soma_V(n,:) = soma_V_last;
    soma_somaGolombK_a(n,:) = soma_somaGolombK_a_last;
    soma_somaGolombK_b(n,:) = soma_somaGolombK_b_last;
    soma_somaGolombKdr_n(n,:) = soma_somaGolombKdr_n_last;
    soma_somaGolombNa_h(n,:) = soma_somaGolombNa_h_last;
    dend_V(n,:) = dend_V_last;
    dend_dendGolombK_a(n,:) = dend_dendGolombK_a_last;
    dend_dendGolombK_b(n,:) = dend_dendGolombK_b_last;
    dend_dendGolombKdr_n(n,:) = dend_dendGolombKdr_n_last;
    dend_dendGolombNa_h(n,:) = dend_dendGolombNa_h_last;
    D1_V(n,:) = D1_V_last;
    D1_naCurrentMSN_m(n,:) = D1_naCurrentMSN_m_last;
    D1_naCurrentMSN_h(n,:) = D1_naCurrentMSN_h_last;
    D1_kCurrentMSN_m(n,:) = D1_kCurrentMSN_m_last;
    D1_mCurrentMSN_m(n,:) = D1_mCurrentMSN_m_last;
    D2_V(n,:) = D2_V_last;
    D2_naCurrentMSN_m(n,:) = D2_naCurrentMSN_m_last;
    D2_naCurrentMSN_h(n,:) = D2_naCurrentMSN_h_last;
    D2_kCurrentMSN_m(n,:) = D2_kCurrentMSN_m_last;
    D2_mCurrentMSN_m(n,:) = D2_mCurrentMSN_m_last;
    soma_soma_somaSomaiSYN_s(n,:) = soma_soma_somaSomaiSYN_s_last;
    D1_soma_somaMSNiSYN_s(n,:) = D1_soma_somaMSNiSYN_s_last;
    D1_D1_gabaRecInputMSN_s(n,:) = D1_D1_gabaRecInputMSN_s_last;
    D2_soma_somaMSNiSYN_s(n,:) = D2_soma_somaMSNiSYN_s_last;
    D2_D2_gabaRecInputMSN_s(n,:) = D2_D2_gabaRecInputMSN_s_last;

    n=n+1;
  end
end

T=T(1:downsample_factor:ntime);

end

