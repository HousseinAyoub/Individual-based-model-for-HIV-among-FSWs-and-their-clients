
function Cost =CostFunction(xx,IDUF,DataFSWHIV,DataFSWHIVPWID,CovPrEPClient,CovPrEPFSW,CovARTClient,CovARTFSW,CoverageCond,CoverageMC,t0,tf,INIT_INF_nregm,INIT_INF_regm,INIT_INF_fsw,TEND,BURN_IN,YearBURNIN,YearEND,removelink_nreg,removelink_reg,transHIV,tpHIVn,tpHIVreg,tpHIVnpercoitalact,tpHIVregpercoitalact,ePrEP,eCond,eMC,eART,nact_nr,nact_r,mu,muHIV,nregm,regm,fsw,Unit,DT)

a=xx(1); 
HIVPWID=(0.099/5)*DT*a; %Likelihood of a FSW acquiring HIV through IDU (Likelihood= 9.9% prevalence/5 years [assuming 10 years duration of infection and 10 years before exiting IDU, that is 1/10+1/10=1/5years] *1/12 for the time step per month)

c=xx(2); 

mean_reg        = c*2.0/10 ;  % mean acquistion probabiltiy of new partner for regular client, this corresponds to each sex worker having 2 regular clients per month
var_reg         = 0.25*c*2.0/10;%0.01;  % variance acquistion probabiltiy of new partner for regular client 
mean_nreg       = c*3.0/10;  % mean acquistion probabiltiy of new partner for non-regular client, this corresponds to each sex worker having 6 partnerships with non-regular clients per year.
var_nreg        = 0.25*c*3.0/10; %0.01;  % variance acquistion probabiltiy of new partner for non-regular client
PAR_reg         = [mean_reg^2/var_reg, var_reg/mean_reg];     % parameters for acquistion probabiltiy of new partner for regular client
PAR_nreg        = [mean_nreg^2/var_nreg, var_nreg/mean_nreg];    % parameters for acquistion probabiltiy of new partner for regular client
acq_regm    = zeros(1,regm); % acquistion probability of new partner for regular clients
acq_nregm   = zeros(1,nregm); % acquistion probability of new partner for non-regular clients
%% initialization
for cnt = 1:regm
    acq_regm(cnt) = gamrnd(PAR_reg(1), PAR_reg(2));
end

for cnt = 1:nregm
    acq_nregm(cnt) = gamrnd(PAR_nreg(1), PAR_nreg(2));
end

[PrevalenceYear1_Avg,PrevalenceYear1IDU_Avg]=Model_Interv_IDUFitting(IDUF,HIVPWID,DT,Unit,fsw,regm,nregm,tpHIVreg,tpHIVn,transHIV,removelink_reg,removelink_nreg,YearEND,YearBURNIN,BURN_IN,TEND,INIT_INF_fsw,INIT_INF_regm,INIT_INF_nregm,CoverageMC,CoverageCond,CovARTFSW,CovARTClient,CovPrEPFSW,CovPrEPClient,eART,eMC,eCond,ePrEP,acq_regm,acq_nregm,mu,muHIV);

Cost=0.8*((PrevalenceYear1_Avg(2020-t0+1)-DataFSWHIV)/DataFSWHIV)^2+0.2*((PrevalenceYear1IDU_Avg(2020-t0+1)-DataFSWHIVPWID)/DataFSWHIVPWID)^2;
