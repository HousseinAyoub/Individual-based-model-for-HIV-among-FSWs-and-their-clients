function [MeannumberofContactreg,MeannumberofContactnreg,PrevalenceYear1,PrevalenceYear2,PrevalenceYear3,PrevalenceYear23,PrevalenceYear4,PrevalenceYear1IDU,PrevalenceYear1notIDU,PrevalenceYear1_Avg,PrevalenceYear2_Avg,PrevalenceYear3_Avg,PrevalenceYear23_Avg,PrevalenceYear4_Avg,PrevalenceYear1IDU_Avg,PrevalenceYear1notIDU_Avg,IncidenceR1Year,IncidenceR2Year,IncidenceR3Year,IncidenceR23Year,IncidenceR4Year,IncidenceR1YearIDU,IncidenceR1YearnotIDU,IncidenceR1Year_Avg,IncidenceR2Year_Avg,IncidenceR3Year_Avg,IncidenceR23Year_Avg,IncidenceR4Year_Avg,IncidenceR1YearIDU_Avg,IncidenceR1YearnotIDU_Avg,IncContYearS,IncContYearIDU,Incidence1YearS,Incidence1YearIDU,Incidence1YearnotIDU,Incidence1Year,Incidence2YearS,Incidence3YearS,Incidence23YearS,Incidence4YearS,IncidenceALL,Incidence1YearS_Avg,Incidence1YearIDU_Avg,Incidence1YearnotIDU_Avg,Incidence1Year_Avg,Incidence2YearS_Avg,Incidence3YearS_Avg,Incidence23YearS_Avg,Incidence4YearS_Avg,IncidenceALL_Avg,Number_I_FSW,Number_I_FSWIDU,Number_I_FSWnotIDU,Ratio,Inc_sp,IncR_sp,Inc_sp_Avg,IncR_sp_Avg,IncR_spr,IncR_spnr,Inc_spr, Inc_spnr,IncR_spr_Avg,IncR_spnr_Avg,Inc_spr_Avg,Inc_spnr_Avg,IncidenceContSexual,IncidenceContInjecting,NPWID,NnotPWID,Susc,Susc_IDU,Susc_notIDU,SuscRClient,SuscNRClient,SuscClient,Susc_Avg,Susc_IDU_Avg,Susc_notIDU_Avg,SuscRClient_Avg,SuscNRClient_Avg,SuscClient_Avg,InfectedFSW,InfectedFSWIDU,InfectedFSWnotIDU,InfectedRClients,InfectedNRClients,InfectedClients,InfectedALL,InfectedFSW_Avg,InfectedFSWIDU_Avg,InfectedFSWnotIDU_Avg,InfectedRClients_Avg,InfectedNRClients_Avg,InfectedClients_Avg,InfectedALL_Avg,Ndisc,Ndiscr,Ndiscnr,Ndisc_Avg,Ndiscr_Avg,Ndiscnr_Avg,Susc_spr,Susc_spnr,Susc_sp,Susc_spr_Avg,Susc_spnr_Avg,Susc_sp_Avg]= Multiple_Runs_for_Average(XXX)

tic
XXX;

load('job15_output.mat','job15_output')
EstimatedParameters=job15_output{1,1};
MeanEstimatedParameters=zeros(2,1);
    
%figure
[BinHeight,BinCenter]=createFit(EstimatedParameters(1,:));
[M,I]=max(BinHeight); %I is the bin number, M is the density
MeanEstimatedParameters(1)=0.75.*BinCenter(I);

%figure
[BinHeight,BinCenter]=createFit(EstimatedParameters(2,:));
[M,I]=max(BinHeight);
MeanEstimatedParameters(2)=1.3.*BinCenter(I);  



%%FIXED POPULATION VERSION
%%parameter sets
DT      = 1/12; % definition of unit time (delta t)
Unit=12; %% months

fsw     =600; % the number of FSW (FSW:clients = 1:10)
regm    =2000; % the number of regular clients (assuming 30% of clients)
nregm   =4000; % the number of non-regular clients
Marital=0.564; %fraction of clients who are married

npartr=regm*Marital; %number of partnerships for regular clients
npartnr=nregm*Marital; %number of partnerships for non-regular clients
npart=(regm+nregm)*Marital; %total number of partnerships for regular and non-reglar clients

IDUF=0.136; %Assuming that 1% of FSWs are IDUs (from Systematic review)

muHIV   = 1/(2*Unit); % death probability by HIV in last stage per day
mu      = 1.0/(35*Unit);  % death by non-HIV 

nact_r=3; % number of coital acts per month by regular clients
nact_nr=1; % number of coital acts per month by non-regular clients
nact_sp=25; %number of coital acts per year with spouse
Cond_sp=0.029; %1.5% of acts with spouses are coovered by condom use

eART=0.57; %Pooled efficacy of ART
eMC=0.58;  %Pooled efficacy of MC
eCond=0.8; %Pooled efficacy of cond
ePrEP=0.51; %Pooled efficacy of PrEP 

tpHIVregpercoitalact= [0.0360,0.0008,0.0042]; % *Unit *DT transmission probability of HIV per coital act in AIDS for regular clients (assuming 3 time per week)
tpHIVnpercoitalact= [0.0360,0.0008,0.0042]; % *Unit *DT transmission probability of HIV per coital act in AIDS for non-regular clients (!!SHOULD CONFIRM NON-REGULAR CLIENTS do one time coital act!!)


tpHIVspouse=((0.036*49/365)+(0.0008*9)+(0.0042*2))/((49/365)+9+2);
tpsp=1-((1-tpHIVspouse).^(nact_sp*(1-Cond_sp))*((1-(1-eCond)*tpHIVspouse).^(nact_sp*Cond_sp)));

tpHIVreg=(ones(1,3)-(ones(1,3)-tpHIVregpercoitalact).^(nact_r*1));%*DT;
tpHIVn=(ones(1,3)-(ones(1,3)-tpHIVnpercoitalact).^(nact_nr*1));%*DT;

transHIV = [365.0/49.0*DT, 1.0/9.0*DT]; % transition probability of HIV from 3 to 2 per day

removelink_reg  = 1.0/((3.0*30.0/365.0))*DT;             % removal rate of partnership for regular clients (Partnership duration regular)
removelink_nreg = 1.0/((1.0*30.0/365.0))*DT;             % removal rate of partnership for non-regular clients (Partnership duration non-regular)

% Let us assume that the simulation start at 1890 until 1940 to stabilize networks. The burn in is assumed at 1940, then we must to wait until 1990 to reach equilibirm point.
YearEND=400;
YearBURNIN=50;
BURN_IN = YearBURNIN*Unit; % the length of "burn-in" period  (year)  % Let us say that the burn in is on 1940, then we must to wait until 1990 to reach equilibirm point.
TEND    = YearEND*Unit; %the length of simulation (year)          %% We similulate 50 years after 1990.

% seeds of HIV epidemics at initial time
INIT_INF_fsw = 0.03*fsw; % seeds of HIV epidemics at initial time among fsw
INIT_INF_regm = 0.01*regm; % seeds of HIV epidemics at initial time among regm
INIT_INF_nregm = 0.01*nregm; % seeds of HIV epidemics at initial time among nregm

t0=2050-YearEND+1;
tf=2050;

%%These need to be updated as per each country
DataCoverageMC=0.997; 
DataCoverageCond=0.571;
DataCovARTFSW=0.20;
DataCovARTClient=0.20;
DataCovPrEPFSW=0;
DataCovPrEPClient=0;

NumberofScenarios=12; 
%NumberofScenarios=1; 

%%%Time trend for intervention inclsuion
CoverageMC=zeros(TEND,1);
CovARTFSW=zeros(TEND,1);
CovARTClient=zeros(TEND,1);
CoverageCond=zeros(TEND,1);
CovPrEPFSW=zeros(TEND,1);
CovPrEPClient=zeros(TEND,1);

eARTT=zeros(TEND,1);

for n=1:NumberofScenarios 
     for t=t0+YearBURNIN:DT:tf+1-DT
              CovARTClient(round((t-t0)*Unit+1))=DataCovARTClient; %%No intervention for client
              CovPrEPClient(round((t-t0)*Unit+1))=DataCovPrEPClient; %%No intervention for client
              CoverageMC(round((t-t0)*Unit+1))=DataCoverageMC;  %%Not applicable as coverage MC in Morocco is already 99.9%
 
               if n==1 %%Baseline
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW; 
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
              
               
               %%%%%%%%%%%%%%%Condom Scenarios%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               elseif n==2 %%Increase Cond to 80%
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=0.8;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               end  
                              
               %%%%%%%%%%%%%%%ART Scenarios%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               elseif n==3 %%ART=25% eART=0.57
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=0.25;  
               end 
               
               elseif n==4 %%ART=50% eART=0.57
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=0.50;  
               end
               
               elseif n==5 %%ART=81% eART=0.57
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=0.81;  
               end
               
               elseif n==6 %%ART=25% eART=0.96
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=0.96;
               CovARTFSW(round((t-t0)*Unit+1))=0.25;  
               end 
               
               elseif n==7 %%ART=50% eART=0.96
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=0.96;
               CovARTFSW(round((t-t0)*Unit+1))=0.50;  
               end 
               
               elseif n==8 %%ART=81% eART=0.96
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=0.96;
               CovARTFSW(round((t-t0)*Unit+1))=0.81;  
               end 
               
               %%%%%%%%%%%%%%%Prep Scenarios%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               elseif n==9 %%PREP=25%
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=0.25;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               end 
               
               elseif n==10 %%PREP=50%
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=0.50;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               end 
               
               %%%%%Moderately optimistic: Cond at baseline,
               %%%%%PREP=25%, ART=50% (effecacy=96%)               
               elseif n==11
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW;
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=0.25;
               eARTT(round((t-t0)*Unit+1))=0.96;
               CovARTFSW(round((t-t0)*Unit+1))=0.50;
               end 
              
              %%%%%%%%%%%%%%%Most Scenario%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              else %%%%%Most optimistic: Cond increase to 80%,
               %%%%%PREP=50%, ART=81% (effecacy=96%)
               if t<2020
               CoverageCond(round((t-t0)*Unit+1))=DataCoverageCond;   
               CovPrEPFSW(round((t-t0)*Unit+1))=DataCovPrEPFSW; 
               eARTT(round((t-t0)*Unit+1))=eART;
               CovARTFSW(round((t-t0)*Unit+1))=DataCovARTFSW;  
               else
               CoverageCond(round((t-t0)*Unit+1))=0.8;   
               CovPrEPFSW(round((t-t0)*Unit+1))=0.5;
               eARTT(round((t-t0)*Unit+1))=0.96;
               CovARTFSW(round((t-t0)*Unit+1))=0.81;  
               end    
             end
     end
      
a=MeanEstimatedParameters(1);   

HIVPWID=(0.099/5)*DT*a; %Likelihood of a FSW acquiring HIV through IDU (Likelihood= 9.9% prevalence/5 years [assuming 10 years duration of infection and 10 years before exiting IDU, that is 1/10+1/10=5years] *1/12 for the time step per month)

c=MeanEstimatedParameters(2);   

mean_reg        = c*2.0/10 ;  % mean acquistion probabiltiy of new partner for regular client, this corresponds to each sex worker having 2 regular clients per month
var_reg         = 0.25*c*2.0/10;  % variance acquistion probabiltiy of new partner for regular client 
mean_nreg       = c*3.0/10;  % mean acquistion probabiltiy of new partner for non-regular client, this corresponds to each sex worker having 3 partnerships with non-regular clients per year.
var_nreg        = 0.25*c*3.0/10;% variance acquistion probabiltiy of new partner for non-regular client
PAR_reg         = [mean_reg.^2/var_reg, var_reg/mean_reg];     % parameters for acquistion probabiltiy of new partner for regular client
PAR_nreg        = [mean_nreg.^2/var_nreg, var_nreg/mean_nreg];    % parameters for acquistion probabiltiy of new partner for regular client

PAR_reg1=PAR_reg(1);
PAR_reg2=PAR_reg(2);
PAR_nreg1=PAR_nreg(1);
PAR_nreg2=PAR_nreg(2);

parfor m=1:500
    [MeannumberofContactreg(:,m,n),MeannumberofContactnreg(:,m,n),PrevalenceYear1(:,m,n),PrevalenceYear2(:,m,n),PrevalenceYear3(:,m,n),PrevalenceYear23(:,m,n),PrevalenceYear4(:,m,n),PrevalenceYear1IDU(:,m,n),PrevalenceYear1notIDU(:,m,n),PrevalenceYear1_Avg(:,m,n),PrevalenceYear2_Avg(:,m,n),PrevalenceYear3_Avg(:,m,n),PrevalenceYear23_Avg(:,m,n),PrevalenceYear4_Avg(:,m,n),PrevalenceYear1IDU_Avg(:,m,n),PrevalenceYear1notIDU_Avg(:,m,n),IncidenceR1Year(:,m,n),IncidenceR2Year(:,m,n),IncidenceR3Year(:,m,n),IncidenceR23Year(:,m,n),IncidenceR4Year(:,m,n),IncidenceR1YearIDU(:,m,n),IncidenceR1YearnotIDU(:,m,n),IncidenceR1Year_Avg(:,m,n),IncidenceR2Year_Avg(:,m,n),IncidenceR3Year_Avg(:,m,n),IncidenceR23Year_Avg(:,m,n),IncidenceR4Year_Avg(:,m,n),IncidenceR1YearIDU_Avg(:,m,n),IncidenceR1YearnotIDU_Avg(:,m,n),IncContYearS(:,m,n),IncContYearIDU(:,m,n),Incidence1YearS(:,m,n),Incidence1YearIDU(:,m,n),Incidence1YearnotIDU(:,m,n),Incidence1Year(:,m,n),Incidence2YearS(:,m,n),Incidence3YearS(:,m,n),Incidence23YearS(:,m,n),Incidence4YearS(:,m,n),IncidenceALL(:,m,n),Incidence1YearS_Avg(:,m,n),Incidence1YearIDU_Avg(:,m,n),Incidence1YearnotIDU_Avg(:,m,n),Incidence1Year_Avg(:,m,n),Incidence2YearS_Avg(:,m,n),Incidence3YearS_Avg(:,m,n),Incidence23YearS_Avg(:,m,n),Incidence4YearS_Avg(:,m,n),IncidenceALL_Avg(:,m,n),Number_I_FSW(:,m,n),Number_I_FSWIDU(:,m,n),Number_I_FSWnotIDU(:,m,n),Ratio(:,m,n),Inc_sp(:,m,n),IncR_sp(:,m,n),Inc_sp_Avg(:,m,n),IncR_sp_Avg(:,m,n),IncR_spr(:,m,n),IncR_spnr(:,m,n),Inc_spr(:,m,n), Inc_spnr(:,m,n),IncR_spr_Avg(:,m,n),IncR_spnr_Avg(:,m,n),Inc_spr_Avg(:,m,n),Inc_spnr_Avg(:,m,n),IncidenceContSexual(:,m,n),IncidenceContInjecting(:,m,n),NPWID(:,m,n),NnotPWID(:,m,n),Susc(:,m,n),Susc_IDU(:,m,n),Susc_notIDU(:,m,n),SuscRClient(:,m,n),SuscNRClient(:,m,n),SuscClient(:,m,n),Susc_Avg(:,m,n),Susc_IDU_Avg(:,m,n),Susc_notIDU_Avg(:,m,n),SuscRClient_Avg(:,m,n),SuscNRClient_Avg(:,m,n),SuscClient_Avg(:,m,n),InfectedFSW(:,m,n),InfectedFSWIDU(:,m,n),InfectedFSWnotIDU(:,m,n),InfectedRClients(:,m,n),InfectedNRClients(:,m,n),InfectedClients(:,m,n),InfectedALL(:,m,n),InfectedFSW_Avg(:,m,n),InfectedFSWIDU_Avg(:,m,n),InfectedFSWnotIDU_Avg(:,m,n),InfectedRClients_Avg(:,m,n),InfectedNRClients_Avg(:,m,n),InfectedClients_Avg(:,m,n),InfectedALL_Avg(:,m,n),Ndisc(:,m,n),Ndiscr(:,m,n),Ndiscnr(:,m,n),Ndisc_Avg(:,m,n),Ndiscr_Avg(:,m,n),Ndiscnr_Avg(:,m,n),Susc_spr(:,m,n),Susc_spnr(:,m,n),Susc_sp(:,m,n),Susc_spr_Avg(:,m,n),Susc_spnr_Avg(:,m,n),Susc_sp_Avg(:,m,n)]=Model_Interv_IDU(IDUF,HIVPWID,DT,Unit,fsw,regm,nregm,tpHIVreg,tpHIVn,transHIV,removelink_reg,removelink_nreg,YearEND,YearBURNIN,BURN_IN,TEND,INIT_INF_fsw,INIT_INF_regm,INIT_INF_nregm,CoverageMC,CoverageCond,CovARTFSW,CovARTClient,CovPrEPFSW,CovPrEPClient,eART,eARTT,eMC,eCond,ePrEP,PAR_reg1,PAR_reg2,PAR_nreg1,PAR_nreg2,mu,muHIV,npartr,npartnr,npart,nact_sp,tpHIVspouse,Cond_sp,tpsp);
end
end
