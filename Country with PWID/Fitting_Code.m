%Authors: Hiam Chemaitelly and Houssein H. Ayoub
%Date: 14 February 2022
%Code for FSWs in Iran

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

IDUF=0.136; %Assuming that 4% of FSWs are IDUs (from Systematic review)

muHIV   = 1/(2*Unit); % death probability by HIV in last stage per day
mu      = 1.0/(35*Unit);  % death by non-HIV 

nact_r=3; % number of coital acts per month by regular clients
nact_nr=1; % number of coital acts per month by non-regular clients
nact_sp=25; %number of coital acts per year with spouse
Cond_sp=0.029; %2.9% of acts with spouses are coovered by condom use

eART=0.57;
eMC=0.58;
eCond=0.8;
ePrEP=0.51; 

tpHIVregpercoitalact= [0.0360,0.0008,0.0042]; % *Unit *DT transmission probability of HIV per coital act in AIDS for regular clients (assuming 3 time per week)
tpHIVnpercoitalact= [0.0360,0.0008,0.0042]; % *Unit *DT transmission probability of HIV per coital act in AIDS for non-regular clients (!!SHOULD CONFIRM NON-REGULAR CLIENTS do one time coital act!!)

tpHIVspouse=((0.036*49/365)+(0.0008*9)+(0.0042*2))/((49/365)+9+2);
tpsp=1-((1-tpHIVspouse).^(nact_sp*(1-Cond_sp))*((1-(1-eCond)*tpHIVspouse).^(nact_sp*Cond_sp)));

tpHIVreg=(ones(1,3)-(ones(1,3)-tpHIVregpercoitalact).^(nact_r*1));%*DT;
tpHIVn=(ones(1,3)-(ones(1,3)-tpHIVnpercoitalact).^(nact_nr*1));%*DT;

transHIV = [365.0/49.0*DT, 1.0/9.0*DT]; % transition probability of HIV from 3 to 2 per day

removelink_reg  = 1.0/((3.0*30.0/365.0))*DT;             % removal rate of partnership for regular clients (Partnership duration regular)
removelink_nreg = 1.0/((1.0*30.0/365.0))*DT;             % removal rate of partnership for non-regular clients (Partnership duration non-regular)

YearEND=380;
YearBURNIN=50;
BURN_IN = YearBURNIN*Unit; % the length of "burn-in" period  (year)  % Let us say that the burn in is on 1940, then we must to wait until 1990 to reach equilibirm point.
TEND    = YearEND*Unit; %the length of simulation (year)          %% We similulate 50 years after 1990.

% seeds of HIV epidemics at initial time
INIT_INF_fsw = 0.03*fsw; % seeds of HIV epidemics at initial time among fsw
INIT_INF_regm = 0.01*regm; % seeds of HIV epidemics at initial time among regm
INIT_INF_nregm = 0.01*nregm; % seeds of HIV epidemics at initial time among nregm


%%%Time trend for intervention inclsuion
CoverageMC=zeros(TEND,1);
CovARTFSW=zeros(TEND,1);
CovARTClient=zeros(TEND,1);
CoverageCond=zeros(TEND,1);
CovPrEPFSW=zeros(TEND,1);
CovPrEPClient=zeros(TEND,1);

t0=2030-YearEND+1;
tf=2030;
for t=t0+YearBURNIN:DT:tf+1-DT
    %%%Coverage MC
    if t<2020
      CoverageMC(round((t-t0)*Unit+1))=0.997; % This is the baseline coverage for MC for South-Sudan
    else
      CoverageMC(round((t-t0)*Unit+1))=0.997;  %%This needs to be updated according to the intervention scenario
    end
    
    %%%Coverage Cond
    if t<2020
      CoverageCond(round((t-t0)*Unit+1))=0.571; % This is the baseline coverage for Condome use for South-Sudan
    else
      CoverageCond(round((t-t0)*Unit+1))=0.571;  %%This needs to be updated according to the intervention scenario (fox example 0.5 or 0.81)
    end
    
    %%%Coverage ART
    if t<2020
      CovARTFSW(round((t-t0)*Unit+1))=0.20; % This is the baseline coverage for ART in FSWs for South-Sudan
    else
      CovARTFSW(round((t-t0)*Unit+1))=0.20;  %%This needs to be updated according to the intervention scenario (fox example 0.25 or 0.5 or 0.81)
    end
    
    if t<2020
      CovARTClient(round((t-t0)*Unit+1))=0.20; % This is the baseline coverage for ART for South-Sudan for general population
    else
      CovARTClient(round((t-t0)*Unit+1))=0.20;  %%This needs to be updated according to the intervention scenario (fox example 0.25 or 0.5 or 0.81)
    end
    
    %%%Coverage PrEP
    if t<2020
      CovPrEPFSW(round((t-t0)*Unit+1))=0; % This is the baseline coverage for PrEP for South-Sudan-No PrEP in South Sudan
    else
      CovPrEPFSW(round((t-t0)*Unit+1))=0;  %%This needs to be updated according to the intervention scenario (fox example 0.25 or 0.5 or 0.81)
    end
    
    if t<2020
      CovPrEPClient(round((t-t0)*Unit+1))=0; % This is the baseline coverage for PrEP for South-Sudan-No PrEP in South Sudan
    else
      CovPrEPClient(round((t-t0)*Unit+1))=0;  %%This needs to be updated according to the intervention scenario (fox example 0.25 or 0.5 or 0.81)
    end
end


xx=zeros(2,1);

xx(1)=0.3; %a
lb(1)=0.08;
ub(1)=0.6;

xx(2)=1.0; %c
lb(2)=0.05;
ub(2)=3.0;

DataFSWHIV=3.3;
DataFSWHIVPWID=9.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=optimset('Display','iter','TolFun',1e-3,'MaxIter',50,'MaxFunEvals',50);
tic
NumberofBestFits=50;
NumberofParameters=2;
EstimatedParameters=zeros(NumberofParameters,NumberofBestFits);
Cost=zeros(NumberofBestFits,1);
parfor n=1:NumberofBestFits
    [EstimatedParameters(:,n),Cost(n)]=fminsearchbnd(@CostFunction,xx,lb,ub,options,IDUF,DataFSWHIV,DataFSWHIVPWID,CovPrEPClient,CovPrEPFSW,CovARTClient,CovARTFSW,CoverageCond,CoverageMC,t0,tf,INIT_INF_nregm,INIT_INF_regm,INIT_INF_fsw,TEND,BURN_IN,YearBURNIN,YearEND,removelink_nreg,removelink_reg,transHIV,tpHIVn,tpHIVreg,tpHIVnpercoitalact,tpHIVregpercoitalact,ePrEP,eCond,eMC,eART,nact_nr,nact_r,mu,muHIV,nregm,regm,fsw,Unit,DT);
    if Cost(n)>0.03
       [EstimatedParameters(:,n),Cost(n)]=fminsearchbnd(@CostFunction,xx,lb,ub,options,IDUF,DataFSWHIV,DataFSWHIVPWID,CovPrEPClient,CovPrEPFSW,CovARTClient,CovARTFSW,CoverageCond,CoverageMC,t0,tf,INIT_INF_nregm,INIT_INF_regm,INIT_INF_fsw,TEND,BURN_IN,YearBURNIN,YearEND,removelink_nreg,removelink_reg,transHIV,tpHIVn,tpHIVreg,tpHIVnpercoitalact,tpHIVregpercoitalact,ePrEP,eCond,eMC,eART,nact_nr,nact_r,mu,muHIV,nregm,regm,fsw,Unit,DT);
    end   
end
toc  
