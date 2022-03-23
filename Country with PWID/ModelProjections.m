    clear all
    clc
    eART=0.57;
    DataCovARTClient=0.20;
    FactorART=((1-DataCovARTClient)+(1-eART)*DataCovARTClient);

    load('job33_output.mat','job33_output') 
    IranScenario500=job33_output;
    
    PrevalenceYear1=IranScenario500{1,3};
    PrevalenceYear23=IranScenario500{1,6};
    PrevalenceYear1IDU=IranScenario500{1,8};
    
    IncidenceR1Year=IranScenario500{1,17};
    IncidenceR23Year=IranScenario500{1,20};
    IncidenceR1YearIDU=IranScenario500{1,22};

    Incidence1YearS=IranScenario500{1,33};
    Incidence1YearIDU=IranScenario500{1,34};
    Incidence1Year=IranScenario500{1,36};
    Incidence23YearS=IranScenario500{1,39};
    IncidenceALL=IranScenario500{1,41};

    
    Inc_sp=IranScenario500{1,55};
    IncR_sp=IranScenario500{1,56};
    
    Inc_spr=IranScenario500{1,61};
    IncR_spr=IranScenario500{1,59};
    
    Inc_spnr=IranScenario500{1,62};
    IncR_spnr=IranScenario500{1,60};

    Susc=IranScenario500{1,71};
    Susc_IDU=IranScenario500{1,72};
    SuscClient=IranScenario500{1,76};
    Susc_sp=IranScenario500{1,105};
   

NumberofScenarios=1; %we simulate here only the baseline, but his can be changed to simulate the other scenarios
  
DT      = 1/12; 
Unit=12; %% months
YearEND=400;
YearBURNIN=50;
BURN_IN = YearBURNIN*Unit; % the length of "burn-in" period  (year)  % Let us say that the burn in is on 1940, then we must to wait until 1990 to reach equilibirm point.
TEND    = YearEND*Unit; %the length of simulation (year)          %% We similulate 50 years after 1990.
t0=2050-YearEND+1;
tf=2050; 

dtt=1651:2050; 
ConsecutativeYear=3;
for m=1:500
  for t=YearEND-50:YearEND-ConsecutativeYear 
    if (sum(Incidence1Year(t:t+ConsecutativeYear,m,1))==0)||(sum(Incidence23YearS(t:t+ConsecutativeYear,m,1))==0)
       PrevalenceYear1N(:,m)=zeros(YearEND,1);
       PrevalenceYear23N(:,m)=zeros(YearEND,1);
       PrevalenceYear1IDUN(:,m)=zeros(YearEND,1);
       
       Incidence1YearN(:,m)=zeros(YearEND,1);
       Incidence1YearIDUN(:,m)=zeros(YearEND,1);
       Incidence23YearSN(:,m)=zeros(YearEND,1);
       Inc_spN(:,m)=zeros(YearEND,1);
       
       IncidenceR1YearN(:,m)=zeros(YearEND,1);
       IncidenceR1YearIDUN(:,m)=zeros(YearEND,1);
       IncidenceR23YearN(:,m)=zeros(YearEND,1);
       IncR_spN(:,m)=zeros(YearEND,1);
    else
       PrevalenceYear1N(:,m)=PrevalenceYear1(:,m,1);
       PrevalenceYear23N(:,m)=PrevalenceYear23(:,m,1);
       PrevalenceYear1IDUN(:,m)=PrevalenceYear1IDU(:,m,1);
       
       Incidence1YearN(:,m)=Incidence1Year(:,m,1);
       Incidence1YearIDUN(:,m)=Incidence1YearIDU(:,m,1);
       Incidence23YearSN(:,m)=Incidence23YearS(:,m,1);
       Inc_spN(:,m)=Inc_sp(:,m,1);
       
       IncidenceR1YearN(:,m)=IncidenceR1Year(:,m,1);
       IncidenceR1YearIDUN(:,m)=IncidenceR1YearIDU(:,m,1);
       IncidenceR23YearN(:,m)=IncidenceR23Year(:,m,1);
       IncR_spN(:,m)=IncR_sp(:,m,1);
    end
  end
    PrevalenceYear1N(:,~any(PrevalenceYear1N,1))=[];  %columns
    PrevalenceYear1N=sort(PrevalenceYear1N(:,:),2);
   PrevalenceYear1IDUN(:,~any(PrevalenceYear1IDUN,1))=[];  %columns
    PrevalenceYear1IDUN=sort(PrevalenceYear1IDUN(:,:),2);
    PrevalenceYear23N(:,~any(PrevalenceYear23N,1))=[];  %columns
    PrevalenceYear23N=sort(PrevalenceYear23N(:,:),2);
    Incidence1YearN(:,~any(Incidence1YearN,1))=[];  %columns
    Incidence1YearN=sort(Incidence1YearN(:,:),2);
    Incidence1YearIDUN(:,~any(Incidence1YearIDUN,1))=[];  %columns
    Incidence1YearIDUN=sort(Incidence1YearIDUN(:,:),2);
    Incidence23YearSN(:,~any(Incidence23YearSN,1))=[];  %columns
    Incidence23YearSN=sort(Incidence23YearSN(:,:),2);
    Inc_spN(:,~any(Inc_spN,1))=[];  %columns            
    Inc_spN=sort(Inc_spN(:,:),2);
    IncidenceR1YearN(:,~any(IncidenceR1YearN,1))=[];  %columns
    IncidenceR1YearN=sort(IncidenceR1YearN(:,:),2);
    IncidenceR1YearIDUN(:,~any(IncidenceR1YearIDUN,1))=[];  %columns
    IncidenceR1YearIDUN=sort(IncidenceR1YearIDUN(:,:),2);
    IncidenceR23YearN(:,~any(IncidenceR23YearN,1))=[];  %columns
    IncidenceR23YearN=sort(IncidenceR23YearN(:,:),2);
    IncR_spN(:,~any(IncR_spN,1))=[];  %columns
    IncR_spN=sort(IncR_spN(:,:),2);
end

for n=1:1 
for t=1:YearEND
    PrevalenceYear1PPPCell{t,n}(:)=PrevalenceYear1N(t,:,n);
    PrevalenceYear1PPCell{t,n}(:)=sort(PrevalenceYear1(t,:,n)); 
    PrevalenceYear1forCI{t,n}(:)=PrevalenceYear1PPPCell{t,n}(round(0.025*length(PrevalenceYear1PPPCell{t,n})):round((1-0.025)*length(PrevalenceYear1PPPCell{t,n})));   
    meanPrevalenceYear1(t,n)=mean(PrevalenceYear1PPCell{t,n}(:)); 
    ci_lowerPrevalenceYear1(t,n)=PrevalenceYear1forCI{t,n}(1);   
    ci_upperPrevalenceYear1(t,n)=PrevalenceYear1forCI{t,n}(end);  

    PrevalenceYear1IDUPPPCell{t,n}(:)=sort(PrevalenceYear1IDUN(t,:,n)); 
    PrevalenceYear1IDUPPCell{t,n}(:)=sort(PrevalenceYear1IDU(t,:,n));  
    PrevalenceYear1IDUforCI{t,n}(:)=PrevalenceYear1IDUPPPCell{t,n}(round(0.025*length(PrevalenceYear1IDUPPPCell{t,n})):round((1-0.025)*length(PrevalenceYear1IDUPPPCell{t,n})));   
    meanPrevalenceYear1IDU(t,n)=mean(PrevalenceYear1IDUPPCell{t,n}(:)); 
    ci_lowerPrevalenceYear1IDU(t,n)=PrevalenceYear1IDUforCI{t,n}(1);   
    ci_upperPrevalenceYear1IDU(t,n)=PrevalenceYear1IDUforCI{t,n}(end); 

    PrevalenceYear23PPPCell{t,n}(:)=PrevalenceYear23N(t,:,n);
    PrevalenceYear23PPCell{t,n}(:)=sort(PrevalenceYear23(t,:,n));  
    PrevalenceYear23forCI{t,n}(:)=PrevalenceYear23PPPCell{t,n}(round(0.025*length(PrevalenceYear23PPPCell{t,n})):round((1-0.025)*length(PrevalenceYear23PPPCell{t,n})));  
    meanPrevalenceYear23(t,n)=mean(PrevalenceYear23PPCell{t,n}(:)); 
    ci_lowerPrevalenceYear23(t,n)=PrevalenceYear23forCI{t,n}(1);   
    ci_upperPrevalenceYear23(t,n)=PrevalenceYear23forCI{t,n}(end); 

    Incidence1YearPPPCell{t,n}(:)=Incidence1YearN(t,:,n);
    Incidence1YearPPCell{t,n}(:)=sort(Incidence1Year(t,:,n));  
    Incidence1YearforCI{t,n}(:)=Incidence1YearPPPCell{t,n}(round(0.025*length(Incidence1YearPPPCell{t,n})):round((1-0.025)*length(Incidence1YearPPPCell{t,n})));   
    meanIncidence1Year(t,n)=mean(Incidence1YearPPCell{t,n}(:)); 
    ci_lowerIncidence1Year(t,n)=Incidence1YearforCI{t,n}(1); 
    ci_upperIncidence1Year(t,n)=Incidence1YearforCI{t,n}(end); 
    
    Incidence1YearIDUPPPCell{t,n}(:)=sort(Incidence1YearIDUN(t,:,n)); 
    Incidence1YearIDUPPCell{t,n}(:)=sort(Incidence1YearIDU(t,:,n)); 
    Incidence1YearIDUforCI{t,n}(:)=Incidence1YearIDUPPPCell{t,n}(round(0.025*length(Incidence1YearIDUPPPCell{t,n})):round((1-0.025)*length(Incidence1YearIDUPPPCell{t,n})));  
    meanIncidence1YearIDU(t,n)=mean(Incidence1YearIDUPPCell{t,n}(:)); 
    ci_lowerIncidence1YearIDU(t,n)=Incidence1YearIDUforCI{t,n}(1);   
    ci_upperIncidence1YearIDU(t,n)=Incidence1YearIDUforCI{t,n}(end);  
    
    Incidence23YearSPPPCell{t,n}(:)=Incidence23YearSN(t,:,n);
    Incidence23YearSPPCell{t,n}(:)=sort(Incidence23YearS(t,:,n));
    Incidence23YearSforCI{t,n}(:)=Incidence23YearSPPPCell{t,n}(round(0.025*length(Incidence23YearSPPPCell{t,n})):round((1-0.025)*length(Incidence23YearSPPPCell{t,n})));
    meanIncidence23YearS(t,n)=mean(Incidence23YearSPPCell{t,n}(:));
    ci_lowerIncidence23YearS(t,n)=Incidence23YearSforCI{t,n}(1);
    ci_upperIncidence23YearS(t,n)=Incidence23YearSforCI{t,n}(end);
    
    Inc_spPPPCell{t,n}(:)=Inc_spN(t,:,n);
    Inc_spPPCell{t,n}(:)=sort(Inc_sp(t,:,n));
    Inc_spPPforCI{t,n}(:)=Inc_spPPPCell{t,n}(round(0.025*length(Inc_spPPPCell{t,n})):round((1-0.025)*length(Inc_spPPPCell{t,n})));
    meanInc_sp(t,n)=mean(Inc_spPPCell{t,n}(:));
    ci_lowerInc_sp(t,n)=Inc_spPPforCI{t,n}(1);
    ci_upperInc_sp(t,n)=Inc_spPPforCI{t,n}(end);

    IncidenceR1YearPPPCell{t,n}(:)=IncidenceR1YearN(t,:,n);
    IncidenceR1YearPPCell{t,n}(:)=sort(IncidenceR1Year(t,:,n));
    IncidenceR1YearforCI{t,n}(:)=IncidenceR1YearPPPCell{t,n}(round(0.025*length(IncidenceR1YearPPPCell{t,n})):round((1-0.025)*length(IncidenceR1YearPPPCell{t,n})));
    meanIncidenceR1Year(t,n)=mean(IncidenceR1YearPPCell{t,n}(:));
    ci_lowerIncidenceR1Year(t,n)=IncidenceR1YearforCI{t,n}(1);
    ci_upperIncidenceR1Year(t,n)=IncidenceR1YearforCI{t,n}(end);
    
    IncidenceR1YearIDUPPPCell{t,n}(:)=sort(IncidenceR1YearIDUN(t,:,n));  
    IncidenceR1YearIDUPPCell{t,n}(:)=sort(IncidenceR1YearIDU(t,:,n));
    IncidenceR1YearIDUforCI{t,n}(:)=IncidenceR1YearIDUPPPCell{t,n}(round(0.025*length(IncidenceR1YearIDUPPPCell{t,n})):round((1-0.025)*length(IncidenceR1YearIDUPPPCell{t,n})));
    meanIncidenceR1YearIDU(t,n)=mean(IncidenceR1YearIDUPPCell{t,n}(:));
    ci_lowerIncidenceR1YearIDU(t,n)=IncidenceR1YearIDUforCI{t,n}(1);
    ci_upperIncidenceR1YearIDU(t,n)=IncidenceR1YearIDUforCI{t,n}(end);
    
    IncidenceR23YearPPPCell{t,n}(:)=IncidenceR23YearN(t,:,n);
    IncidenceR23YearPPCell{t,n}(:)=sort(IncidenceR23Year(t,:,n));
    IncidenceR23YearforCI{t,n}(:)=IncidenceR23YearPPPCell{t,n}(round(0.025*length(IncidenceR23YearPPPCell{t,n})):round((1-0.025)*length(IncidenceR23YearPPPCell{t,n})));
    meanIncidenceR23Year(t,n)=mean(IncidenceR23YearPPCell{t,n}(:));
    ci_lowerIncidenceR23Year(t,n)=IncidenceR23YearforCI{t,n}(1);
    ci_upperIncidenceR23Year(t,n)=IncidenceR23YearforCI{t,n}(end);
   
    IncR_spPPPCell{t,n}(:)=IncR_spN(t,:,n);
    IncR_spPPCell{t,n}(:)=sort(IncR_sp(t,:,n));
    IncR_spforCI{t,n}(:)=IncR_spPPPCell{t,n}(round(0.025*length(IncR_spPPPCell{t,n})):round((1-0.025)*length(IncR_spPPPCell{t,n})));
    meanIncR_sp(t,n)=mean(IncR_spPPCell{t,n}(:));
    ci_lowerIncR_sp(t,n)=IncR_spforCI{t,n}(1);
    ci_upperIncR_sp(t,n)=IncR_spforCI{t,n}(end);
end
end


%Estimates at baseline

meanPrevalenceYear1(2020-t0+1,1)
ci_lowerPrevalenceYear1(2020-t0+1,1)
ci_upperPrevalenceYear1(2020-t0+1,1)

meanPrevalenceYear1IDU(2020-t0+1,1)
ci_lowerPrevalenceYear1IDU(2020-t0+1,1)
ci_upperPrevalenceYear1IDU(2020-t0+1,1)

meanPrevalenceYear23(2020-t0+1,1)
ci_lowerPrevalenceYear23(2020-t0+1,1)
ci_upperPrevalenceYear23(2020-t0+1,1)

meanIncidence1Year(2020-t0+1,1)
ci_lowerIncidence1Year(2020-t0+1,1)
ci_upperIncidence1Year(2020-t0+1,1)

A=meanIncidence1Year(2020-t0+1,1)*91500/600
AL=ci_lowerIncidence1Year(2020-t0+1,1)*91500/600
AU=ci_upperIncidence1Year(2020-t0+1,1)*91500/600

meanIncidence1YearIDU(2020-t0+1,1)
ci_lowerIncidence1YearIDU(2020-t0+1,1)
ci_upperIncidence1YearIDU(2020-t0+1,1)

B=meanIncidence1YearIDU(2020-t0+1,1)*91500/600
BL=ci_lowerIncidence1YearIDU(2020-t0+1,1)*91500/600
BU=ci_upperIncidence1YearIDU(2020-t0+1,1)*91500/600

meanIncidence23YearS(2020-t0+1,1)
ci_lowerIncidence23YearS(2020-t0+1,1)
ci_upperIncidence23YearS(2020-t0+1,1)

C=meanIncidence23YearS(2020-t0+1,1)*915000/6000
CL=ci_lowerIncidence23YearS(2020-t0+1,1)*915000/6000
CU=ci_upperIncidence23YearS(2020-t0+1,1)*915000/6000

FactorART*meanInc_sp(2020-t0+1,1)
FactorART*ci_lowerInc_sp(2020-t0+1,1)
FactorART*ci_upperInc_sp(2020-t0+1,1)

D=FactorART*meanInc_sp(2020-t0+1,1)*915000/6000
DL=FactorART*ci_lowerInc_sp(2020-t0+1,1)*915000/6000
DU=FactorART*ci_upperInc_sp(2020-t0+1,1)*915000/6000

E=meanIncidenceR1Year(2020-t0+1,1)*10
EL=ci_lowerIncidenceR1Year(2020-t0+1,1)*10
EU=ci_upperIncidenceR1Year(2020-t0+1,1)*10

F=meanIncidenceR1YearIDU(2020-t0+1,1)*10
FL=ci_lowerIncidenceR1YearIDU(2020-t0+1,1)*10
FU=ci_upperIncidenceR1YearIDU(2020-t0+1,1)*10

G=meanIncidenceR23Year(2020-t0+1,1)*10
GL=ci_lowerIncidenceR23Year(2020-t0+1,1)*10
GU=ci_upperIncidenceR23Year(2020-t0+1,1)*10

H=FactorART*meanIncR_sp(2020-t0+1,1)*10
HL=FactorART*ci_lowerIncR_sp(2020-t0+1,1)*10
HU=FactorART*ci_upperIncR_sp(2020-t0+1,1)*10
    