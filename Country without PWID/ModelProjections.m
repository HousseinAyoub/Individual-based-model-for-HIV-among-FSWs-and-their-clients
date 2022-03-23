   clear all
    clc
   
 eART=0.57;
 DataCovARTClient=0.57;
 FactorART=((1-DataCovARTClient)+(1-eART)*DataCovARTClient);  
    
     load('job32_output.mat','job32_output') 
    MoroccoScenario500=job32_output;
    
    PrevalenceYear1=MoroccoScenario500{1,3};
    PrevalenceYear23=MoroccoScenario500{1,6};
    
    IncidenceR1Year=MoroccoScenario500{1,13};
    IncidenceR23Year=MoroccoScenario500{1,16};

    Incidence1YearS=MoroccoScenario500{1,23};
    Incidence1Year=MoroccoScenario500{1,24};
    Incidence23YearS=MoroccoScenario500{1,27};
    IncidenceALL=MoroccoScenario500{1,29};

    Inc_sp=MoroccoScenario500{1,38};
    IncR_sp=MoroccoScenario500{1,39};
    
    Inc_spr=MoroccoScenario500{1,44};
    IncR_spr=MoroccoScenario500{1,42};
    
    Inc_spnr=MoroccoScenario500{1,45};
    IncR_spnr=MoroccoScenario500{1,43};

    Susc=MoroccoScenario500{1,50};
    SuscClient=MoroccoScenario500{1,53};
    Susc_sp=MoroccoScenario500{1,76};
    

NumberofScenarios=1; %we simulate here only the baseline, but his can be changed to simulate the other scenarios
YearEND=400; %this is 2050
YearBURNIN=50;
DT      = 1/12; % definition of unit time (delta t)
Unit=12; %% months
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
       
       Incidence1YearN(:,m)=zeros(YearEND,1);
       Incidence23YearSN(:,m)=zeros(YearEND,1);
       Inc_spN(:,m)=zeros(YearEND,1);
       
       IncidenceR1YearN(:,m)=zeros(YearEND,1);
       IncidenceR23YearN(:,m)=zeros(YearEND,1);
       IncR_spN(:,m)=zeros(YearEND,1);
    else
       PrevalenceYear1N(:,m)=PrevalenceYear1(:,m,1);
       PrevalenceYear23N(:,m)=PrevalenceYear23(:,m,1);
       
       Incidence1YearN(:,m)=Incidence1Year(:,m,1);
       Incidence23YearSN(:,m)=Incidence23YearS(:,m,1);
       Inc_spN(:,m)=Inc_sp(:,m,1);
       
       IncidenceR1YearN(:,m)=IncidenceR1Year(:,m,1);
       IncidenceR23YearN(:,m)=IncidenceR23Year(:,m,1);
       IncR_spN(:,m)=IncR_sp(:,m,1);
    end
  end
    PrevalenceYear1N(:,~any(PrevalenceYear1N,1))=[];  %columns
    PrevalenceYear1N=sort(PrevalenceYear1N(:,:),2);
    PrevalenceYear23N(:,~any(PrevalenceYear23N,1))=[];  %columns
    PrevalenceYear23N=sort(PrevalenceYear23N(:,:),2);
    Incidence1YearN(:,~any(Incidence1YearN,1))=[];  %columns
    Incidence1YearN=sort(Incidence1YearN(:,:),2);
    Incidence23YearSN(:,~any(Incidence23YearSN,1))=[];  %columns
    Incidence23YearSN=sort(Incidence23YearSN(:,:),2);
    Inc_spN(:,~any(Inc_spN,1))=[];  %columns            
    Inc_spN=sort(Inc_spN(:,:),2);
    IncidenceR1YearN(:,~any(IncidenceR1YearN,1))=[];  %columns
    IncidenceR1YearN=sort(IncidenceR1YearN(:,:),2);
    IncidenceR23YearN(:,~any(IncidenceR23YearN,1))=[];  %columns
    IncidenceR23YearN=sort(IncidenceR23YearN(:,:),2);
    IncR_spN(:,~any(IncR_spN,1))=[];  %columns
    IncR_spN=sort(IncR_spN(:,:),2);
end


    
for n=1:1 
for t=1:YearEND
    PrevalenceYear1PPPCell{t,n}(:)=PrevalenceYear1N(t,:);
    PrevalenceYear1PPCell{t,n}(:)=sort(PrevalenceYear1(t,:,n));  
    PrevalenceYear1forCI{t,n}(:)=PrevalenceYear1PPPCell{t,n}(round(0.025*length(PrevalenceYear1PPPCell{t,n})):round((1-0.025)*length(PrevalenceYear1PPPCell{t,n})));   
    meanPrevalenceYear1(t,n)=mean(PrevalenceYear1PPCell{t,n}(:)); 
    ci_lowerPrevalenceYear1(t,n)=PrevalenceYear1forCI{t,n}(1);   
    ci_upperPrevalenceYear1(t,n)=PrevalenceYear1forCI{t,n}(end);  

    
    PrevalenceYear23PPPCell{t,n}(:)=PrevalenceYear23N(t,:);
    PrevalenceYear23PPCell{t,n}(:)=sort(PrevalenceYear23(t,:,n));  
    PrevalenceYear23forCI{t,n}(:)=PrevalenceYear23PPPCell{t,n}(round(0.025*length(PrevalenceYear23PPPCell{t,n})):round((1-0.025)*length(PrevalenceYear23PPPCell{t,n})));   
    meanPrevalenceYear23(t,n)=mean(PrevalenceYear23PPCell{t,n}(:)); 
    ci_lowerPrevalenceYear23(t,n)=PrevalenceYear23forCI{t,n}(1);   
    ci_upperPrevalenceYear23(t,n)=PrevalenceYear23forCI{t,n}(end); 

    Incidence1YearPPPCell{t,n}(:)=Incidence1YearN(t,:);
    Incidence1YearPPCell{t,n}(:)=sort(Incidence1Year(t,:,n));  
    Incidence1YearforCI{t,n}(:)=Incidence1YearPPPCell{t,n}(round(0.025*length(Incidence1YearPPPCell{t,n})):round((1-0.025)*length(Incidence1YearPPPCell{t,n})));   
    meanIncidence1Year(t,n)=mean(Incidence1YearPPCell{t,n}(:)); 
    ci_lowerIncidence1Year(t,n)=Incidence1YearforCI{t,n}(1);   
    ci_upperIncidence1Year(t,n)=Incidence1YearforCI{t,n}(end);  
    
    Incidence23YearSPPPCell{t,n}(:)=Incidence23YearSN(t,:);
    Incidence23YearSPPCell{t,n}(:)=sort(Incidence23YearS(t,:,n));
    Incidence23YearSforCI{t,n}(:)=Incidence23YearSPPPCell{t,n}(round(0.025*length(Incidence23YearSPPPCell{t,n})):round((1-0.025)*length(Incidence23YearSPPPCell{t,n})));
    meanIncidence23YearS(t,n)=mean(Incidence23YearSPPCell{t,n}(:));
    ci_lowerIncidence23YearS(t,n)=Incidence23YearSforCI{t,n}(1);
    ci_upperIncidence23YearS(t,n)=Incidence23YearSforCI{t,n}(end);
    
    Inc_spPPPCell{t,n}(:)=Inc_spN(t,:);
    Inc_spPPCell{t,n}(:)=sort(Inc_sp(t,:,n));
    Inc_spPPforCI{t,n}(:)=Inc_spPPPCell{t,n}(round(0.025*length(Inc_spPPPCell{t,n})):round((1-0.025)*length(Inc_spPPPCell{t,n})));
    meanInc_sp(t,n)=mean(Inc_spPPCell{t,n}(:));
    ci_lowerInc_sp(t,n)=Inc_spPPforCI{t,n}(1);
    ci_upperInc_sp(t,n)=Inc_spPPforCI{t,n}(end);

    IncidenceR1YearPPPCell{t,n}(:)=IncidenceR1YearN(t,:);
    IncidenceR1YearPPCell{t,n}(:)=sort(IncidenceR1Year(t,:,n));
    IncidenceR1YearforCI{t,n}(:)=IncidenceR1YearPPPCell{t,n}(round(0.025*length(IncidenceR1YearPPPCell{t,n})):round((1-0.025)*length(IncidenceR1YearPPPCell{t,n})));
    meanIncidenceR1Year(t,n)=mean(IncidenceR1YearPPCell{t,n}(:));
    ci_lowerIncidenceR1Year(t,n)=IncidenceR1YearforCI{t,n}(1);
    ci_upperIncidenceR1Year(t,n)=IncidenceR1YearforCI{t,n}(end);
    
    IncidenceR23YearPPPCell{t,n}(:)=IncidenceR23YearN(t,:);
    IncidenceR23YearPPCell{t,n}(:)=sort(IncidenceR23Year(t,:,n));
    IncidenceR23YearforCI{t,n}(:)=IncidenceR23YearPPPCell{t,n}(round(0.025*length(IncidenceR23YearPPPCell{t,n})):round((1-0.025)*length(IncidenceR23YearPPPCell{t,n})));
    meanIncidenceR23Year(t,n)=mean(IncidenceR23YearPPCell{t,n}(:));
    ci_lowerIncidenceR23Year(t,n)=IncidenceR23YearforCI{t,n}(1);
    ci_upperIncidenceR23Year(t,n)=IncidenceR23YearforCI{t,n}(end);
   
    IncR_spPPPCell{t,n}(:)=IncR_spN(t,:);
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

meanPrevalenceYear23(2020-t0+1,1)
ci_lowerPrevalenceYear23(2020-t0+1,1)
ci_upperPrevalenceYear23(2020-t0+1,1)

meanIncidence1Year(2020-t0+1,1)
ci_lowerIncidence1Year(2020-t0+1,1)
ci_upperIncidence1Year(2020-t0+1,1)

A=meanIncidence1Year(2020-t0+1,1)*72000/600
AL=ci_lowerIncidence1Year(2020-t0+1,1)*72000/600
AU=ci_upperIncidence1Year(2020-t0+1,1)*72000/600

meanIncidence23YearS(2020-t0+1,1)
ci_lowerIncidence23YearS(2020-t0+1,1)
ci_upperIncidence23YearS(2020-t0+1,1)

B=meanIncidence23YearS(2020-t0+1,1)*720000/6000
BL=ci_lowerIncidence23YearS(2020-t0+1,1)*720000/6000
BU=ci_upperIncidence23YearS(2020-t0+1,1)*720000/6000

FactorART*meanInc_sp(2020-t0+1,1)
FactorART*ci_lowerInc_sp(2020-t0+1,1)
FactorART*ci_upperInc_sp(2020-t0+1,1)

C=FactorART*meanInc_sp(2020-t0+1,1)*720000/6000
CL=FactorART*ci_lowerInc_sp(2020-t0+1,1)*720000/6000
CU=FactorART*ci_upperInc_sp(2020-t0+1,1)*720000/6000

E=meanIncidenceR1Year(2020-t0+1,1)*10
EL=ci_lowerIncidenceR1Year(2020-t0+1,1)*10
EU=ci_upperIncidenceR1Year(2020-t0+1,1)*10

F=meanIncidenceR23Year(2020-t0+1,1)*10
FL=ci_lowerIncidenceR23Year(2020-t0+1,1)*10
FU=ci_upperIncidenceR23Year(2020-t0+1,1)*10

G=FactorART*meanIncR_sp(2020-t0+1,1)*10
GL=FactorART*ci_lowerIncR_sp(2020-t0+1,1)*10
GU=FactorART*ci_upperIncR_sp(2020-t0+1,1)*10


