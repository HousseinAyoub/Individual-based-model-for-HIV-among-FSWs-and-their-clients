
function [MeannumberofContactreg,MeannumberofContactnreg,PrevalenceYear1,PrevalenceYear2,PrevalenceYear3,PrevalenceYear23,PrevalenceYear4,PrevalenceYear1_Avg,PrevalenceYear2_Avg,PrevalenceYear3_Avg,PrevalenceYear23_Avg,PrevalenceYear4_Avg,IncidenceR1Year,IncidenceR2Year,IncidenceR3Year,IncidenceR23Year,IncidenceR4Year,IncidenceR1Year_Avg,IncidenceR2Year_Avg,IncidenceR3Year_Avg,IncidenceR23Year_Avg,IncidenceR4Year_Avg,Incidence1YearS,Incidence1Year,Incidence2YearS,Incidence3YearS,Incidence23YearS,Incidence4YearS,IncidenceALL,Incidence1YearS_Avg,Incidence1Year_Avg,Incidence2YearS_Avg,Incidence3YearS_Avg,Incidence23YearS_Avg,Incidence4YearS_Avg,IncidenceALL_Avg,Number_I_FSW,Inc_sp,IncR_sp,Inc_sp_Avg,IncR_sp_Avg,IncR_spr,IncR_spnr,Inc_spr,Inc_spnr,IncR_spr_Avg,IncR_spnr_Avg,Inc_spr_Avg,Inc_spnr_Avg,Susc,SuscRClient,SuscNRClient,SuscClient,Susc_Avg,SuscRClient_Avg,SuscNRClient_Avg,SuscClient_Avg,InfectedFSW,InfectedRClients,InfectedNRClients,InfectedClients,InfectedALL,InfectedFSW_Avg,InfectedRClients_Avg,InfectedNRClients_Avg,InfectedClients_Avg,InfectedALL_Avg,Ndisc,Ndiscr,Ndiscnr,Ndisc_Avg,Ndiscr_Avg,Ndiscnr_Avg,Susc_spr,Susc_spnr,Susc_sp,Susc_spr_Avg,Susc_spnr_Avg,Susc_sp_Avg]=Model_Interv_IDU(DT,Unit,fsw,regm,nregm,tpHIVreg,tpHIVn,transHIV,removelink_reg,removelink_nreg,YearEND,YearBURNIN,BURN_IN,TEND,INIT_INF_fsw,INIT_INF_regm,INIT_INF_nregm,CoverageMC,CoverageCond,CovARTFSW,CovARTClient,CovPrEPFSW,CovPrEPClient,eARTT,eMC,eCond,ePrEP,PAR_reg1, PAR_reg2,PAR_nreg1, PAR_nreg2,mu,muHIV,npartr,npartnr,npart,nact_sp,tpHIVspouse,Cond_sp,tpsp)


%% initialization
acq_regm    = zeros(1,regm); % acquistion probability of new partner for regular clients
acq_nregm   = zeros(1,nregm); % acquistion probability of new partner for non-regular clients
    for cnt = 1:regm
        acq_regm(cnt) = gamrnd(PAR_reg1, PAR_reg2);
    end

    for cnt = 1:nregm
        acq_nregm(cnt) = gamrnd(PAR_nreg1, PAR_nreg2);
    end


%% prepare array and matrix
adjmatrix   = zeros(fsw+regm+nregm); % adjacency matrix (who connects with whom, 0: not connected, 1: connected)
   
state_ART_regm=zeros(1,regm); % ART status
state_ART_nregm=zeros(1,nregm); % ART status
state_ART_fsw=zeros(1,fsw); % ART status

state_fsw   = zeros(1,fsw); %HIV state for FSW (state varies from 0 (Suseptible) to 3 (advanced) for each person
state_regm  = zeros(1,regm); %HIV state for reg
state_nregm = zeros(1,nregm); %HIV state for nonreg

%%%%%%%%%%%%%%%%%%%%Injecting drug use (Not applicable)%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fsw_PWID=0;
%   fsw_notPWID=0;
%   IDU_fsw   = zeros(1,fsw); %IDU state for FSW (state varies from 0 (non-injecting) 1 (injecting) for each person
% 
%   for cnt= 1:fsw
%           if (rand<IDUF)
%               IDU_fsw(cnt)=1;
%               fsw_PWID = fsw_PWID +1;
%           else 
%               fsw_notPWID=fsw_notPWID +1;
%           end
%   end    

%% main rourtine

for t = 1:TEND
    if mod(t,12)==0  
        disp(t)
    end
    %%%% FSWs death
    for cnt = 1:fsw 
        % by AIDS and by non-AIDS
       if state_ART_fsw(cnt) == 0
          if ((state_fsw(cnt) == 3) && (rand < muHIV))||(rand < mu) % if a FSW is advance stage
                state_fsw(cnt)= 0; % for constant population size
                adjmatrix(cnt,1:fsw+regm+nregm) = zeros(1,fsw+regm+nregm); % reset all paths in network involved by her
                adjmatrix(1:(fsw+regm+nregm),cnt) = zeros((fsw+regm+nregm),1); % reset all paths in network involved by her
          end
       else
          if ((state_fsw(cnt) == 3) && (rand < muHIV/3))||(rand < mu) % if a FSW is advance stage
                state_fsw(cnt)= 0; % for constant population size
                adjmatrix(cnt,1:fsw+regm+nregm) = zeros(1,fsw+regm+nregm); % reset all paths in network involved by her
                adjmatrix(1:(fsw+regm+nregm),cnt) = zeros((fsw+regm+nregm),1); % reset all paths in network involved by her
          end
       end
    end
    
    %%%% for regular clients death
    for cnt = 1:regm 
        % by AIDS and by non-AIDS 
       if state_ART_regm(cnt) == 0
          if ((state_regm(cnt) == 3) && (rand < muHIV))||(rand < mu) % if a regular client is advance stage
                state_regm(cnt)= 0; % for constant population size
                adjmatrix(cnt+fsw,1:(fsw+regm+nregm))= zeros(1,(fsw+regm+nregm)); % reset all paths in network involved by him
                adjmatrix(1:(fsw+regm+nregm),cnt+fsw)= zeros((fsw+regm+nregm),1); % reset all paths in network involved by him  
          end
       else
          if ((state_regm(cnt) == 3) && (rand < muHIV/3))||(rand < mu) % if a regular client is advance stage
                state_regm(cnt)= 0; % for constant population size
                adjmatrix(cnt+fsw,1:(fsw+regm+nregm))= zeros(1,(fsw+regm+nregm)); % reset all paths in network involved by him
                adjmatrix(1:(fsw+regm+nregm),cnt+fsw)= zeros((fsw+regm+nregm),1); % reset all paths in network involved by him  
          end
       end
    end
    
    % for non-regular clients
    fswregm=fsw+regm;
    for cnt = 1:nregm 
        % by AIDS and by non-AIDS 
       if state_ART_nregm(cnt) == 0
          if ((state_nregm(cnt) == 3) && (rand < muHIV))||(rand < mu)  % if a regular client is advance stage
               state_nregm(cnt)= 0; % for constant population size
               adjmatrix(cnt+fsw+regm,1:(fsw+regm+nregm))= zeros(1,(fsw+regm+nregm)); % reset all paths in network involved by him
               adjmatrix(1:(fsw+regm+nregm),cnt+fsw+regm)= zeros((fsw+regm+nregm),1); % reset all paths in network involved by him
          end
       else
          if ((state_nregm(cnt) == 3) && (rand < muHIV/3))||(rand < mu)  % if a regular client is advance stage
               state_nregm(cnt)= 0; % for constant population size
               adjmatrix(cnt+fsw+regm,1:(fsw+regm+nregm))= zeros(1,(fsw+regm+nregm)); % reset all paths in network involved by him
               adjmatrix(1:(fsw+regm+nregm),cnt+fsw+regm)= zeros((fsw+regm+nregm),1); % reset all paths in network involved by him
          end
       end
    end
 
    %% dissolution of partnership (we need to focus on only males, dissolution connetion by male is equivalent with that by female)
    % for regular client
    for cnt = 1:regm % stochastic dissolution for regular client partnership
        for cnt2 = 1:fsw
            if (adjmatrix(cnt2,cnt+fsw) == 1)&& (rand < removelink_reg)
                    adjmatrix(cnt2,cnt+fsw) = 0;
                    adjmatrix(cnt+fsw,cnt2) = 0;
            end
        end
    end
    % for non-regular client
    for cnt = 1:nregm % dissolution for non-regular client partnership
        for cnt2 = 1:fsw
            if (adjmatrix(cnt2,cnt+fswregm) == 1)&& (rand < removelink_nreg)
                    adjmatrix(cnt2,cnt+fswregm) = 0;
                    adjmatrix(cnt+fsw+regm,cnt2) = 0;
            end
        end
    end

    %% forming new partnership
    % forming partnership for regular client
    for cnt = 1:regm % stochastic dissolution for regular client partnership
        for cnt2 = 1:fsw
            if (adjmatrix(cnt2,cnt+fsw) == 0)&& (rand <(acq_regm(cnt)/fsw)) %fsw regm
                    adjmatrix(cnt2,cnt+fsw) = 1;
                    adjmatrix(cnt+fsw,cnt2) = 1;
            end
        end
    end
    % forming partnership for non-regular client
    for cnt = 1:nregm
        for cnt2 = 1:fsw
            if (adjmatrix(cnt2,cnt+fswregm) == 0)&& (rand < (acq_nregm(cnt)/fsw))%fsw nregm
                    adjmatrix(cnt2,cnt+fswregm) = 1;
                    adjmatrix(cnt+fsw+regm,cnt2) = 1;
            end
        end
    end

k=0;
    for i=1:fsw
        for j=fsw+1:regm+fsw
            if adjmatrix(i,j)==1
                k=k+1;
            end
        end
    end
    MeannumberofContactreg(t)=k;
    
h=0;
     for i=1:fsw
        for j=fsw+regm+1:nregm+regm+fsw
            if adjmatrix(i,j)==1
                h=h+1;
            end
        end
    end
    MeannumberofContactnreg(t)=h;
    %% introduction of HIV into the population at t= BURN_IN
    if t == BURN_IN
        for cnt = 1:INIT_INF_fsw
            tmp = rand*(fsw); % for fixed population size
            state_fsw(ceil(tmp)) = 2;
        end
        for cnt = 1:INIT_INF_regm
            tmp = rand*(regm); % for fixed population size
            state_regm(ceil(tmp)) = 2;
        end
        for cnt = 1:INIT_INF_nregm
            tmp = rand*(nregm); % for fixed population size
            state_nregm(ceil(tmp)) = 2;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% update HIV infection state (starts after "burn-in" to create stable network)
    %infection
        incidence_fsw = 0;      % reset counter for incidence for sexual transmission in FSW
        incidence_regm = 0;     % reset counter for incidence for sexual transmission in regular client
        incidence_nregm = 0;    % reset counter for incidence for sexual transmission in non-regular client
%         incidence_fsw_IDU = 0;    % reset counter for incidence for IDU in FSW
%         incidence_fsw_notIDU =0;  % reset counter for incidence in FSW who are not IDU 

    if t > BURN_IN
        
%%%%%%%%%%%%%%%%%%% HIV transmission due IDUs%%%%%%%%%%%%%%%%%%%%%%%%HIV among PWID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%   for cnt= 1:fsw
%       if (IDU_fsw(cnt)==1) && (state_fsw(cnt)==0) 
%           if rand<HIVPWID
%              state_fsw(cnt)=1;
%              incidence_fsw_IDU = incidence_fsw_IDU + 1;
%           end
%       end
%   end 

%%%intervention
eART=eARTT(t);

%I renamed this coverages for technical (Parallel) purpose
CovMC=CoverageMC(t); %This is for MC coverage
CovCond=CoverageCond(t);  %This is for Cond coverage
CovARTF=CovARTFSW(t); %This is for ART coverage in FSWs
CovARTM=CovARTClient(t); %This is for ART coverage in Clients
CovPrEPF=CovPrEPFSW(t); %This is for PrEP coverage in FSWs
CovPrEPM=CovPrEPClient(t); %This is for PrEP coverage in Clients

        %% HIV trasmission due sexual acts
        % Trasmission from FSW to regular and non-regular 
        for cnt = 1:fsw  % look all FSWs
            if state_fsw(cnt) ~= 0 % if a FSW is infected
               %%%I renamed this for technical (Parallel) purpose 
               tpreg=tpHIVreg(state_fsw(cnt));
               tpnreg=tpHIVn(state_fsw(cnt)); 
                if rand<CovARTF  % if a FSW had ART 
                   state_ART_fsw(cnt)=1;
                   %%%% Transmission to Regular 
                   for cnt2 = 1:regm % look all Reg
                        if rand<CovMC  % if he is circumcised
                            %state_MC_regm(cnt2)=1;
                            if rand<CovCond % if he had condom
                               %state_Cond_regm(cnt2)=1;
                                if rand<CovPrEPM % if he is on PrEP (ART+ & MC+ & Cond+ & PrEP+)
                                   %state_Prep_regm(cnt2)=1;
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*(1-eART)*(1-ePrEP)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                else %(ART+ & MC+ & Cond+ & PrEP-)
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eART)*(1-eCond)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                end
                            else %if he did not have condoms (ART+ & MC+ & Cond- & PrEP+)
                                if rand<CovPrEPM % if he did not have condoms but is on PrEP
                                   %state_Prep_regm(cnt2)=1; 
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eART)*(1-ePrEP)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                else %(ART+ & MC+ & Cond- & PrEP-)
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eART)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                end
                            end
                       else % if he is not circumcised             
                            if rand<CovCond %when he is not cicumcised but uses condoms (ART+ & MC- & Cond+ & PrEP+)
                               %state_Cond_regm(cnt2)=1;
                                if rand<CovPrEPM
                                    %state_Prep_regm(cnt2)=1;
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eCond)*(1-eART)*(1-ePrEP)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                else %(ART+ & MC- & Cond+ & PrEP-)
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eCond)*(1-eART)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                end
                            else %when he is not cicumcised and does not use condoms (ART+ & MC- & Cond- & PrEP+)
                                if rand<CovPrEPM
                                    %state_Prep_regm(cnt2)=1;
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eART)*(1-ePrEP)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                else %(ART+ & MC- & Cond- & PrEP-)
                                    if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eART)*tpreg)
                                    state_regm(cnt2) = 1;
                                    incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                    end
                                end
                            end
                        end
                    end
                   
                %%%% Transmission to non Regular 
                   for cnt2 = 1:nregm % look all Reg
                      if rand<CovMC %if non-regular client is circumcised
                         %state_MC_nregm(cnt2)=1; 
                         if rand<CovCond % if non-regular client uses condoms
                             %state_Cond_nregm(cnt2)=1;
                             if rand<CovPrEPM % if he is on PrEP (ART+ & MC+ & Cond+ & PrEP+)
                               %state_Prep_nregm(cnt2)=1;
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*(1-eART)*(1-ePrEP)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of nonregular clients
                                end
                             else %(ART+ & MC+ & Cond+ & PrEP-)
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*(1-eART)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of nonregular clients
                                end
                             end
                         else %if non-regular client is circumcised but not using condoms
                             if rand<CovPrEPM % if he is on PrEP (ART+ & MC+ & Cond- & PrEP+)
                                %state_Prep_nregm(cnt2)=1;
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eART)*(1-ePrEP)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of nonregular clients
                                end
                             else %(ART+ & MC+ & Cond- & PrEP-)
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eART)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of nonregular clients
                                end
                             end
                         end
                      else % if non-regular client is not circumcised
                        if rand<CovCond %but non-regular client uses condoms
                            %state_Cond_nregm(cnt2)=1;
                             if rand<CovPrEPM % if he is on PrEP (ART+ & MC- & Cond+ & PrEP+)
                                %state_Prep_nregm(cnt2)=1; 
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eCond)*(1-eART)*(1-ePrEP)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                                end
                             else % if he is not on PrEP (ART+ & MC- & Cond+ & PrEP-)
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eCond)*(1-eART)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                                end
                             end
                        else % if non-regular client does not use condoms
                            if rand<CovPrEPM % if he is on PrEP (ART+ & MC- & Cond- & PrEP+)
                               %state_Prep_nregm(cnt2)=1;
                               if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eART)*(1-ePrEP)*tpnreg)
                               state_nregm(cnt2) = 1;
                               incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                               end
                            else % if he is not on PrEP (ART+ & MC- & Cond- & PrEP-)
                               if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eART)*tpnreg)
                               state_nregm(cnt2) = 1;
                               incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                               end
                            end
                         end
                      end 
                   end
                else % If FSW is not on ART
                %%%% Transmission to Regular
                   for cnt2 = 1:regm % look all Reg
                      if rand<CovMC %if regular client is circumcised
                         %state_MC_regm(cnt2)=1; 
                         if rand<CovCond % if he is using comdoms
                            %state_Cond_regm(cnt2)=1;  
                             if rand<CovPrEPM % if he is on PrEP (ART- & MC+ & Cond+ & PrEP+)
                                %state_Prep_regm(cnt2)=1; 
                                if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*(1-ePrEP)*tpreg)
                                state_regm(cnt2) = 1;
                                incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                end
                             else % if he is not on PrEP (ART- & MC+ & Cond+ & PrEP-)
                                if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*tpreg)
                                state_regm(cnt2) = 1;
                                incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                end
                             end
                         else % if he is not using condoms
                             if rand<CovPrEPM % if he is on PrEP (ART- & MC+ & Cond- & PrEP+)
                                %state_Prep_regm(cnt2)=1; 
                                if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*(1-ePrEP)*tpreg)
                                state_regm(cnt2) = 1;
                                incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                end
                             else % if he is not on PrEP (ART- & MC+ & Cond- & PrEP-)
                                if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eMC)*tpreg)
                                state_regm(cnt2) = 1;
                                incidence_regm = incidence_regm + 1; % count incidence of regular clients
                                end 
                             end
                         end 
                      else %if regular client is not circumcised
                        if rand<CovCond %but regular client uses condoms
                           %state_Cond_regm(cnt2)=1;  
                           if rand<CovPrEPM % if he is on PrEP (ART- & MC- & Cond+ & PrEP+)
                              %state_Prep_regm(cnt2)=1;  
                              if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eCond)*(1-ePrEP)*tpreg)
                              state_regm(cnt2) = 1;
                              incidence_regm = incidence_regm + 1; % count incidence of regular clients
                              end
                           else % if he is not on PrEP (ART- & MC- & Cond+ & PrEP-)
                              if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-eCond)*tpreg)
                              state_regm(cnt2) = 1;
                              incidence_regm = incidence_regm + 1; % count incidence of regular clients
                              end
                           end     
                        else % if regular client is not using condoms
                            if rand<CovPrEPM % if he is on PrEP (ART- & MC- & Cond- & PrEP+)
                               %state_Prep_regm(cnt2)=1;  
                               if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < (1-ePrEP)*tpreg)
                               state_regm(cnt2) = 1;
                               incidence_regm = incidence_regm + 1; % count incidence of regular clients
                               end
                            else % if he is not on PrEP (ART- & MC- & Cond- & PrEP-)
                               if (state_regm(cnt2) == 0)&&(adjmatrix(cnt,fsw+cnt2) == 1)&& (rand < tpreg)
                               state_regm(cnt2) = 1;
                               incidence_regm = incidence_regm + 1; % count incidence of regular clients
                               end
                            end
                        end
                      end
                   end
                  %%%% Transmission to nonRegular
                   for cnt2 = 1:nregm % look all nonReg
                      if rand<CovMC %if non-regular client is not on ART but is circumcised
                         %state_MC_nregm(cnt2)=1;  
                         if rand<CovCond % if he uses condoms
                            %state_Cond_nregm(cnt2)=1; 
                             if rand<CovPrEPM % (ART- & MC+ & Cond+ & PrEP+)
                                %state_Prep_nregm(cnt2)=1; 
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*(1-ePrEP)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1;
                                end
                             else %(ART- & MC+ & Cond+ & PrEP-)
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-eCond)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1;
                                end
                             end
                         else % if non-regular client is not using condoms
                             if rand<CovPrEPM % (ART- & MC+ & Cond- & PrEP+)
                                %state_Prep_nregm(cnt2)=1;  
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*(1-ePrEP)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; 
                                end
                             else % (ART- & MC+ & Cond- & PrEP-)
                                if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eMC)*tpnreg)
                                state_nregm(cnt2) = 1;
                                incidence_nregm = incidence_nregm + 1; 
                                end
                             end
                         end
                      else % if non-regular client is not circumcised
                        if rand<CovCond %but non-regular client uses condoms
                           %state_Cond_nregm(cnt2)=1;  
                           if rand<CovPrEPM % (ART- & MC- & Cond+ & PrEP+)
                              %state_Prep_nregm(cnt2)=1;  
                              if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eCond)*(1-ePrEP)*tpnreg)
                              state_nregm(cnt2) = 1;
                              incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                              end
                           else % (ART- & MC- & Cond+ & PrEP-)
                              if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-eCond)*tpnreg)
                              state_nregm(cnt2) = 1;
                              incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                              end
                           end
                        else % if non-regular client is not using condoms
                            if rand<CovPrEPM % (ART- & MC- & Cond- & PrEP+)
                              %state_Prep_nregm(cnt2)=1;  
                               if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < (1-ePrEP)*tpnreg)
                               state_nregm(cnt2) = 1;
                               incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                               end
                            else %(ART- & MC- & Cond- & PrEP-)
                               if (state_nregm(cnt2) == 0)&&(adjmatrix(cnt,fswregm+cnt2) == 1)&& (rand < tpnreg)
                               state_nregm(cnt2) = 1;
                               incidence_nregm = incidence_nregm + 1; % count incidence of regular clients
                               end
                            end
                        end
                      end
                   end   
                   
                end
            end
        end
                           
        % trasmission from regular clients to FSW
        for cnt = 1:regm % look all regular clients
            if state_regm(cnt) ~= 0 % if a regular client is infected...
             %%%I renamed this for technical (Parallel) purpose 
             tpregFSW=tpHIVreg(state_regm(cnt));  
               if rand<CovARTM  % if this regular client had ART
               state_ART_regm(cnt)=1;  
                    if rand<CovCond % if this regular client used condoms  
                       %state_Cond_regm(cnt)=1;  
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW is on PrEP (ART+ & Cond+ &PrEP+)
                              %state_Prep_fsw(cnt2)=1; 
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eART)*(1-eCond)*(1-ePrEP)*tpregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                              end
                           else % if FSW is not on PrEP (ART+ & Cond+ &PrEP-)
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eART)*(1-eCond)*tpregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                              end
                           end
                       end
                   else % if the regular client did not use condoms
                      for cnt2 = 1:fsw % for FSWs
                          if rand<CovPrEPF %if the FSW is on PrEP (ART+ & Cond- &PrEP+)
                             %state_Prep_fsw(cnt2)=1; 
                             if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eART)*(1-ePrEP)*tpregFSW)
                             state_fsw(cnt2) = 1;
                             incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                             end
                          else %if the FSW is not on PrEP (ART+ & Cond- &PrEP-)
                             if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eART)*tpregFSW)
                             state_fsw(cnt2) = 1;
                             incidence_fsw = incidence_fsw + 1; % count incidence of fsw 
                             end
                          end
                      end
                    end
               else %if the regular client was not on ART
                    if rand<CovCond % but this regular client used condoms
                       %state_Cond_regm(cnt)=1; 
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW is on PrEP (ART- & Cond+ &PrEP+)
                              %state_Prep_fsw(cnt2)=1;  
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eCond)*(1-ePrEP)*tpregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                              end
                           else %if the FSW is not on PrEP (ART- & Cond+ &PrEP-)
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-eCond)*tpregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                              end
                           end
                       end
                    else %if the regular client did not use condoms
                        for cnt2 = 1:fsw % for FSWs
                            if rand<CovPrEPF %if the FSW is on PrEP (ART- & Cond- &PrEP+)
                               %state_Prep_fsw(cnt2)=1; 
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < (1-ePrEP)*tpregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                               end
                            else  %if the FSW is not on PrEP (ART- & Cond- &PrEP-)
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fsw+cnt,cnt2) == 1) && (rand < tpregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                               end
                            end
                        end
                    end
               end
            end
        end
                   
        % trasmission from non-regular clients to FSW
        for cnt = 1:nregm % look all regular clients
            if state_nregm(cnt) ~= 0 % if a nregular client is infected...
             %%%I renamed this for technical (Parallel) purpose 
             tpnregFSW=tpHIVn(state_nregm(cnt));     
                if rand<CovARTM  % if this nreg did not had ART
                   state_ART_nregm(cnt)=1; 
                   if rand<CovCond % if this nonreg used condoms
                      %state_Cond_nregm(cnt)=1; 
                      for cnt2 = 1:fsw % for FSWs
                          if rand<CovPrEPF %if the FSW was on PrEP (ART+ & Cond+ & PreP+)
                             %state_Prep_fsw(cnt2)=1; 
                             if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eART)*(1-eCond)*(1-ePrEP)*tpnregFSW)
                             state_fsw(cnt2) = 1;
                             incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                             end
                          else %if FSW was not on PrEP (ART+ & Cond+ & PrEP-)
                             if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eART)*(1-eCond)*tpnregFSW)
                             state_fsw(cnt2) = 1;
                             incidence_fsw = incidence_fsw + 1; % count incidence of fsw   
                             end
                          end
                       end  
                   else %if non-regular client is not using condoms
                        for cnt2 = 1:fsw % for FSWs
                            if rand<CovPrEPF %if the FSW was on PrEP (ART+ & Cond- & PreP+)
                               %state_Prep_fsw(cnt2)=1;
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eART)*(1-ePrEP)*tpnregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                               end
                            else %if the FSW was not on PrEP (ART+ & Cond- & PreP-)
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eART)*tpnregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                               end
                            end
                        end
                    end
               else %if non-regular client was not on ART 
                    if rand<CovCond % if this nreg had cond  
                       %state_Cond_nregm(cnt)=1; 
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW was on PrEP (ART- & Cond+ & PreP+)
                              %state_Prep_fsw(cnt2)=1; 
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eCond)*(1-ePrEP)*tpnregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                              end
                           else %if the FSW was not on PrEP (ART- & Cond+ & PreP-)
                              if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-eCond)*tpnregFSW)
                              state_fsw(cnt2) = 1;
                              incidence_fsw = incidence_fsw + 1; % count incidence of fsw                               
                              end
                           end
                       end
                    else % if this nreg did not use condoms
                        for cnt2 = 1:fsw % for FSWs
                            if rand<CovPrEPF %if the FSW was on PrEP (ART- & Cond- & PreP+)
                               %state_Prep_fsw(cnt2)=1; 
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < (1-ePrEP)*tpnregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw  
                               end
                            else %if the FSW was not on PrEP (ART- & Cond- & PreP-)
                               if (state_fsw(cnt2) == 0) &&(adjmatrix(fswregm+cnt,cnt2) == 1) && (rand < tpnregFSW)
                               state_fsw(cnt2) = 1;
                               incidence_fsw = incidence_fsw + 1; % count incidence of fsw
                               end
                            end
                        end
                    end
                end
            end
         end
        
        %% HIV infection state transition
        %%%I renamed this for technical (Parallel) purpose 
        transHIV1=transHIV(1);
        transHIV2=transHIV(2);
        % for FSW
        for cnt = 1:fsw
          if state_ART_fsw(cnt) == 0
            if (state_fsw(cnt) == 2) && (rand < transHIV2)
                    state_fsw(cnt) = 3;
            elseif (state_fsw(cnt) == 1) && (rand < transHIV1)
                    state_fsw(cnt) = 2;
            end
          else
            if (state_fsw(cnt) == 2) && (rand < transHIV2/3) 
                    state_fsw(cnt) = 3;
            elseif (state_fsw(cnt) == 1) && (rand < transHIV1)
                    state_fsw(cnt) = 2;
            end
          end
        end

        % for regular clients
        for cnt = 1:regm
          if state_ART_regm(cnt) == 0
            if (state_regm(cnt) == 2) && (rand < transHIV2)
                    state_regm(cnt) = 3;
            elseif (state_regm(cnt) == 1) && (rand < transHIV1)
                    state_regm(cnt) = 2;
            end
          else
            if (state_regm(cnt) == 2) && (rand < transHIV2/3) %
                    state_regm(cnt) = 3;
            elseif (state_regm(cnt) == 1) && (rand < transHIV1)
                    state_regm(cnt) = 2;
            end
          end
        end

        % for non-regular clients
        for cnt = 1:nregm
          if state_ART_nregm(cnt) == 0
            if (state_nregm(cnt) == 2) && (rand < transHIV2)
                    state_nregm(cnt) = 3;
            elseif (state_nregm(cnt) == 1) && (rand < transHIV1)
                    state_nregm(cnt) = 2;
            end
          else
            if (state_nregm(cnt) == 2) && (rand < transHIV2/3) %
                    state_nregm(cnt) = 3;
            elseif (state_nregm(cnt) == 1) && (rand < transHIV1)
                    state_nregm(cnt) = 2;
            end
          end
        end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prevalence_regm = 0;    % reset counter for prevalence of regular client
prevalence_nregm = 0;   % reset counter for prevalence of non-regular client
prevalence_fsw = 0; % reset counter for prevalence in FSWs who inject drugs

susceptible_regm = 0;    % reset counter for susceptible of regular client
susceptible_nregm = 0;   % reset counter for susceptible of non-regular client
susceptible_fsw = 0;   % reset counter for susceptible FSWs who inject drugs

    for cnt = 1:regm
        if state_regm(cnt) ~= 0  % if he is infected (not susceptible)
           prevalence_regm = prevalence_regm + 1; % count prevalence for regular clients
        else
           susceptible_regm=susceptible_regm+1; 
        end
    end
    for cnt = 1:nregm
        if state_nregm(cnt) ~= 0  % if he is infected (not susceptible)
           prevalence_nregm = prevalence_nregm + 1; % count prevalence for non-regular clients
        else
           susceptible_nregm=susceptible_nregm+1; 
        end
    end

    for cnt = 1:fsw
           if state_fsw(cnt) ~= 0
              prevalence_fsw = prevalence_fsw + 1;
           else
              susceptible_fsw = susceptible_fsw + 1;
           end
     end
        
     
    for cnt = 1:regm
        if state_regm(cnt) ~= 0  % if he is infected (not susceptible)
           prevalence_regm = prevalence_regm + 1; % count prevalence for regular clients
        else
           susceptible_regm=susceptible_regm+1; 
        end
    end
    for cnt = 1:nregm
        if state_nregm(cnt) ~= 0  % if he is infected (not susceptible)
           prevalence_nregm = prevalence_nregm + 1; % count prevalence for non-regular clients
        else
           susceptible_nregm=susceptible_nregm+1; 
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
       Number_I_FSW(t,1)=prevalence_fsw;     %Number of infected FSWs (injecting and non-injecting)
       Number_I_regm(t,2)=prevalence_regm;                             %Number of infected regular clients
       Number_I_nregm(t,3)=prevalence_nregm;                           %Number of infected non-regular clients
       Number_I_ALL(t,4)=(Number_I_FSW(t,1)+prevalence_regm+prevalence_nregm); %Total number of infected individuals in the population

       Inc_sexual(t,1)=incidence_fsw;                                  %Total number of new infections among FSWs due to sexual transmission
       Inc_sexual(t,2)=incidence_regm;                                 %Total number of new infections among regular clients due to sexual transmission
       Inc_sexual(t,3)=incidence_nregm;                                %Total number of new infections among non-regular clients due to sexual transmission
       Inc_sexual(t,4)=(incidence_fsw+incidence_regm+incidence_nregm); %Total number of new infections in the population due to sexual transmission
       Incidence(t,1)=incidence_fsw;               %Total number of new infections among FSWs
       
       Susceptible(t,1)=susceptible_fsw;    %Number of susceptible FSWs (injecting and non-injecting)
       Susceptible(t,2)=susceptible_regm;                              %Number of susceptible regular clients
       Susceptible(t,3)=susceptible_nregm;                             %Number of susceptible non-regular clients
       Susceptible(t,4)=(Susceptible(t,1)+susceptible_regm+susceptible_nregm); %Number of susceptible individuals in the population

end

for t=1:YearEND
    
    InfectedFSW(t)=Number_I_FSW((t*Unit),1);                                %Number of infected FSWs (injecting and non-injecting) at 12 months, 24 months...
     
    InfectedRClients(t)=Number_I_regm((t*Unit),2);                          %Number of infected regular clients at 12 months, 24 months...
    InfectedNRClients(t)=Number_I_nregm((t*Unit),3);                        %Number of infected non-regular clients at 12 months, 24 months...
    InfectedClients(t)=InfectedRClients(t)+InfectedNRClients(t);            %Number of infected clients at 12 months, 24 months...
    InfectedALL(t)=Number_I_ALL((t*Unit),4);                                %Number of infected FSWs (all) and clients (all) at 12 months, 24 months...
        
    Susc(t)=Susceptible((t*Unit),1);                                        %Number of susceptible FSWs (injecting and non-injecting) at 12 months, 24 months..
    SuscRClient(t)=Susceptible((t*Unit),2);                                 %Number of susceptible regular clients at 12 months, 24 months..
    SuscNRClient(t)=Susceptible((t*Unit),3);                                %Number of susceptible non-regular clients at 12 months, 24 months..
    SuscClient(t)=SuscRClient(t)+SuscNRClient(t);                           %Number of susceptible clients (regular and non-regular) at 12 months, 24 months..
    
    Incidence1YearS(t)=sum(Inc_sexual(1+(t-1)*Unit:t*Unit,1));              %Total number of new infections per year for xxx years among FSWs due to sexual transmission, that is incidence from 0:12 months, 12:24 months ....
    Incidence2YearS(t)=sum(Inc_sexual(1+(t-1)*Unit:t*Unit,2));              %Total number of new infections per year due to sexual transmission for xxx years among regular clients, that is incidence from 0:12 months, 12:24 months ....
    Incidence3YearS(t)=sum(Inc_sexual(1+(t-1)*Unit:t*Unit,3));              %Total number of new infections per year due to sexual transmission for xxx years among non-regular clients, that is incidence from 0:12 months, 12:24 months ....
    Incidence23YearS(t)=sum(Inc_sexual(1+(t-1)*Unit:t*Unit,2))+sum(Inc_sexual(1+(t-1)*Unit:t*Unit,3)); %Total number of new infections per year due to sexual transmission for xxx years among all clients, that is incidence from 0:12 months, 12:24 months ....
    Incidence4YearS(t)=sum(Inc_sexual(1+(t-1)*Unit:t*Unit,4));              %Total number of new infections per year due to sexual transmission for xxx years in the population, that is incidence from 0:12 months, 12:24 months ....
    
    Incidence1Year(t)=sum(Incidence(1+(t-1)*Unit:t*Unit,1));                %Total number of new infections per year for xxx years among FSWs, that is incidence from 0:12 months, 12:24 months ....
    IncidenceALL(t)=Incidence1Year(t)+Incidence23YearS(t);                  %Total number of new infections per year among all FSWs and all clients
    
    PrevalenceYear1(t)=100*InfectedFSW(t)/fsw;                              %HIV prevalence among FSWs at 12 months, 24 months...
    PrevalenceYear2(t)=100*InfectedRClients(t)/regm;                        %HIV prevalence among regular clients at 12 months, 24 months...
    PrevalenceYear3(t)=100*InfectedNRClients(t)/nregm;                      %HIV prevalence among non-regular clients at 12 months, 24 months...
    PrevalenceYear23(t)=100*(InfectedRClients(t)+InfectedNRClients(t))/(nregm+regm); %HIV prevalence among clients (regular and non-regular) at 12 months, 24 months...
    PrevalenceYear4(t)=100*InfectedALL(t)/(fsw+regm+nregm);                 %HIV prevalence in the total population (FSWs and clients) at 12 months, 24 months...   
    
    IncidenceR1Year(t)=100*Incidence1Year(t)/Susc(t);                       %Incidence rate among FSWs (injecting and non-injecting)
    IncidenceR2Year(t)=100*Incidence2YearS(t)/SuscRClient(t);               %Incidence rate among regular clients 
    IncidenceR3Year(t)=100*Incidence3YearS(t)/SuscNRClient(t);              %Incidence rate among non-regular clients 
    IncidenceR23Year(t)=100*Incidence23YearS(t)/(SuscRClient(t)+SuscNRClient(t));   %Incidence rate among clients (regular and non-regular)
    IncidenceR4Year(t)=100*Incidence4YearS(t)/(Susc(t)+SuscRClient(t)+SuscNRClient(t)); %Incidence rate among FSWs and clients 

    eps=1.0E-10;
 
    Ndiscr(t)=npartr*(PrevalenceYear2(t)/100)*(1-(1/3*PrevalenceYear2(t)/100));     %number of discordant couples for regular clients assuming that the prevalence among stable partners is 1/3 that in clients
    Ndiscnr(t)=npartnr*(PrevalenceYear3(t)/100)*(1-(1/3*PrevalenceYear3(t)/100));   %number of discordant couples for non-regular clients assuming that the prevalence among stable partners is 1/3 that in clients
    Ndisc(t)=npart*(PrevalenceYear23(t)/100)*(1-(1/3*PrevalenceYear23(t)/100));     %number of discordant couples for clients assuming that the prevalence among stable partners is 1/3 that in clients
    
    Inc_spr(t)=Ndiscr(t)*tpsp;                                                      %Number of incident infections in spouses of regular clients
    Inc_spnr(t)=Ndiscnr(t)*tpsp;                                                    %Number of incident infections in spouses of non-regular clients
    Inc_sp(t)=Ndisc(t)*tpsp;                                                        %Number of incident infections in spouses of clients

    Susc_spr(t)=(1-(1/3*PrevalenceYear2(t)/100))*npartr;                             %Number of susceptible spouses of regular clients in all marital partnerships
    Susc_spnr(t)=(1-(1/3*PrevalenceYear3(t)/100))*npartnr;                          %Number of susceptible spouses of non-regular clients in all marital partnerships
    Susc_sp(t)=(1-(1/3*PrevalenceYear23(t)/100))*npart;                             %Number of susceptible spouses of all clients in all marital partnerships

    IncR_spr(t)=100*Inc_spr(t)/((1-(1/3*PrevalenceYear2(t)/100))*npartr);           %Incidence rate among spouses of regular clients (number of new infections/number of suscpetible spouses in all marital relations)
    IncR_spnr(t)=100*Inc_spnr(t)/((1-(1/3*PrevalenceYear3(t)/100))*npartnr);        %Incidence rate among spouses of non-regular clients (number of new infections/number of suscpetible spouses in all marital relations)
    IncR_sp(t)=100*Inc_sp(t)/((1-(1/3*PrevalenceYear23(t)/100))*npart);             %Incidence rate among spouses of clients (number of new infections/number of suscpetible spouses in all marital relations)

end

Susc_Avg=zeros(YearEND,1);
SuscRClient_Avg=zeros(YearEND,1);
SuscNRClient_Avg=zeros(YearEND,1);
SuscClient_Avg=zeros(YearEND,1);

Incidence1YearS_Avg=zeros(YearEND,1); 
Incidence2YearS_Avg=zeros(YearEND,1); 
Incidence3YearS_Avg=zeros(YearEND,1);
Incidence23YearS_Avg=zeros(YearEND,1);
Incidence4YearS_Avg=zeros(YearEND,1);
 
Incidence1YearIDU_Avg=zeros(YearEND,1);
Incidence1YearnotIDU_Avg=zeros(YearEND,1);
Incidence1Year_Avg=zeros(YearEND,1);
IncidenceALL_Avg=zeros(YearEND,1);

IncidenceR1Year_Avg=zeros(YearEND,1);
IncidenceR2Year_Avg=zeros(YearEND,1);
IncidenceR3Year_Avg=zeros(YearEND,1);
IncidenceR23Year_Avg=zeros(YearEND,1);
IncidenceR4Year_Avg=zeros(YearEND,1);
    
InfectedFSW_Avg=zeros(YearEND,1);                           

InfectedRClients_Avg=zeros(YearEND,1);                       
InfectedNRClients_Avg=zeros(YearEND,1);                 
InfectedClients_Avg=zeros(YearEND,1);
InfectedALL_Avg=zeros(YearEND,1);                           
    
PrevalenceYear1_Avg=zeros(YearEND,1);                         
PrevalenceYear2_Avg=zeros(YearEND,1);
PrevalenceYear3_Avg=zeros(YearEND,1);
PrevalenceYear23_Avg=zeros(YearEND,1);
PrevalenceYear4_Avg=zeros(YearEND,1); 
      
Ndiscr_Avg=zeros(YearEND,1);
Ndiscnr_Avg=zeros(YearEND,1);
Ndisc_Avg=zeros(YearEND,1);

Inc_spr_Avg=zeros(YearEND,1);
Inc_spnr_Avg=zeros(YearEND,1);
Inc_sp_Avg=zeros(YearEND,1);

Susc_spr_Avg=zeros(YearEND,1);
Susc_spnr_Avg=zeros(YearEND,1);
Susc_sp_Avg=zeros(YearEND,1);
    
IncR_spr_Avg=zeros(YearEND,1);
IncR_spnr_Avg=zeros(YearEND,1);
IncR_sp_Avg=zeros(YearEND,1);
    
for t=50:YearEND
    Susc_Avg(t)=sum(Susc(t-49:t))/50;                                       %average number of susceptible FSWs (injecting and non-injecting)
    SuscRClient_Avg(t)=sum(SuscRClient(t-49:t))/50;                         %average number of susceptible regular clients 
    SuscNRClient_Avg(t)=sum(SuscNRClient(t-49:t))/50;                       %average number of susceptible non-regular clients 
    SuscClient_Avg(t)=sum(SuscClient(t-49:t))/50;                           %average number of susceptible clients (regular and non-regular)

    Incidence1YearS_Avg(t)=sum(Incidence1YearS(t-49:t))/50;                 %average number of incident infections per year among FSWs due to sexual transmission
    Incidence2YearS_Avg(t)=sum(Incidence2YearS(t-49:t))/50;                 %average number of incident infections per year among regular clients due to sexual transmission
    Incidence3YearS_Avg(t)=sum(Incidence3YearS(t-49:t))/50;                 %average number of incident infections per year among non-regular clients due to sexual transmission
    Incidence23YearS_Avg(t)=sum(Incidence23YearS(t-49:t))/50;               %average number of incident infections per year among clients (regular or non-regular) due to sexual transmission
    Incidence4YearS_Avg(t)=sum(Incidence4YearS(t-49:t))/50;                 %average number of incident infections per year among FSWs and clients due to sexual transmission
    
    Incidence1Year_Avg(t)=sum(Incidence1Year(t-49:t))/50;                   %average number of new infections per year among FSWs (injecting and non-injecting)
    IncidenceALL_Avg(t)=sum(Incidence1Year(t-49:t))/50;                     %average number of new infections per year in all FSWs and all clients

    IncidenceR1Year_Avg(t)=100*Incidence1Year_Avg(t)/Susc_Avg(t);               %average incidence rate among FSWs (injecting and non-injecting)
    IncidenceR2Year_Avg(t)=100*Incidence2YearS_Avg(t)/SuscRClient_Avg(t);       %average incidence rate among regular clients
    IncidenceR3Year_Avg(t)=100*Incidence3YearS_Avg(t)/SuscNRClient_Avg(t);      %average incidence rate among non-regular clients
    IncidenceR23Year_Avg(t)=100*Incidence23YearS_Avg(t)/SuscClient_Avg(t);      %average incidence rate among clients (regular and non-regular)
    IncidenceR4Year_Avg(t)=100*Incidence4YearS_Avg(t)/(Susc_Avg(t)+SuscClient_Avg(t)); %average incidence rate due to sexual transmission among FSWs and clients

    InfectedFSW_Avg(t)=sum(InfectedFSW(t-49:t))/50;                             %average number of infected FSWs (injecting and non-injecting) per year
 
    InfectedRClients_Avg(t)=sum(InfectedRClients(t-49:t))/50;                       %average number of infected regular clients 
    InfectedNRClients_Avg(t)=sum(InfectedNRClients(t-49:t))/50;                     %average number of infected non-regular clients 
    InfectedClients_Avg(t)=(InfectedRClients_Avg(t)+InfectedNRClients_Avg(t));    %average number of infected clients 
    InfectedALL_Avg(t)=sum(InfectedALL(t-49:t))/50;                             %average number of infected FSWs (all) and clients (all) 
    
    PrevalenceYear1_Avg(t)=100*InfectedFSW_Avg(t)/fsw;                          %average HIV prevalence among FSWs (injecting and non-injecting)
    PrevalenceYear2_Avg(t)=100*InfectedRClients_Avg(t)/regm;                    %average HIV prevalence among regular clients
    PrevalenceYear3_Avg(t)=100*InfectedNRClients_Avg(t)/nregm;                  %average HIV prevalence among non-regular clients
    PrevalenceYear23_Avg(t)=100*(InfectedRClients_Avg(t)+InfectedNRClients_Avg(t))/(nregm+regm); %average HIV prevalence among clients (regular and non-regular) 
    PrevalenceYear4_Avg(t)=100*InfectedALL_Avg(t)/(fsw+regm+nregm);             %average HIV prevalence in the total population (FSWs and clients) 
    

    Ndiscr_Avg(t)=sum(Ndiscr(t-49:t))/50;                                       %average number of discordant partnerships among regular clients
    Ndiscnr_Avg(t)=sum(Ndiscnr(t-49:t))/50;                                     %average number of discordant partnerships among non-regular clients
    Ndisc_Avg(t)=sum(Ndisc(t-49:t))/50;                                         %average number of discordant partnerships among clients (regular and non-regular)

    Inc_spr_Avg(t)=Ndiscr_Avg(t)*tpsp;                                          %average number of incident infections among spouses of regular clients
    Inc_spnr_Avg(t)=Ndiscnr_Avg(t)*tpsp;                                        %average number of incident infections among spouses of non-regular clients
    Inc_sp_Avg(t)=Ndisc_Avg(t)*tpsp;                                            %average number of incident infections among spouses of all clients
     
    Susc_spr_Avg(t)=(1-(1/3*PrevalenceYear2_Avg(t)/100))*npartr;                %average number of susceptible spouses of regular clients in all marital partnerships
    Susc_spnr_Avg(t)=(1-(1/3*PrevalenceYear3_Avg(t)/100))*npartnr;              %average number of susceptible spouses of non-regular clients in all marital partnerships
    Susc_sp_Avg(t)=(1-(1/3*PrevalenceYear23_Avg(t)/100))*npart;                 %average number of susceptible spouses of all clients in all marital partnerships

    IncR_spr_Avg(t)=100*Inc_spr_Avg(t)/Susc_spr_Avg(t);                         %average incident rate among spouses of regular clients
    IncR_spnr_Avg(t)=100*Inc_spnr_Avg(t)/Susc_spnr_Avg(t);                      %average incident rate among spouses of non-regular clients
    IncR_sp_Avg(t)=100*Inc_sp_Avg(t)/Susc_sp_Avg(t);                            %average incident rate among spouses of all clients
end


