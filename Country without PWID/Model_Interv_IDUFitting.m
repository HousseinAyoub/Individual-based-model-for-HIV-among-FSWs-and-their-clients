
function [PrevalenceYear1_Avg]=Model_Interv_IDUFitting(DT,Unit,fsw,regm,nregm,tpHIVreg,tpHIVn,transHIV,removelink_reg,removelink_nreg,YearEND,YearBURNIN,BURN_IN,TEND,INIT_INF_fsw,INIT_INF_regm,INIT_INF_nregm,CoverageMC,CoverageCond,CovARTFSW,CovARTClient,CovPrEPFSW,CovPrEPClient,eART,eMC,eCond,ePrEP,acq_regm,acq_nregm,mu,muHIV)

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
%     if mod(t,12)==0  
%         disp(t);
%     end
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
        
%%%%%%%%%%%%%%%%%%% HIV transmission due IDUs
%%%%%%%%%%%%%%%%%%%%%%%%HIV among PWID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%   for cnt= 1:fsw
%       if (IDU_fsw(cnt)==1) && (state_fsw(cnt)==0) 
%           if rand<HIVPWID
%              state_fsw(cnt)=1;
%              incidence_fsw_IDU = incidence_fsw_IDU + 1;
%           end
%       end
%   end


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
                            if rand<CovCond % if he had condom
                                if rand<CovPrEPM % if he is on PrEP (ART+ & MC+ & Cond+ & PrEP+)
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
                                if rand<CovPrEPM
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
                         if rand<CovCond % if non-regular client uses condoms
                             if rand<CovPrEPM % if he is on PrEP (ART+ & MC+ & Cond+ & PrEP+)
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
                             if rand<CovPrEPM % if he is on PrEP (ART+ & MC- & Cond+ & PrEP+)
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
                         if rand<CovCond % if he is using comdoms
                             if rand<CovPrEPM % if he is on PrEP (ART- & MC+ & Cond+ & PrEP+)
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
                           if rand<CovPrEPM % if he is on PrEP (ART- & MC- & Cond+ & PrEP+)
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
                         if rand<CovCond % if he uses condoms
                             if rand<CovPrEPM % (ART- & MC+ & Cond+ & PrEP+)
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
                           if rand<CovPrEPM % (ART- & MC- & Cond+ & PrEP+)
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
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW is on PrEP (ART+ & Cond+ &PrEP+)
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
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW is on PrEP (ART- & Cond+ &PrEP+)
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
                      for cnt2 = 1:fsw % for FSWs
                          if rand<CovPrEPF %if the FSW was on PrEP (ART+ & Cond+ & PreP+)
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
                       for cnt2 = 1:fsw % for FSWs
                           if rand<CovPrEPF %if the FSW was on PrEP (ART- & Cond+ & PreP+)
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


    for cnt = 1:fsw
           if state_fsw(cnt) ~= 0
              prevalence_fsw = prevalence_fsw + 1;
           else
              susceptible_fsw= susceptible_fsw + 1;
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

end

for t=1:YearEND 
    InfectedFSW(t)=Number_I_FSW((t*Unit),1);                                %Number of infected FSWs (injecting and non-injecting) at 12 months, 24 months...
end
    
InfectedFSW_Avg=zeros(YearEND,1);                           
PrevalenceYear1_Avg=zeros(YearEND,1); 

for t=50:YearEND
    InfectedFSW_Avg(t)=sum(InfectedFSW(t-49:t))/50;                             %average number of infected FSWs (injecting and non-injecting) per year
    PrevalenceYear1_Avg(t)=100*InfectedFSW_Avg(t)/fsw;                          %average HIV prevalence among FSWs (injecting and non-injecting)
end


