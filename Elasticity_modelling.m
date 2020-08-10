%% Last update by YQZ, 10/08/2020



%% data source
clear;close all; clc

load ('Inputs','Monthly_matrix','Yearly_matrix','GloCat','CatInfo','Headers','Var_finals');

% Values for all drivers must be not less than zero
Monthly_matrix(Monthly_matrix<0) = nan;
Yearly_matrix(Yearly_matrix<0)   = nan;

%% Standardise (Xi - mean (X)) / mean (X), yearly, monthly, and seasonal
for ii=1:size(Yearly_matrix,3)
  % yearly
       DATA                        = Yearly_matrix(:,:,ii); 
       DATA(isnan(DATA(:,1))==1,:) = [];
       Yearly_matrix_stand(:,:,ii) = (Yearly_matrix(:,:,ii) - ones(size(Yearly_matrix,1),1)*nanmean(DATA,1))./...
                                     (ones(size(Yearly_matrix,1),1)*nanmean(DATA,1));
       TEMP                        = (Yearly_matrix(:,:,ii) - ones(size(Yearly_matrix,1),1)*nanmean(DATA,1));
       Yearly_matrix_delta(:,ii)   = TEMP(:,1);
  % monthly
       DATA                        = Monthly_matrix(:,:,ii);
       DATA(isnan(DATA(:,1))==1,:) = [];
       Monthly_matrix_stand(:,:,ii) = (Monthly_matrix(:,:,ii) - ones(size(Monthly_matrix,1),1)*nanmean(DATA,1))./...
                                    (ones(size(Monthly_matrix,1),1)*nanmean(DATA,1));
       TEMP                         = (Monthly_matrix(:,:,ii) - ones(size(Monthly_matrix,1),1)*nanmean(DATA,1));
       Monthly_matrix_delta(:,ii)   = TEMP(:,1);

  % seasonal
    DATA  = Monthly_matrix(:,:,ii);
    for jj=1:size(DATA,2)
        Temp = reshape(DATA(:,jj),3,size(DATA,1)/3);
        if jj>1 & jj<5    
        Seasons(:,jj) = Temp(1,:)'; % for Qb, to select the start of the month
        else
        Seasons(:,jj) = mean(Temp,1)'; % for other doing block average at each 3-months block
        end
    end
    Season_matrix(:,:,ii)       = Seasons ;
    Season_matrix_stand(:,:,ii) = (Seasons - ones(size(Seasons,1),1)*nanmean(Seasons,1))./ ...
                                  (ones(size(Seasons,1),1)*nanmean(Seasons,1));
    TEMP                        = (Seasons - ones(size(Seasons,1),1)*nanmean(Seasons,1));
    Season_matrix_delta(:,ii)   = TEMP(:,1);
    
    clear Seasons;
end


%% Criteria and parameters
% 1. CatAreaThresh = 10000; % km2 
% 2. YearThresh    = 10;
% 3. Irrigation    < 10%
% 4. Water body    < 10%
% 5. Reservior cap < 10%

% P               - 2
% P, ETp          - 2 + 8
% P, ETp, LAI     - 10 + 16
% P, ETp, S       - 26 + 24
% P, ETp, LAI, S  - 50 + 48 = 98

P_column   = [5 8];
ETp_column = [6 7 9 10];
S_column   = [2:4];
LAI_column = [11 12];

Var_finals ={ 'Q (mm/d)','Qb_LH (mm/d)','Qb_CM (mm/d)','Qb_Eck (mm/d)','WFDEI_P (mm/d)' ...
           'WFDEI_PET_Penman (mm/d)','WFDEI_PET_PT (mm/d)','Princeton_P (mm/d)', ...
           'Princeton_PET_Penman (mm/d)','Princeton_PET_PT (mm/d)','Globalmap_LAI','GIMMS_LAI'};

%% load modelling

 % 4 variables
start = 0;   
for jj=1:size(P_column,2)
     for kk=1:size(ETp_column,2)
         for ll=1:size(LAI_column,2)
             for mm=1:size(S_column,2)
             start = start +1;
             X_4var {start,1} =  [P_column(jj) ETp_column(kk) LAI_column(ll) S_column(mm)];
         end
         end 
     end
end
 
X = X_4var; 
%% Monthly modelling 

for ii = 1:size(Yearly_matrix,3)
    
     DATA                   = Monthly_matrix_stand(:,:,ii);
     DATA_original          = Monthly_matrix(:,:,ii);
     DeltaY                 = Monthly_matrix_delta(:,ii);
          
   for jj = 1:length(X)
            Columns                = X{jj};
            INPUTS                 = [DATA(:,Columns) ones(size(DATA,1),1)];
            [B,BINT,R,RINT,STATS]  = regress(DATA(:,1), INPUTS);
            if mean(B)~=0 % making sure the MLR is valid
            Monthly_outputs.Parameters_matrix{jj,ii} = B';
            Monthly_outputs.R2_matrix(jj,ii)         = STATS(1); % R^2
            else
            Monthly_outputs.Parameters_matrix{jj,ii} = B' + nan;
            Monthly_outputs.R2_matrix(jj,ii)         = STATS(1) +nan; % R^2
            end
            
        Nums = find(isnan(INPUTS(:,1))==0 & isnan(INPUTS(:,2))==0 & ...
            isnan(INPUTS(:,3))==0 & isnan(INPUTS(:,4))==0);   
        Monthly_outputs.Importance(ii,:,jj) = [sum(abs(INPUTS(Nums,1)*abs(B(1)))) ...
                                      sum(abs(INPUTS(Nums,2)*abs(B(2)))) ...
                                      sum(abs(INPUTS(Nums,3)*abs(B(3)))) ...
                                      sum(abs(INPUTS(Nums,4)*abs(B(4))))]; % importance, P, ETp, LAI, S
%         Monthly_outputs.MetaData{ii,jj-FourVar_Threshold}  = single([DeltaY DATA(:,1), INPUTS]); % matadata


   end
   
    Nums                                        = find(isnan(DATA_original(:,1))==1);
    Monthly_outputs.P(ii,1:length(P_column))    = nanmean(DATA_original(Nums,P_column));
    Monthly_outputs.Q(ii,1)                     = nanmean(DATA_original(Nums,1));
    Monthly_outputs.S(ii,1:length(S_column))    = nanmean(DATA_original(Nums,S_column));
    Monthly_outputs.LAI(ii,1:length(LAI_column))= nanmean(DATA_original(Nums,LAI_column));
    Monthly_outputs.PET(ii,1:length(ETp_column))= nanmean(DATA_original(Nums,ETp_column));
 
    Monthly_outputs.AI(ii,1:length(P_column))   = (nanmean(DATA_original(:,ETp_column([1,3])))./nanmean(DATA_original(:,P_column )) + ...
                                                nanmean(DATA_original(:,ETp_column([2,4])))./nanmean(DATA_original(:,P_column )))/2 ;

end

%% Seasonal modelling 
for ii = 1:size(Yearly_matrix,3)
    
     DATA                   = Season_matrix_stand(:,:,ii);
     DATA_original          = Season_matrix(:,:,ii);
     DeltaY                 = Season_matrix_delta(:,ii);
          
   for jj = 1:length(X)
            Columns                = X{jj};
            INPUTS                 = [DATA(:,Columns) ones(size(DATA,1),1)];
            [B,BINT,R,RINT,STATS]  = regress(DATA(:,1), INPUTS);
            
            if mean(B)~=0 % making sure the MLR is valid
            Season_outputs.Parameters_matrix{jj,ii} = B';
            Season_outputs.R2_matrix(jj,ii)         = STATS(1); % R^2
            else
            Season_outputs.Parameters_matrix{jj,ii} = B' + nan;
            Season_outputs.R2_matrix(jj,ii)         = STATS(1) +nan; % R^2
            end
            
        Nums = find(isnan(INPUTS(:,1))==0 & isnan(INPUTS(:,2))==0 & ...
            isnan(INPUTS(:,3))==0 & isnan(INPUTS(:,4))==0);   
        Season_outputs.Importance(ii,:,jj) = [sum(abs(INPUTS(Nums,1)*abs(B(1)))) ...
                                      sum(abs(INPUTS(Nums,2)*abs(B(2)))) ...
                                      sum(abs(INPUTS(Nums,3)*abs(B(3)))) ...
                                      sum(abs(INPUTS(Nums,4)*abs(B(4))))]; % importance, P, ETp, LAI, S

   end
    Nums                                       = find(isnan(DATA_original(:,1))==1);
    Season_outputs.P(ii,1:length(P_column))    = nanmean(DATA_original(Nums,P_column));
    Season_outputs.Q(ii,1)                     = nanmean(DATA_original(Nums,1));
    Season_outputs.S(ii,1:length(S_column))    = nanmean(DATA_original(Nums,S_column));
    Season_outputs.LAI(ii,1:length(LAI_column))= nanmean(DATA_original(Nums,LAI_column));
    Season_outputs.PET(ii,1:length(ETp_column))= nanmean(DATA_original(Nums,ETp_column));
 
    Season_outputs.AI(ii,1:length(P_column)) = (nanmean(DATA_original(Nums,ETp_column([1,3])))./nanmean(DATA_original(Nums,P_column )) + ...
                                                nanmean(DATA_original(Nums,ETp_column([2,4])))./nanmean(DATA_original(Nums,P_column )))/2;

end
%% Yearly modelling 
for ii = 1:size(Yearly_matrix,3)
    
     DATA                   = Yearly_matrix_stand(:,:,ii);
     DATA_original          = Yearly_matrix(:,:,ii);
     DeltaY                 = Yearly_matrix_delta(:,ii);
           
   for jj = 1:length(X)
            Columns                = X{jj};
            INPUTS                 = [DATA(:,Columns) ones(size(DATA,1),1)];
            [B,BINT,R,RINT,STATS]  = regress(DATA(:,1), INPUTS);
          
%          for kk=1:size(INPUTS,2)  
%          end
            if mean(B)~=0 
            Yearly_outputs.Parameters_matrix{jj,ii} = B';
            Yearly_outputs.R2_matrix(jj,ii)         = STATS(1); % R^2
            else
            Yearly_outputs.Parameters_matrix{jj,ii} = B' + nan;
            Yearly_outputs.R2_matrix(jj,ii)         = STATS(1) +nan; % R^2
            end
            
        Nums = find(isnan(INPUTS(:,1))==0 & isnan(INPUTS(:,2))==0 & ...
            isnan(INPUTS(:,3))==0 & isnan(INPUTS(:,4))==0);   
        Yearly_outputs.Importance(ii,:,jj)    = [sum(abs(INPUTS(Nums,1)*abs(B(1)))) ...
                                      sum(abs(INPUTS(Nums,2)*abs(B(2)))) ...
                                      sum(abs(INPUTS(Nums,3)*abs(B(3)))) ...
                                      sum(abs(INPUTS(Nums,4)*abs(B(4))))]; % importance, P, ETp, LAI, S

   end
   
    Nums                                       = find(isnan(DATA_original(:,1))==1);
    Yearly_outputs.P(ii,1:length(P_column))    = nanmean(DATA_original(Nums,P_column));
    Yearly_outputs.Q(ii,1)                     = nanmean(DATA_original(Nums,1));
    Yearly_outputs.S(ii,1:length(S_column))    = nanmean(DATA_original(Nums,S_column));
    Yearly_outputs.LAI(ii,1:length(LAI_column)) = nanmean(DATA_original(Nums,LAI_column));
    Yearly_outputs.PET(ii,1:length(ETp_column)) = nanmean(DATA_original(Nums,ETp_column));

    Yearly_outputs.AI(ii,1:length(P_column)) = (nanmean(DATA_original(Nums,ETp_column([1,3])))./nanmean(DATA_original(Nums,P_column )) + ...
                                                nanmean(DATA_original(Nums,ETp_column([2,4])))./nanmean(DATA_original(Nums,P_column )))/2;

end

   
save ('Sensitivity_MLR_Month_Season_Year','GloCat','CatInfo','Headers','Season_outputs','Monthly_outputs','Yearly_outputs','-v7.3');


