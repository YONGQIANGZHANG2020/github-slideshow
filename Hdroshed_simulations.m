%% YQZ, 07 Aug 2020

clc;clear;close all;

%% load hydroshed catchment attributes
load('E:\GlobalData\hydroshed_attribute\Hydroshed_Inputs','Table_Inputs','Table_Final');
load (['RTE_TrainingCrossVal'], 'Mdl_Matrix');


%% Prediction variable order
VariableNames ={'catchsize','mean_elev','mean_slope','permeability','forest_ratio',...
           'clay','gravel','sand','silt','mean_Tmean','mean_Tmax','mean_Tmin','mean_P','mean_PET','mean_LAI',...
           'std_Tmean','std_Tmax','std_Tmin','std_P','std_PET','std_LAI',...
          'seasonality_Tmean','seasonality_Tmax','seasonality_Tmin','seasonality_P','seasonality_PET','seasonality_LAI'};      

Elasticities_names = {'Monthly P', 'Season_P', 'Yearly P', ...
    'Monthly PET' ,'Season PET', 'Yearly PET', ...
    'Monthly LAI', 'Season LAI', 'Yearly LAI', ...
    'Monthly S' ,'Season S', 'Yearly S'};
                 
for ii = 1:size(Mdl_Matrix,1) 
  fprintf('ii = %.0d\n',ii);  
 Elasticities(:,ii)      = predict(Mdl_Matrix{ii},[table2array(Table_Final)]);
%  Elasticities(:,ii)      = predict(Mdl_Matrix{ii},[table2array(Table_Final) zeros(size(Table_Final,1),3)]);
end

save('Hydroshed_predictions','Elasticities','Elasticities_names');

