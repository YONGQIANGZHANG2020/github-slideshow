%% YQZ, 06 August 2020
% Using new catchment attributes, new permeability dataset:
% Gleeson, T., Smith. L., Jansen, N., Hartmann, J., D¨¹rr, H., Manning, A.H.,van Beek, R. and A.M. Jellinek (2011)


clear;clc;close all;

%% Load Elasticities and attributes

load ('Sensitivity_MLR_Month_Season_Year','GloCat','CatInfo','Headers','Season_outputs','Monthly_outputs','Yearly_outputs');


load ('CatchAttributes_postpr_newPerm.mat', 'CatchAttributes_Table','IDs');
CatchAttributes_YH = CatchAttributes_Table(IDs,:);
CatchAttributes_YH.Forests(CatchAttributes_YH.Forests>1) = 1;
Soil_Attributes = [nanmean(table2array(CatchAttributes_YH(:,[32:38])),2) ...
                   nanmean(table2array(CatchAttributes_YH(:,[39:45])),2) ...
                   nanmean(table2array(CatchAttributes_YH(:,[46:52])),2) ...
                   nanmean(table2array(CatchAttributes_YH(:,[53:59])),2)];

load ('CatchAttributes_Hylke.mat', 'CatchAttributes_Table','IDs');
CatchAttributes_Table_Hylke = CatchAttributes_Table(IDs,:);


TEMP = [CatchAttributes_Table_Hylke.catchsize CatchAttributes_Table_Hylke.mean_elev CatchAttributes_Table_Hylke.mean_slope ...
        CatchAttributes_YH.PermeabilityNew CatchAttributes_YH.Forests Soil_Attributes ...
        CatchAttributes_Table_Hylke.Mean_Tmean CatchAttributes_Table_Hylke.Mean_Tmax CatchAttributes_Table_Hylke.Mean_Tmin...
        nanmean([CatchAttributes_Table_Hylke.Mean_WFDEI_P CatchAttributes_Table_Hylke.Mean_Princeton_P],2) ...
        nanmean([CatchAttributes_Table_Hylke.Mean_WFDEI_PET_Penman CatchAttributes_Table_Hylke.Mean_WFDEI_PET_PT CatchAttributes_Table_Hylke.Mean_Princeton_PET_Penman CatchAttributes_Table_Hylke.Mean_Princeton_PET_PT],2) ...
        nanmean([CatchAttributes_Table_Hylke.Mean_GIMMS_LAI CatchAttributes_Table_Hylke.Mean_Globalmap_LAI],2) ...
        CatchAttributes_Table_Hylke.Std_Tmean CatchAttributes_Table_Hylke.Std_Tmax CatchAttributes_Table_Hylke.Std_Tmin...
        nanmean([CatchAttributes_Table_Hylke.Std_WFDEI_P CatchAttributes_Table_Hylke.Std_Princeton_P],2) ...
        nanmean([CatchAttributes_Table_Hylke.Std_WFDEI_PET_Penman CatchAttributes_Table_Hylke.Std_WFDEI_PET_PT CatchAttributes_Table_Hylke.Std_Princeton_PET_Penman CatchAttributes_Table_Hylke.Std_Princeton_PET_PT],2) ...
        nanmean([CatchAttributes_Table_Hylke.Std_GIMMS_LAI CatchAttributes_Table_Hylke.Std_Globalmap_LAI],2) ...
        CatchAttributes_Table_Hylke.Season_Tmean CatchAttributes_Table_Hylke.Season_Tmax CatchAttributes_Table_Hylke.Season_Tmin...
        nanmean([CatchAttributes_Table_Hylke.Season_WFDEI_P CatchAttributes_Table_Hylke.Season_Princeton_P],2) ...
        nanmean([CatchAttributes_Table_Hylke.Season_WFDEI_PET_Penman CatchAttributes_Table_Hylke.Season_WFDEI_PET_PT CatchAttributes_Table_Hylke.Season_Princeton_PET_Penman CatchAttributes_Table_Hylke.Season_Princeton_PET_PT],2) ...
        nanmean([CatchAttributes_Table_Hylke.Season_GIMMS_LAI CatchAttributes_Table_Hylke.Season_Globalmap_LAI],2) ...
        ];

Table_Final = array2table(TEMP,...
    'VariableNames',{'catchsize','mean_elev','mean_slope','permeability','forest_ratio',...
                      'clay','gravel','sand','silt','mean_Tmean','mean_Tmax','mean_Tmin','mean_P','mean_PET','mean_LAI',...
                     'std_Tmean','std_Tmax','std_Tmin','std_P','std_PET','std_LAI',...
                      'seasonality_Tmean','seasonality_Tmax','seasonality_Tmin','seasonality_P','seasonality_PET','seasonality_LAI'});        

Table_Final.std_PET(isnan(Table_Final.mean_PET)==1)         = nan;
Table_Final.seasonality_PET(isnan(Table_Final.mean_PET)==1) = nan;

%% Prediction variables:

% monthly - ensemble median of 48 replicates
Parameters_matrix = Monthly_outputs.Parameters_matrix;
Parameters_matrix = cell2mat (Parameters_matrix);
Parameters_matrix = nanmedian(Parameters_matrix,1);

Monthly_P_Elas    = Parameters_matrix(:,[1:5:size(Parameters_matrix,2)]); 
Monthly_P_Elas    = reshape(Monthly_P_Elas,size(Monthly_P_Elas,1).*size(Monthly_P_Elas,2),1);
Monthly_PET_Elas    = Parameters_matrix(:,[2:5:size(Parameters_matrix,2)]); 
Monthly_PET_Elas    = reshape(Monthly_PET_Elas,size(Monthly_PET_Elas,1).*size(Monthly_PET_Elas,2),1);
Monthly_LAI_Elas    = Parameters_matrix(:,[3:5:size(Parameters_matrix,2)]); 
Monthly_LAI_Elas    = reshape(Monthly_LAI_Elas,size(Monthly_LAI_Elas,1).*size(Monthly_LAI_Elas,2),1);
Monthly_S_Elas    = Parameters_matrix(:,[4:5:size(Parameters_matrix,2)]); 
Monthly_S_Elas    = reshape(Monthly_S_Elas,size(Monthly_S_Elas,1).*size(Monthly_S_Elas,2),1);

% seasonal - ensemble median of 48 replicates
Parameters_matrix = Season_outputs.Parameters_matrix;
Parameters_matrix = cell2mat (Parameters_matrix);
Parameters_matrix = nanmedian(Parameters_matrix,1);
Season_P_Elas    = Parameters_matrix(:,[1:5:size(Parameters_matrix,2)]); 
Season_P_Elas    = reshape(Season_P_Elas,size(Season_P_Elas,1).*size(Season_P_Elas,2),1);
Season_PET_Elas    = Parameters_matrix(:,[2:5:size(Parameters_matrix,2)]); 
Season_PET_Elas    = reshape(Season_PET_Elas,size(Season_PET_Elas,1).*size(Season_PET_Elas,2),1);
Season_LAI_Elas    = Parameters_matrix(:,[3:5:size(Parameters_matrix,2)]); 
Season_LAI_Elas    = reshape(Season_LAI_Elas,size(Season_LAI_Elas,1).*size(Season_LAI_Elas,2),1);
Season_S_Elas    = Parameters_matrix(:,[4:5:size(Parameters_matrix,2)]); 
Season_S_Elas    = reshape(Season_S_Elas,size(Season_S_Elas,1).*size(Season_S_Elas,2),1);

% yearly - ensemble median of 48 replicates
Parameters_matrix = Yearly_outputs.Parameters_matrix;
Parameters_matrix = cell2mat (Parameters_matrix);
Parameters_matrix = nanmedian(Parameters_matrix,1);
Yearly_P_Elas    = Parameters_matrix(:,[1:5:size(Parameters_matrix,2)]); 
Yearly_P_Elas    = reshape(Yearly_P_Elas,size(Yearly_P_Elas,1).*size(Yearly_P_Elas,2),1);
Yearly_PET_Elas    = Parameters_matrix(:,[2:5:size(Parameters_matrix,2)]); 
Yearly_PET_Elas    = reshape(Yearly_PET_Elas,size(Yearly_PET_Elas,1).*size(Yearly_PET_Elas,2),1);
Yearly_LAI_Elas    = Parameters_matrix(:,[3:5:size(Parameters_matrix,2)]); 
Yearly_LAI_Elas    = reshape(Yearly_LAI_Elas,size(Yearly_LAI_Elas,1).*size(Yearly_LAI_Elas,2),1);
Yearly_S_Elas      = Parameters_matrix(:,[4:5:size(Parameters_matrix,2)]); 
Yearly_S_Elas      = reshape(Yearly_S_Elas,size(Yearly_S_Elas,1).*size(Yearly_S_Elas,2),1);

Elasticities = [Monthly_P_Elas Season_P_Elas Yearly_P_Elas ...
    Monthly_PET_Elas Season_PET_Elas Yearly_PET_Elas ...
    Monthly_LAI_Elas Season_LAI_Elas Yearly_LAI_Elas ...
    Monthly_S_Elas Season_S_Elas Yearly_S_Elas];
%% remove outliers
% if the elasticity is larger than 10, taken as the outliers
Elasticities (abs(Elasticities)>10) = nan;

%% putting together

Elasticities_names = {'Monthly_P_Elas', 'Season_P_Elas', 'Yearly_P_Elas', ...
    'Monthly_PET_Elas' ,'Season_PET_Elas', 'Yearly_PET_Elas', ...
    'Monthly_LAI_Elas', 'Season_LAI_Elas', 'Yearly_LAI_Elas', ...
    'Monthly_S_Elas' ,'Season_S_Elas', 'Yearly_S_Elas'};

plot_names = {'monthly P', 'seasonal P', 'annual P', ...
    'monthly PET' ,'seasonal PET', 'annual PET', ...
    'monthly LAI', 'seasonal LAI', 'annual LAI', ...
    'monthly S' ,'season S', 'annual S'};

for ii=1:size(Elasticities,2)
    
       YX                   = [Elasticities(:,ii) table2array(Table_Final)]; 
       [ResponseV,Mdl]      = RandomForest(YX);
       imp                  = oobPermutedPredictorImportance(Mdl);
       [a,b]                = sort(imp,'descend','MissingPlacement','last');
       CatchAttributes_sort = b;
        Mdl_Matrix{ii,1}    = Mdl;
        Imp_Matrix{ii,1}    = imp;
        CatchAttributes_sort_Matrix{ii,1} = CatchAttributes_sort;
        ResponseV_Matrix{:,:,ii} = ResponseV;

      for jj = 1:size(ResponseV,2)-1
        NSE_HylkeFinal(ii,jj)   = CE(ResponseV(:,[1 jj+1]));
        KGE_HylkeFinal(ii,jj)   = klinggupta(ResponseV(:,jj+1), ResponseV(:,1));
        TEMP                    = corrcoef(ResponseV(:,[1 jj+1]));
        r2_HylkeFinal(ii,jj)    = TEMP(2,1).^2;
      end
      
      fprintf('ii = %.0d, r2 = %.2f\n', ii,TEMP(2,1).^2);



end
save(['RTE_TrainingCrossVal'], 'NSE_HylkeFinal','KGE_HylkeFinal','r2_HylkeFinal','ResponseV_Matrix','Mdl_Matrix','Imp_Matrix','CatchAttributes_sort_Matrix','Elasticities_names','Table_Final','-v7.3');

