% Import GloRiSe NEW FOR ALL ME
load SedimentDatabase

% creat mean function omitting empty entries
meanOmitNan = @(x) mean(x,'omitNan');
medOmitNan = @(x) median(x,'omitNan');

% Merge ME and Min data with date
% Mineral and petrographic data is included to provide the merging strategy by sample_ID, but not used here
ME.Sample_ID = categorical(ME.Sample_ID);
Mineral.Sample_ID = categorical(Mineral.Sample_ID);
ID.Sample_ID = categorical(ID.Sample_ID);

%ME
[commonBasinID,~] = intersect(ME.Sample_ID,ID.Sample_ID);
for i=1:size(commonBasinID,1)
    
    ME_idx=find(ME.Sample_ID==commonBasinID(i)); %find position of similar entries in RiverNames table
    ID_idx=find(ID.Sample_ID==commonBasinID(i)); %find position of similar entries in Basin-wise median table
    
    if size(ID_idx,1) ~=1 %this if-loop is here because there are rows in the ID table that have the same ID#
        ID_idx=ID_idx(1,1);
    end
    
    if ~isempty(ME_idx) && ~isempty(ID_idx)
    ME(ME_idx,52) = ID(ID_idx,4);
    ME(ME_idx,53) = ID(ID_idx,5);
    ME(ME_idx,54) = ID(ID_idx,6);
    ME(ME_idx,55) = ID(ID_idx,3);
    end  
    
end
%rename added variable
ME.Properties.VariableNames{52} = 'd';
ME.Properties.VariableNames{53} = 'm';
ME.Properties.VariableNames{54} = 'y';
ME.Properties.VariableNames{55} = 'datetime';


%Mineral
[commonBasinID2,~] = intersect(Mineral.Sample_ID,ID.Sample_ID);
for i=1:size(commonBasinID2,1)
    
    Mineral_idx=find(Mineral.Sample_ID==commonBasinID2(i)); 
    ID_idx=find(ID.Sample_ID==commonBasinID2(i)); 
    
    if size(ID_idx,1) ~=1 
        ID_idx=ID_idx(1,1);
    end
    
    if ~isempty(Mineral_idx) && ~isempty(ID_idx)
    Mineral(Mineral_idx,80) = ID(ID_idx,4);
    Mineral(Mineral_idx,81) = ID(ID_idx,5);
    Mineral(Mineral_idx,82) = ID(ID_idx,6);
    Mineral(Mineral_idx,83) = ID(ID_idx,3);
    end  
    
end
%rename added variable
Mineral.Properties.VariableNames{80} = 'D';
Mineral.Properties.VariableNames{81} = 'M';
Mineral.Properties.VariableNames{82} = 'Y';
Mineral.Properties.VariableNames{83} = 'Datetime';

%merge me & Min
Mineral.Properties.VariableNames{1} = 'sample_ID';
Mineral.Properties.VariableNames{2} = 'location_ID';
Mineral.Properties.VariableNames{3} = 'seaCat';
Mineral.Properties.VariableNames{4} = 'observationtype';
Mineral.Properties.VariableNames{5} = 'Sampletype';
Mineral.Properties.VariableNames{6} = 'basin_ID';
Mineral.Properties.VariableNames{7} = 'original_UnitMin';
Mineral.Properties.VariableNames{9} = 'rep_ID';
Mineral.Properties.VariableNames{10} = 'Filtersizem';
Mineral.Properties.VariableNames{11} = 'Sievesizem';
Mineral.Properties.VariableNames{12} = 'dischargem3_s';
Mineral.Properties.VariableNames{13} = 'TSS_mgL';
Mineral.Properties.VariableNames{14} = 'Sand';
Mineral.Properties.VariableNames{15} = 'Silt';
Mineral.Properties.VariableNames{16} = 'Clay';
Mineral.Properties.VariableNames{17} = 'Avggrainsizem';
ME.Sample_ID = string(ME.Sample_ID);
Mineral.sample_ID = string(Mineral.sample_ID);

%preallocate table spaces
Omin = nan(size(ME,1),size(Mineral,2));
Omin = array2table(Omin);
Omin = splitvars(Omin);
Omin.Properties.VariableNames = Mineral.Properties.VariableNames;
MeMin = [ME Omin];

Ome = nan(size(Mineral,1),size(ME,2));
Ome = array2table(Ome);
Ome = splitvars(Ome);
Ome.Properties.VariableNames = ME.Properties.VariableNames;
Min2 = [Mineral Ome];

[commonSample,~] = intersect(MeMin.Sample_ID,Min2.sample_ID);
for i=1:size(commonSample,1)
    
    MeMin_idx=find(MeMin.Sample_ID == commonSample(i)); 
    Min2_idx=find(Min2.sample_ID == commonSample(i)); 
    
    if size(MeMin_idx,1) ~=1 %this if-loop is here because there are rows in the ID table that have the same ID#
        MeMin_idx=MeMin_idx(1,1);
    end
    
    if size(Min2_idx,1) ~=1 %this if-loop is here because there are rows in the ID table that have the same ID#
        Min2_idx=Min2_idx(1,1);
    end
    
    
    if ~isempty(MeMin_idx) && ~isempty(Min2_idx)
    MeMin(MeMin_idx,(size(ME,2)+1):end) = Min2(Min2_idx,1:size(Mineral,2));
    end  
    
end

MeMin.sample_ID = string(MeMin.sample_ID);
[rest, restPos] = setdiff(Min2.sample_ID,MeMin.sample_ID);
MeMin = [MeMin; Min2(restPos,:)];

MeMin.Sample_ID(size(ME,1):end) = MeMin.sample_ID(size(ME,1):end);
MeMin.Location_ID(size(ME,1):end) = MeMin.location_ID(size(ME,1):end);
MeMin.Basin_ID(size(ME,1):end) = MeMin.basin_ID(size(ME,1):end);
MeMin.SeaCat(size(ME,1):end) = MeMin.seaCat(size(ME,1):end);
MeMin.Observationtype(size(ME,1):end) = MeMin.observationtype(size(ME,1):end);
MeMin.Original_UnitME(size(ME,1):end) = MeMin.original_UnitMin(size(ME,1):end);
MeMin.Method(size(ME,1):end) = MeMin.MethodMin(size(ME,1):end);
MeMin.filtersize_m(size(ME,1):end) = MeMin.Filtersizem(size(ME,1):end);
MeMin.sievesize_m(size(ME,1):end) = MeMin.Sievesizem(size(ME,1):end);
MeMin.Dischargem3_s(size(ME,1):end) = MeMin.dischargem3_s(size(ME,1):end);
MeMin.TSS_mg_L(size(ME,1):end) = MeMin.TSS_mgL(size(ME,1):end);
MeMin.sand(size(ME,1):end) = MeMin.Sand(size(ME,1):end);
MeMin.silt(size(ME,1):end) = MeMin.Silt(size(ME,1):end);
MeMin.clay(size(ME,1):end) = MeMin.Clay(size(ME,1):end);
MeMin.avggrainsizem(size(ME,1):end) = MeMin.Avggrainsizem(size(ME,1):end);

%% River with information on sediment flux
%subset with only reduced variables
MeMint = MeMin;
MeMint = removevars(MeMint, {'SeaCat','Observationtype','sampletype','Original_UnitME','Treatment','Method','SO3_wt','N_orgwt','C_orgwt','C_inorgwt','C_totwt','Camol_L','Mgmol_L','Kmol_L','Namol_L','Clmol_L','Simol_L','DICmol_L','DOCmol_L','SO42mol_L','HCO3mol_L','T_waterDegC','pH','Alk','Cond','SICal','sample_ID','location_ID','seaCat','observationtype','Sampletype','basin_ID','original_UnitMin','MethodMin','rep_ID','Filtersizem','Sievesizem','dischargem3_s','TSS_mgL','Sand','Silt','Clay','Avggrainsizem','Quartz','Albite','Plg','Anorthite','Kfeldspar','Feldspar','Ampibole','Pyroxene','basGlass','totalmafic','Olivine','Mica','Smectite','Illite','Kaolinite','Chlorite','Gibbsite','Vermiculite','bulkclay','Calcite','Dolomite','bulkcarbonate','Casulfate','Fesulfide','FeOxide','FeHydroxide','Lvf','Lvm','Lv','Lp','Lcd','Lcc','Lch','Ls','Lms','Lmb','Lmf','Lm','Lu','L','totHM','Opq','Ap','Ep','Grt','Ttn','Zirc','Tourm','Ru','TiOx','Px','Amp','Ol','Sp','Mon','Ala','Cltd','Staur','And','Ky','Sil','trHM'});

%SiO2
%The procedure is commented in detail here and repeated similarly for the following elements
%Isolate values for which discharge, total suspended sediment concentration and SiO2 concentration are available
Sitqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.SiO2_wt));
SiO2tqs = MeMint(Sitqs,:);
%calculate the SiO2 flux for each entry
SiO2tqs.Qs_mgs = (SiO2tqs.Dischargem3_s./1000).*SiO2tqs.TSS_mg_L;
SiO2tqs.fSiO2_mgs = SiO2tqs.Qs_mgs.*(SiO2tqs.SiO2_wt./100);

%get each unique location id
LocationsSi = unique(SiO2tqs.Location_ID);
for i=1:size(LocationsSi,1) %iterate through unique locations
    loc_idx = find(SiO2tqs.Location_ID == LocationsSi(i)); %find all samples from the ith location
    fwm = (nansum(SiO2tqs.fSiO2_mgs(loc_idx))/nansum(SiO2tqs.Qs_mgs(loc_idx)))*100; %calculate sediment flux weighted mean (fwm)
    n = size(loc_idx,1); %how many observations used
    AbsMaxDiff =  max(abs(fwm - SiO2tqs.SiO2_wt(loc_idx))); %calculate maximum difference between fwm and any observation in the series
    RelMaxDiff =  100*(AbsMaxDiff/fwm); % convert to % in respect to fwm
    AbsmeanDiff =  meanOmitNan(abs(fwm - SiO2tqs.SiO2_wt(loc_idx))); %calculate the mean difference between fwm and any observation in the series
    AbsMedDiff =  medOmitNan(abs(fwm - SiO2tqs.SiO2_wt(loc_idx))); %calculate the median difference between fwm and any observation in the series
    RelMeanDiff =  100*(AbsMedDiff/fwm); % convert to % in respect to fwm
    FWM_SiO2(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; %summarize in table
end
%remove locations with less than 10 observations per location (small series and single values)
sgl = find(FWM_SiO2(:,7) <= 10);
FWM_SiO2(sgl,:) = [];
LocationsSi(sgl) = [];
%format to table
FWM_SiO2 = array2table(FWM_SiO2);
LocationsSi = array2table(LocationsSi);
%merge IDs and uncertainty maximum and median
FWM_SiO2 = [LocationsSi FWM_SiO2];
FWM_SiO2.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];

%Al2O3
Altqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.Al2O3_wt));
Al2O3tqs = MeMint(Altqs,:);
Al2O3tqs.Qs_mgs = (Al2O3tqs.Dischargem3_s./1000).*Al2O3tqs.TSS_mg_L;
Al2O3tqs.fAl2O3_mgs = Al2O3tqs.Qs_mgs.*(Al2O3tqs.Al2O3_wt./100);

LocationsAl = unique(Al2O3tqs.Location_ID);
for i=1:size(LocationsAl,1) 
    loc_idx = find(Al2O3tqs.Location_ID == LocationsAl(i)); 
    fwm = (nansum(Al2O3tqs.fAl2O3_mgs(loc_idx))/nansum(Al2O3tqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - Al2O3tqs.Al2O3_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - Al2O3tqs.Al2O3_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - Al2O3tqs.Al2O3_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_Al2O3(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_Al2O3(:,7) <= 10);
FWM_Al2O3(sgl,:) = [];
LocationsAl(sgl) = [];
FWM_Al2O3 = array2table(FWM_Al2O3);
LocationsAl = array2table(LocationsAl);
FWM_Al2O3 = [LocationsAl FWM_Al2O3];
FWM_Al2O3.Properties.VariableNames =  ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


%Fe2O3T
Fetqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.Fe2O3T_wt));
Fe2O3Ttqs = MeMint(Fetqs,:);
Fe2O3Ttqs.Qs_mgs = (Fe2O3Ttqs.Dischargem3_s./1000).*Fe2O3Ttqs.TSS_mg_L;
Fe2O3Ttqs.fFe2O3T_mgs = Fe2O3Ttqs.Qs_mgs.*(Fe2O3Ttqs.Fe2O3T_wt./100);

LocationsFe = unique(Fe2O3Ttqs.Location_ID);
for i=1:size(LocationsFe,1) 
    loc_idx = find(Fe2O3Ttqs.Location_ID == LocationsFe(i)); 
    fwm = (nansum(Fe2O3Ttqs.fFe2O3T_mgs(loc_idx))/nansum(Fe2O3Ttqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - Fe2O3Ttqs.Fe2O3T_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - Fe2O3Ttqs.Fe2O3T_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - Fe2O3Ttqs.Fe2O3T_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_Fe2O3T(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_Fe2O3T(:,7) <= 10);
FWM_Fe2O3T(sgl,:) = [];
LocationsFe(sgl) = [];
FWM_Fe2O3T = array2table(FWM_Fe2O3T);
LocationsFe = array2table(LocationsFe);
FWM_Fe2O3T = [LocationsFe FWM_Fe2O3T];
FWM_Fe2O3T.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


%MnO
Mntqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.MnO_wt));
MnOtqs = MeMint(Mntqs,:);
MnOtqs.Qs_mgs = (MnOtqs.Dischargem3_s./1000).*MnOtqs.TSS_mg_L;
MnOtqs.fMnO_mgs = MnOtqs.Qs_mgs.*(MnOtqs.MnO_wt./100);

LocationsMn = unique(MnOtqs.Location_ID);
for i=1:size(LocationsMn,1) 
    loc_idx = find(MnOtqs.Location_ID == LocationsMn(i)); 
    fwm = (nansum(MnOtqs.fMnO_mgs(loc_idx))/nansum(MnOtqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - MnOtqs.MnO_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - MnOtqs.MnO_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - MnOtqs.MnO_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_MnO(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_MnO(:,7) <= 10);
FWM_MnO(sgl,:) = [];
LocationsMn(sgl) = [];
FWM_MnO = array2table(FWM_MnO);
LocationsMn = array2table(LocationsMn);
FWM_MnO = [LocationsMn FWM_MnO];
FWM_MnO.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];

%CaO
Catqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.CaO_wt));
CaOtqs = MeMint(Catqs,:);
CaOtqs.Qs_mgs = (CaOtqs.Dischargem3_s./1000).*CaOtqs.TSS_mg_L;
CaOtqs.fCaO_mgs = CaOtqs.Qs_mgs.*(CaOtqs.CaO_wt./100);

LocationsCa = unique(CaOtqs.Location_ID);
for i=1:size(LocationsCa,1) 
    loc_idx = find(CaOtqs.Location_ID == LocationsCa(i)); 
    fwm = (nansum(CaOtqs.fCaO_mgs(loc_idx))/nansum(CaOtqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - CaOtqs.CaO_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - CaOtqs.CaO_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - CaOtqs.CaO_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_CaO(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_CaO(:,7) <= 10);
FWM_CaO(sgl,:) = [];
LocationsCa(sgl) = [];
FWM_CaO = array2table(FWM_CaO);
LocationsCa = array2table(LocationsCa);
FWM_CaO = [LocationsCa FWM_CaO];
FWM_CaO.Properties.VariableNames =  ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];

%MgO
Mgtqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.MgO_wt));
MgOtqs = MeMint(Mgtqs,:);
MgOtqs.Qs_mgs = (MgOtqs.Dischargem3_s./1000).*MgOtqs.TSS_mg_L;
MgOtqs.fMgO_mgs = MgOtqs.Qs_mgs.*(MgOtqs.MgO_wt./100);

LocationsMg = unique(MgOtqs.Location_ID);
for i=1:size(LocationsMg,1) 
    loc_idx = find(MgOtqs.Location_ID == LocationsMg(i)); 
    fwm = (nansum(MgOtqs.fMgO_mgs(loc_idx))/nansum(MgOtqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - MgOtqs.MgO_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - MgOtqs.MgO_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - MgOtqs.MgO_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_MgO(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_MgO(:,7) <= 10);
FWM_MgO(sgl,:) = [];
LocationsMg(sgl) = [];
FWM_MgO = array2table(FWM_MgO);
LocationsMg = array2table(LocationsMg);
FWM_MgO = [LocationsMg FWM_MgO];
FWM_MgO.Properties.VariableNames =  ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


%K2O
K2tqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.K2O_wt));
K2Otqs = MeMint(K2tqs,:);
K2Otqs.Qs_mgs = (K2Otqs.Dischargem3_s./1000).*K2Otqs.TSS_mg_L;
K2Otqs.fK2O_mgs = K2Otqs.Qs_mgs.*(K2Otqs.K2O_wt./100);

LocationsK2 = unique(K2Otqs.Location_ID);
for i=1:size(LocationsK2,1) 
    loc_idx = find(K2Otqs.Location_ID == LocationsK2(i)); 
    fwm = (nansum(K2Otqs.fK2O_mgs(loc_idx))/nansum(K2Otqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - K2Otqs.K2O_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - K2Otqs.K2O_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - K2Otqs.K2O_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_K2O(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_K2O(:,7) <= 10);
FWM_K2O(sgl,:) = [];
LocationsK2(sgl) = [];
FWM_K2O = array2table(FWM_K2O);
LocationsK2 = array2table(LocationsK2);
FWM_K2O = [LocationsK2 FWM_K2O];
FWM_K2O.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


%Na2O
Na2tqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.Na2O_wt));
Na2Otqs = MeMint(Na2tqs,:);
Na2Otqs.Qs_mgs = (Na2Otqs.Dischargem3_s./1000).*Na2Otqs.TSS_mg_L;
Na2Otqs.fNa2O_mgs = Na2Otqs.Qs_mgs.*(Na2Otqs.Na2O_wt./100);

LocationsNa2 = unique(Na2Otqs.Location_ID);
for i=1:size(LocationsNa2,1) 
    loc_idx = find(Na2Otqs.Location_ID == LocationsNa2(i)); 
    fwm = (nansum(Na2Otqs.fNa2O_mgs(loc_idx))/nansum(Na2Otqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - Na2Otqs.Na2O_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - Na2Otqs.Na2O_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - Na2Otqs.Na2O_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_Na2O(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_Na2O(:,7) <= 10);
FWM_Na2O(sgl,:) = [];
LocationsNa2(sgl) = [];
FWM_Na2O = array2table(FWM_Na2O);
LocationsNa2 = array2table(LocationsNa2);
FWM_Na2O = [LocationsNa2 FWM_Na2O];
FWM_Na2O.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


%TiO2
Titqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.TiO2_wt));
TiO2tqs = MeMint(Titqs,:);
TiO2tqs.Qs_mgs = (TiO2tqs.Dischargem3_s./1000).*TiO2tqs.TSS_mg_L;
TiO2tqs.fTiO2_mgs = TiO2tqs.Qs_mgs.*(TiO2tqs.TiO2_wt./100);

LocationsTi = unique(TiO2tqs.Location_ID);
for i=1:size(LocationsTi,1) 
    loc_idx = find(TiO2tqs.Location_ID == LocationsTi(i)); 
    fwm = (nansum(TiO2tqs.fTiO2_mgs(loc_idx))/nansum(TiO2tqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - TiO2tqs.TiO2_wt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - TiO2tqs.TiO2_wt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - TiO2tqs.TiO2_wt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_TiO2(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_TiO2(:,7) <= 10);
FWM_TiO2(sgl,:) = [];
LocationsTi(sgl) = [];
FWM_TiO2 = array2table(FWM_TiO2);
LocationsTi = array2table(LocationsTi);
FWM_TiO2 = [LocationsTi FWM_TiO2];
FWM_TiO2.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];



%P2O5
Ptqs = find(~isnan(MeMint.Dischargem3_s) & ~isnan(MeMint.TSS_mg_L) & ~isnan(MeMint.P2O5_totwt));
P2O5tqs = MeMint(Ptqs,:);
P2O5tqs.Qs_mgs = (P2O5tqs.Dischargem3_s./1000).*P2O5tqs.TSS_mg_L;
P2O5tqs.fP2O5_mgs = P2O5tqs.Qs_mgs.*(P2O5tqs.P2O5_totwt./100);

LocationsP = unique(P2O5tqs.Location_ID);
for i=1:size(LocationsP,1) 
    loc_idx = find(P2O5tqs.Location_ID == LocationsP(i)); 
    fwm = (nansum(P2O5tqs.fP2O5_mgs(loc_idx))/nansum(P2O5tqs.Qs_mgs(loc_idx)))*100; 
    n = size(loc_idx,1);
    AbsMaxDiff =  max(abs(fwm - P2O5tqs.P2O5_totwt(loc_idx))); 
    RelMaxDiff =  100*(AbsMaxDiff/fwm); 
    AbsmeanDiff =  meanOmitNan(abs(fwm - P2O5tqs.P2O5_totwt(loc_idx)));
    AbsMedDiff =  medOmitNan(abs(fwm - P2O5tqs.P2O5_totwt(loc_idx))); 
    RelMeanDiff =  100*(AbsMedDiff/fwm); 
    FWM_P2O5(i,1:7) = [fwm AbsMaxDiff RelMaxDiff AbsmeanDiff AbsMedDiff RelMeanDiff n]; 
end
sgl = find(FWM_P2O5(:,7) <= 10);
FWM_P2O5(sgl,:) = [];
LocationsP(sgl) = [];
FWM_P2O5 = array2table(FWM_P2O5);
LocationsP = array2table(LocationsP);
FWM_P2O5 = [LocationsP FWM_P2O5];
FWM_P2O5.Properties.VariableNames = ["Location_ID", "fwm", "AbsMaxDiff", "RelMaxDiff", "AbsMeanDiff", "AbsMedDiff", "RelMeanDiff", "N"];


% calculate median (log-normal) of median differences between fwm and any observation in the series for each element
MedDiffSiO2 = medOmitNan(FWM_SiO2.AbsMedDiff);
MedDiffAl2O3 = medOmitNan(FWM_Al2O3.AbsMedDiff);
MedDiffFe2O3T = medOmitNan(FWM_Fe2O3T.AbsMedDiff);
MedDiffMnO = medOmitNan(FWM_MnO.AbsMedDiff);
MedDiffCaO = medOmitNan(FWM_CaO.AbsMedDiff);
MedDiffMgO = medOmitNan(FWM_MgO.AbsMedDiff);
MedDiffK2O = medOmitNan(FWM_K2O.AbsMedDiff);
MedDiffNa2O = medOmitNan(FWM_Na2O.AbsMedDiff);
MedDiffTiO2 = medOmitNan(FWM_TiO2.AbsMedDiff);
MedDiffP2O5 = medOmitNan(FWM_P2O5.AbsMedDiff);

% calculate mean of mean differences between fwm and any observation in the series for each element
meanDiffSiO2 = meanOmitNan(FWM_SiO2.AbsMeanDiff);
meanDiffAl2O3 = meanOmitNan(FWM_Al2O3.AbsMeanDiff);
meanDiffFe2O3T = meanOmitNan(FWM_Fe2O3T.AbsMeanDiff);
meanDiffMnO = meanOmitNan(FWM_MnO.AbsMeanDiff);
meanDiffCaO = meanOmitNan(FWM_CaO.AbsMeanDiff);
meanDiffMgO = meanOmitNan(FWM_MgO.AbsMeanDiff);
meanDiffK2O = meanOmitNan(FWM_K2O.AbsMeanDiff);
meanDiffNa2O = meanOmitNan(FWM_Na2O.AbsMeanDiff);
meanDiffTiO2 = meanOmitNan(FWM_TiO2.AbsMeanDiff);
meanDiffP2O5 = meanOmitNan(FWM_P2O5.AbsMeanDiff);


% calculate mean of Max differences between fwm and any observation in the series for each element
MaxDiffSiO2 = meanOmitNan(FWM_SiO2.AbsMaxDiff);
MaxDiffAl2O3 = meanOmitNan(FWM_Al2O3.AbsMaxDiff);
MaxDiffFe2O3T = meanOmitNan(FWM_Fe2O3T.AbsMaxDiff);
MaxDiffMnO = meanOmitNan(FWM_MnO.AbsMaxDiff);
MaxDiffCaO = meanOmitNan(FWM_CaO.AbsMaxDiff);
MaxDiffMgO = meanOmitNan(FWM_MgO.AbsMaxDiff);
MaxDiffK2O = meanOmitNan(FWM_K2O.AbsMaxDiff);
MaxDiffNa2O = meanOmitNan(FWM_Na2O.AbsMaxDiff);
MaxDiffTiO2 = meanOmitNan(FWM_TiO2.AbsMaxDiff);
MaxDiffP2O5 = meanOmitNan(FWM_P2O5.AbsMaxDiff);

% calculate mean of rel Max differences between fwm and any observation in the series for each element
rMaxDiffSiO2 = meanOmitNan(FWM_SiO2.RelMaxDiff);
rMaxDiffAl2O3 = meanOmitNan(FWM_Al2O3.RelMaxDiff);
rMaxDiffFe2O3T = meanOmitNan(FWM_Fe2O3T.RelMaxDiff);
rMaxDiffMnO = meanOmitNan(FWM_MnO.RelMaxDiff);
rMaxDiffCaO = meanOmitNan(FWM_CaO.RelMaxDiff);
rMaxDiffMgO = meanOmitNan(FWM_MgO.RelMaxDiff);
rMaxDiffK2O = meanOmitNan(FWM_K2O.RelMaxDiff);
rMaxDiffNa2O = meanOmitNan(FWM_Na2O.RelMaxDiff);
rMaxDiffTiO2 = meanOmitNan(FWM_TiO2.RelMaxDiff);
rMaxDiffP2O5 = meanOmitNan(FWM_P2O5.RelMaxDiff);

%summarize on table
MedDiff = [MedDiffSiO2 MedDiffAl2O3 MedDiffFe2O3T MedDiffMnO MedDiffCaO MedDiffMgO MedDiffK2O MedDiffNa2O MedDiffTiO2 MedDiffP2O5];
MeanDiff = [meanDiffSiO2 meanDiffAl2O3 meanDiffFe2O3T meanDiffMnO meanDiffCaO meanDiffMgO meanDiffK2O meanDiffNa2O meanDiffTiO2 meanDiffP2O5];
MaxDiff = [MaxDiffSiO2 MaxDiffAl2O3 MaxDiffFe2O3T MaxDiffMnO MaxDiffCaO MaxDiffMgO MaxDiffK2O MaxDiffNa2O MaxDiffTiO2 MaxDiffP2O5];
rMaxDiff = [rMaxDiffSiO2 rMaxDiffAl2O3 rMaxDiffFe2O3T rMaxDiffMnO rMaxDiffCaO rMaxDiffMgO rMaxDiffK2O rMaxDiffNa2O rMaxDiffTiO2 rMaxDiffP2O5];




