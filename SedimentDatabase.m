%% Import Data

%save('SedimentDatabase', 'ID', 'Location', 'Ref', 'ME','TE', 'Mineral', 'RiverNames')
%import database
load SedimentDatabase



%% Figure visualising database extent - ROBINSON PROJECTION
Location.Content = categorical(Location.Content);
geobubble(Location.Lat_deg,Location.Lon_deg,Location.Observations,Location.Content)


%% Baisn-wise average concentration representative for export to the ocean
%remove coarse grained bedload samples (> 125 µm sieve size) 
Ssizes = unique(ME.sievesize_m(:));
bed = (find(ME.sampletype == 'bs' & isnan(ME.sievesize_m(:)) | ME.sievesize_m(:) >= 125));
ME(bed,:) = [];
%cut values that are representative for export to the ocean (Rep_ID = 1)
Rep_ME = ME(ME.Rep_ID ==1,:);
AnAvg = Rep_ME(Rep_ME.Observationtype == "an",:); %cut annual averages in separate table
SeaAvg = Rep_ME(Rep_ME.Observationtype == "sea",:); %cut annual averages in separate table
Rep_ME = Rep_ME(Rep_ME.Observationtype == "single" | Rep_ME.Observationtype == "sa",:);

%Make Basin_ID and SeaCat categorical variables
Rep_ME.Basin_ID = categorical(Rep_ME.Basin_ID);
SeaAvg.Basin_ID = categorical(SeaAvg.Basin_ID);
AnAvg.Basin_ID = categorical(AnAvg.Basin_ID);

%format seasonal and annual averages for later use
BasinAn = AnAvg.Basin_ID;
AnAvg = AnAvg(:,19:31);
AnAvg.Basin_ID = BasinAn;

BasinSeaAvg = SeaAvg.Basin_ID;
SeasonSeaAvg = SeaAvg.SeaCat;
SeaAvg = SeaAvg(:,19:31);
SeaAvg.Basin_ID = BasinSeaAvg;
SeaAvg.SeaCat = SeasonSeaAvg;

%define fuctions that calculate median, mean and standard deviation omitting NaN values
%Median is preferred for all elements except for Si, Al and Mn, because data distribution is mostly skewed(non-Gaussian)
%check: open Rep_ME and check histograms via 'Plots' 
medOmitNan = @(x) median(x,'omitNan');
meanOmitNan = @(x) mean(x,'omitNan'); 
SDOmitNan = @(x) std(x,0,'omitNan');

% Average season- and basin-wise using 'splitapply()' workflow
%(https://de.mathworks.com/help/matlab/ref/splitapply.html)


%single observations

% find each unique group in defined by basin and season
[BasinGroup, Basin] = findgroups(Rep_ME.Basin_ID);

% Median
MedianperBasin = [];
for i = 19:31
    y = splitapply(medOmitNan,Rep_ME(:,i),BasinGroup);
    MedianperBasin(:,i) = y;
end
MedianperBasin(:,1:18) = []; %remove zero entries (variable not averaged)
MedianperBasin = array2table(MedianperBasin);
MedianperBasin.Basin_ID = Basin; %recombine with IDs

% Mean
MeanperBasin = [];
for i = 19:31
    y = splitapply(meanOmitNan,Rep_ME(:,i),BasinGroup);
    MeanperBasin(:,i) = y;
end
MeanperBasin(:,1:18) = []; %remove zero entries (variable not averaged)
MeanperBasin = array2table(MeanperBasin);
MeanperBasin.Basin_ID = Basin;%recombine with IDs

% Standard deviation
SDperBasin = [];
for i = 19:31
    y = splitapply(SDOmitNan,Rep_ME(:,i),BasinGroup);
    SDperBasin(:,i) = y;
end
SDperBasin(:,1:18) = []; %remove zero entries (what was not averaged)
SDperBasin = array2table(SDperBasin);
SDperBasin.Basin_ID = Basin;%recombine with IDs

%name variables
MedianperBasin.Properties.VariableNames = [AnAvg.Properties.VariableNames];
MeanperBasin.Properties.VariableNames = [AnAvg.Properties.VariableNames];
SDperBasin.Properties.VariableNames = [AnAvg.Properties.VariableNames];


%seasonal averages
%grouping
[BasinGroupSea, BasinSea] = findgroups(SeaAvg.Basin_ID);

% Median
MedianperBasinSea = [];
for i = 1:13
    y = splitapply(medOmitNan,SeaAvg(:,i),BasinGroupSea);
    MedianperBasinSea(:,i) = y;
end
MedianperBasinSea = array2table(MedianperBasinSea);
MedianperBasinSea.Basin_ID = BasinSea; %recombine with IDs

% Mean
MeanperBasinSea = [];
for i = 1:13
    y = splitapply(meanOmitNan,SeaAvg(:,i),BasinGroupSea);
    MeanperBasinSea(:,i) = y;
end
MeanperBasinSea = array2table(MeanperBasinSea);
MeanperBasinSea.Basin_ID = BasinSea;%recombine with IDs

%name variables
MedianperBasinSea.Properties.VariableNames = [AnAvg.Properties.VariableNames];
MedianperBasinSea = splitvars(MedianperBasinSea);
MeanperBasinSea = splitvars(MeanperBasinSea);
MeanperBasinSea.Properties.VariableNames = [AnAvg.Properties.VariableNames];


%annual averages

%recombine annual scale data
AllMed = [MedianperBasin; MedianperBasinSea; AnAvg];
AllMean = [MeanperBasin; MeanperBasinSea; AnAvg];

%group and average annual values as above
[BasinGroupAn, BasinAn] = findgroups(AllMed.Basin_ID);

% Median
MedianAnnual = [];
for i = 1:13
    y = splitapply(medOmitNan,AllMed(:,i),BasinGroupAn);
    MedianAnnual(:,i) = y;
end
MedianAnnual = array2table(MedianAnnual);
MedianAnnual.Basin_ID = BasinAn; %recombine with IDs

% Mean
MeanAnnual = [];
for i = 1:13
    y = splitapply(meanOmitNan,AllMean(:,i),BasinGroupAn);
    MeanAnnual(:,i) = y;
end
MeanAnnual = array2table(MeanAnnual);
MeanAnnual.Basin_ID = BasinAn;%recombine with IDs

%Rename Vars
MeanAnnual.Properties.VariableNames = [ "SiO2mean" "Al2O3mean" "Fe2O3Tmean" "MnOmean" "CaOmean" "MgOmean" "K2Omean" "Na2Omean" "TiO2mean" "PorgMean" "PinorgMean" "P2O5Tmean" "LOImean" "Basin_ID"];
MedianAnnual.Properties.VariableNames = [ "SiO2med" "Al2O3med" "Fe2O3Tmed" "MnOmed" "CaOmed" "MgOmed" "K2Omed" "Na2Omed" "TiO2med" "Porgmed" "Pinorgmed" "P2O5Tmed" "LOImed" "Basin_ID"];

%merge median, mean and standard deviation
AvgperBasin = [MedianAnnual(:,1:13) MeanAnnual(:,1:14)];
AvgperBasin = removevars(AvgperBasin, {'Porgmed','Pinorgmed','PorgMean','PinorgMean'}); %remove variables only containing Nan, i.e. no values with Rep_ID = 1 are available

%assigning names to each Basin_ID
% make Basin_IDs strings to compare in both tables
AvgperBasin.Basin_ID = string(AvgperBasin.Basin_ID);
RiverNames.Properties.VariableNames{3} = 'Bsn_ID';
RiverNames.Bsn_ID = string(RiverNames.Bsn_ID);

%find common Basin_IDs
[commonBasinID,~] = intersect(AvgperBasin.Basin_ID,RiverNames.Bsn_ID);

for i=1:size(commonBasinID,1)
    
    name_idx=find(RiverNames.Bsn_ID==commonBasinID(i)); %find position of similar entries in RiverNames table
    Avg_idx=find(AvgperBasin.Basin_ID==commonBasinID(i)); %find position of similar entries in Basin-wise median table
    
    if size(name_idx,1) ~=1 %this if-loop is here because there are rows in the TSS table that have the same ID#
    name_idx=name_idx(1,1);
    end
    
    if ~isempty(name_idx) && ~isempty(Avg_idx)
    AvgperBasin(Avg_idx,24) = RiverNames(name_idx,1); %append river name, if Basin_IDs are similar
    end  
    
end
AvgperBasin.Properties.VariableNames{24} = 'RiverNames'; %rename added variable

%note that the standard deviations hold many 0 entries, that is for 
%these basins there is only one value available that is thus equal to the mean
%In the main text, the standard deviation will be evaluated using
%single values representing inter-sample variability the best

%% Merging with Sediment fluxes and calculating exported mass of each major element


%import annual mean sediment fluxes and river names of Milliman & Farnsworth (2011): River
%discharge to the global ocean. Cambridge University Press, DOI:10.1017/CBO9780511781247.
load('TSS_MF')

%find common names
[commonNames,~] = intersect(AvgperBasin.RiverNames,TSS_MF.RiverName); % 129 matches with the same name

%merge both databases by finding and using indices of common names in both table
for i=1:size(commonNames,1)
    
    TSS_idx=find(TSS_MF.RiverName == commonNames(i));
    avg_idx=find(AvgperBasin.RiverNames == commonNames(i));
   
    if size(TSS_idx,1) ~=1 %this if-loop is here because there are rows in the TSS table that have the same ID#
        TSS_idx=TSS_idx(1,1);
    end

    if ~isempty(TSS_idx) && ~isempty(avg_idx)
     AvgperBasin(avg_idx,25) = TSS_MF(TSS_idx,3);
    end
    
end
AvgperBasin.Properties.VariableNames{25} = 'TSS_Mt_a';

%calculate solid oxide fluxes based mean and median concentration (depending on distribution of specific element)
AvgperBasin.fSiO2 = (AvgperBasin.SiO2mean./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fAl2O3 = (AvgperBasin.Al2O3mean./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fFe2O3T = (AvgperBasin.Fe2O3Tmed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fMnO = (AvgperBasin.MnOmean./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fCaO = (AvgperBasin.CaOmed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fMgO = (AvgperBasin.MgOmed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fK2O = (AvgperBasin.K2Omed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fNa2O = (AvgperBasin.Na2Omed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fTiO2 = (AvgperBasin.TiO2med./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fP2O5_tot = (AvgperBasin.P2O5Tmed./100).* AvgperBasin.TSS_Mt_a;
AvgperBasin.fLOI = (AvgperBasin.LOImed./100).* AvgperBasin.TSS_Mt_a;

%Global mean composition
GlobalAvg = table();
GlobalAvg.SiO2_mean = meanOmitNan(AvgperBasin.SiO2mean);
GlobalAvg.Al2O3_mean = meanOmitNan(AvgperBasin.Al2O3mean);
GlobalAvg.Fe2O3_mean = meanOmitNan(AvgperBasin.Fe2O3Tmean);
GlobalAvg.MnO_mean = meanOmitNan(AvgperBasin.MnOmean);
GlobalAvg.CaO_mean = meanOmitNan(AvgperBasin.CaOmean);
GlobalAvg.MgO_mean = meanOmitNan(AvgperBasin.MgOmean);
GlobalAvg.K2O_mean = meanOmitNan(AvgperBasin.K2Omean);
GlobalAvg.Na2O_mean = meanOmitNan(AvgperBasin.Na2Omean);
GlobalAvg.TiO2_mean = meanOmitNan(AvgperBasin.TiO2mean);
GlobalAvg.P2O5_mean = meanOmitNan(AvgperBasin.P2O5Tmean);
GlobalAvg.LOI_mean = meanOmitNan(AvgperBasin.LOImean);
%median
GlobalAvg.SiO2_med = medOmitNan(AvgperBasin.SiO2med);
GlobalAvg.Al2O3_med = medOmitNan(AvgperBasin.Al2O3med);
GlobalAvg.Fe2O3_med = medOmitNan(AvgperBasin.Fe2O3Tmed);
GlobalAvg.MnO_med = medOmitNan(AvgperBasin.MnOmed);
GlobalAvg.CaO_med = medOmitNan(AvgperBasin.CaOmed);
GlobalAvg.MgO_med = medOmitNan(AvgperBasin.MgOmed);
GlobalAvg.K2O_med = medOmitNan(AvgperBasin.K2Omed);
GlobalAvg.Na2O_med = medOmitNan(AvgperBasin.Na2Omed);
GlobalAvg.TiO2_med = medOmitNan(AvgperBasin.TiO2med);
GlobalAvg.P2O5_med = medOmitNan(AvgperBasin.P2O5Tmed);
GlobalAvg.LOI_med = meanOmitNan(AvgperBasin.LOImed);
% sample SD
GlobalAvg.SiO2_SD = SDOmitNan(AvgperBasin.SiO2med);
GlobalAvg.Al2O3_SD = SDOmitNan(AvgperBasin.Al2O3med);
GlobalAvg.Fe2O3_SD = SDOmitNan(AvgperBasin.Fe2O3Tmed);
GlobalAvg.MnO_SD = SDOmitNan(AvgperBasin.MnOmed);
GlobalAvg.CaO_SD = SDOmitNan(AvgperBasin.CaOmed);
GlobalAvg.MgO_SD = SDOmitNan(AvgperBasin.MgOmed);
GlobalAvg.K2O_SD = SDOmitNan(AvgperBasin.K2Omed);
GlobalAvg.Na2O_SD = SDOmitNan(AvgperBasin.Na2Omed);
GlobalAvg.TiO2_SD = SDOmitNan(AvgperBasin.TiO2med);
GlobalAvg.P2O5_SD = SDOmitNan(AvgperBasin.P2O5Tmed);
GlobalAvg.LOI_SD = SDOmitNan(AvgperBasin.LOImed);

% Error estimate with penalty for uncertainty from time-series and sample SD
GlobalAvg.SiO2_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.SiO2med)),1)))^2+(8.6^2));
GlobalAvg.Al2O3_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.Al2O3med)),1)))^2+(5.2^2));
GlobalAvg.Fe2O3_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.Fe2O3Tmed)),1)))^2+(1.5^2));
GlobalAvg.MnO_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.MnOmed)),1)))^2+(0.1^2));
GlobalAvg.CaO_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.CaOmed)),1)))^2+(2.9^2));
GlobalAvg.MgO_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.MgOmed)),1)))^2+(0.5^2));
GlobalAvg.K2O_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.K2Omed)),1)))^2+(0.6^2));
GlobalAvg.Na2O_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.Na2Omed)),1)))^2+(0.2^2));
GlobalAvg.TiO2_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.TiO2med)),1)))^2+(0.2^2));
GlobalAvg.P2O5_ER = sqrt((GlobalAvg.SiO2_SD/sqrt(size(find(~isnan(AvgperBasin.P2O5Tmed)),1)))^2+(0.1^2));


%calculate global sediment-flux weighted mean concentrations
%get indices where which oxide is available
Si_idx = find(~isnan(AvgperBasin.SiO2med));
Al_idx = find(~isnan(AvgperBasin.Al2O3med));
Fe_idx = find(~isnan(AvgperBasin.Fe2O3Tmed));
Mn_idx = find(~isnan(AvgperBasin.MnOmed));
Ca_idx = find(~isnan(AvgperBasin.CaOmed));
Mg_idx = find(~isnan(AvgperBasin.MgOmed));
Na_idx = find(~isnan(AvgperBasin.Na2Omed));
K_idx = find(~isnan(AvgperBasin.K2Omed));
Ti_idx = find(~isnan(AvgperBasin.TiO2med));
P_idx = find(~isnan(AvgperBasin.P2O5Tmed));
LOI_idx = find(~isnan(AvgperBasin.LOImed));

GlobalAvg.SiO2_wgtSed = (nansum(AvgperBasin.fSiO2(Si_idx))/nansum(AvgperBasin.TSS_Mt_a(Si_idx))).*100;
GlobalAvg.Al2O3_wgtSed = (nansum(AvgperBasin.fAl2O3(Al_idx))/nansum(AvgperBasin.TSS_Mt_a(Al_idx))).*100;
GlobalAvg.Fe2O3_wgtSed = (nansum(AvgperBasin.fFe2O3T(Fe_idx))/nansum(AvgperBasin.TSS_Mt_a(Fe_idx))).*100;
GlobalAvg.MnO_wgtSed = (nansum(AvgperBasin.fMnO(Mn_idx))/nansum(AvgperBasin.TSS_Mt_a(Mn_idx))).*100;
GlobalAvg.CaO_wgtSed = (nansum(AvgperBasin.fCaO(Ca_idx))/nansum(AvgperBasin.TSS_Mt_a(Ca_idx))).*100;
GlobalAvg.MgO_wgtSed = (nansum(AvgperBasin.fMgO(Mg_idx))/nansum(AvgperBasin.TSS_Mt_a(Mg_idx))).*100;
GlobalAvg.K2O_wgtSed = (nansum(AvgperBasin.fK2O(K_idx))/nansum(AvgperBasin.TSS_Mt_a(K_idx))).*100;
GlobalAvg.Na2O_wgtSed = (nansum(AvgperBasin.fNa2O(Na_idx))/nansum(AvgperBasin.TSS_Mt_a(Na_idx))).*100;
GlobalAvg.TiO2_wgtSed = (nansum(AvgperBasin.fTiO2(Ti_idx))/nansum(AvgperBasin.TSS_Mt_a(Ti_idx))).*100;
GlobalAvg.P2O5_wgtSed = (nansum(AvgperBasin.fP2O5_tot(P_idx))/nansum(AvgperBasin.TSS_Mt_a(P_idx))).*100;
GlobalAvg.LOI_wgtSed = (nansum(AvgperBasin.fLOI(LOI_idx))/nansum(AvgperBasin.TSS_Mt_a(LOI_idx))).*100;

%The table 'GlobalAvg' holds all global averages, while the table 'AvgperBasin' holds the mean, standard deviation 
%and median of each basin along with its name and annual sediment flux (as given by Milliman and Farnsworth 2011).

writetable(GlobalAvg, 'GlobalAverages_MEOxides.csv')
writetable(AvgperBasin, 'AveragePerBasin_MEOxides.csv')