%Primary calculation and analysis script
%Some elements (noted) are in the accompanying Python notebook
%Written for "Record-setting humid heat" manuscript
%Colin Raymond, 2024

canesm5dir='/Volumes/ExternalDriveD/CanESM5_LE_daily/';
miroc6dir='/Volumes/ExternalDriveD/MIROC6_LE_daily/';
mpigedir='/Volumes/ExternalDriveF/MPI-GE/';
era5dir_1='/Volumes/ExtDriveC/ERA5_Hourly_TTdData_to1999/';
era5dir_2='/Volumes/ExternalDriveZ/ERA5_Hourly_TTdData_since2000/';
jra55dir='/Volumes/ExternalDriveF/JRA55_daily/';
merra2dir='/Volumes/ExternalDriveF/MERRA2/majorvars_global_19802023/';
saveloc='/Volumes/ExternalDriveD/RecordSetting_Humid_Heat_savedarrays/';

icloud='~/Library/Mobile Documents/com~apple~CloudDocs/';
figloc=strcat(icloud,'General_Academics/Research/RecordSetting_Humid_Heat/');
format shortG;

reloadthings=0; %30 sec; must be run on start-up
defineregions=0; %30 sec; must be run on start-up
processcanesm5=0; %3 hr for calcmaxs (for main results); 1 hr for calcmeans
    if processcanesm5==1;calcmaxs_canesm5=0;calcmeans_canesm5=1;end
processmiroc6=0; %5 hr
    if processmiroc6==1;calcmaxs_miroc6=0;calcmeans_miroc6=1;end
processmpige=0; %2 hr 30 min
    if processmpige==1;calcmaxs_mpige=1;calcmeans_mpige=0;end
convertmodeldatatoregsandnormalize=0; %10 sec
rereadcanesm5_savedaily=0;
    canesm5_dopart1=0; %40 min total [then a 40-min Python interlude]
    canesm5_dopart2=1; %25 min total
rereadmiroc6_savedaily=0; %overall time: 7.5 min per ens mem
    miroc6_dopart1=0; %3 min per ens mem, 30 min total [then a 2.5-min-per-ens-mem Python interlude, 25 min total]
    miroc6_dopart2=0; %2 min per ens mem, 20 min total
rereadmpige_savedaily=0; %maxs only; 3 min per ensemble member
    nummpigeensmemstodo=30;
processera5=0; %6 hr
    latwgt_era5=1; %(lat-weighting) default is 1; originally 0
processjra55=0; %1 hr 15 min
    latwgt_jra55=1; %(lat-weighting) default is 1; originally 0
adjustgmst=0; %5 sec; always needed if models have been reprocessed
era5_calcstdanoms=0; %15 sec
jra55_calcstdanoms=0; %10 sec

completelyconsistenttopevents=0; %2 min; Figure S7
preparetable2=0; %15 sec; prepares data for Table 2
morestats=0; %10 sec; stats for Excel; for some SI tables
assessdiurnaleffects=0; %for Figure S8
createfig2=0; %1 min; creates Figure 2
createfig1=1; %30 sec -- creates Figure 1
    %panel a: as in former Fig 1a but with fewer exclusions (exclude region only if magnitudes differ by >=10%)
    %panel b: effect-of-current-record-on-return-period-of-itself map
fig1_siversions=0; %1 min; versions for T and JRA55
createfig3=0; %30 sec -- creates Figure 3
    %panel a: former Fig 1b
fig3_siversions=0; %1 min; versions for T and JRA55
fig4setup_clustering=0; %30 sec
fig4setup_duration=0; %1 sec
fig4setup_intensity=0; %15 sec
createfig4=0; %3.5 min -- creates Figure 4
findmostsingularevents=0; %1 sec; for Table 1


%COMPLEMENTARY ANALYSIS, DIRECTLY EXPANDING ON THOMPSON ET AL. 2023, IS DONE IN PYTHON
%(rshhanalysis.ipynb)

compareagainst='jra55';

numreg=237;
numreg_exclant=216; %excluding the 21 Antarctic regions

numyr=63; %1961-2023; so 9 years from SSP585 run, 54 from historical run
ny_fut=9;ny_hist=numyr-ny_fut; %future: 2015-2023
futfirstyr=1;futlastyr=ny_fut; %default -- listed as relative years, with 2015=1
firstyear=1961;lastyear=2023;
yearlist=1:numyr;

%Option to only add a few years, typically the most recent ones
dohist=1;dofut=1; %default is 1 for both
numyrtodo=numyr; %default is numyr, but can be smaller e.g. if adding only a few more
ny_hist_todo=ny_hist;ny_fut_todo=ny_fut; %defaults are ny_hist and ny_fut respectively
if numyrtodo<numyr;futfirstyr=7;futlastyr=9;end %if doing only a few years, specify exact ones

firstclimoyr=21;lastclimoyr=50; %climatologies are based on 1981-2010, as in Thompson et al. 2022
numclimoyr=lastclimoyr-firstclimoyr+1;

monthstarts=[1;32;60;91;121;152;182;213;244;274;305;335];
monthstops=[31;59;90;120;151;181;212;243;273;304;334;365];
monthlens=[31;28;31;30;31;30;31;31;30;31;30;31];

splabels={'a)','b)','c)','d)','e)','f)'};

modelnames={'canesm5';'miroc6';'mpige'};

%Set up lat/lon arrays for each dataset
exist latsz_era5;
if ans==0
    z_192x288=ncread(strcat(icloud,'General_Academics/Research/KeyFiles/elev_192x288.nc'),'orog');
    z_192x288=recenter(flipud(z_192x288'));
    lat1d=ncread(strcat(icloud,'General_Academics/Research/KeyFiles/elev_192x288.nc'),'lat');
    lon1d=ncread(strcat(icloud,'General_Academics/Research/KeyFiles/elev_192x288.nc'),'lon');
    [lat2d_192x288,lon2d_192x288]=latlon_2dfrom1d(lat1d,lon1d);
    lat2d_192x288=flipud(lat2d_192x288');lon2d_192x288=recenter(flipud(lon2d_192x288'));

    lat1d=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r1i1p2f1_gn_20150101-21001231.nc'),'lat');
    lon1d=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r1i1p2f1_gn_20150101-21001231.nc'),'lon');
    [lat2d_canesm5,lon2d_canesm5]=latlon_2dfrom1d(lat1d,lon1d);
    lat2d_canesm5=recenter(flipud(lat2d_canesm5'));lon2d_canesm5=recenter(flipud(lon2d_canesm5'));

    lat1d=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_historical_r1i1p1f1_gn_19900101-19991231.nc'),'lat');
    lon1d=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_historical_r1i1p1f1_gn_19900101-19991231.nc'),'lon');
    [lat2d_miroc6,lon2d_miroc6]=latlon_2dfrom1d(lat1d,lon1d);
    lat2d_miroc6=recenter(flipud(lat2d_miroc6'));lon2d_miroc6=recenter(flipud(lon2d_miroc6'));

    lat1d=ncread(strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_20100101-20141231.nc'),'lat');
    lon1d=ncread(strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_20100101-20141231.nc'),'lon');
    [lat2d_mpige,lon2d_mpige]=latlon_2dfrom1d(lat1d,lon1d);
    lat2d_mpige=recenter(flipud(lat2d_mpige'));lon2d_mpige=recenter(flipud(lon2d_mpige'));

    lat1d=ncread(strcat(era5dir_2,'t2m_2010.nc'),'latitude');
    lon1d=ncread(strcat(era5dir_2,'t2m_2010.nc'),'longitude');
    [lat2d_era5,lon2d_era5]=latlon_2dfrom1d(lat1d,lon1d);
    lon2d_era5=recenter(lon2d_era5');lat2d_era5=recenter(lat2d_era5');

    lat1d=ncread(strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.1961010100_1961123118.raymond741169.nc'),'g4_lat_1');
    lon1d=ncread(strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.1961010100_1961123118.raymond741169.nc'),'g4_lon_2');
    [lat2d_jra55,lon2d_jra55]=latlon_2dfrom1d(lat1d,lon1d);
    lon2d_jra55=recenter(lon2d_jra55');lat2d_jra55=recenter(lat2d_jra55');

    year=1991;namepart='100';thismon=1;thisdom=1;thismon_str=strcat('0',num2str(thismon));thisdom_str=strcat('0',num2str(thisdom));
    thisfilename=strcat(merra2dir,'MERRA2_',namepart,'.tavg1_2d_slv_Nx.',num2str(year),thismon_str,thisdom_str,'.SUB.nc');
    clear lat2d_merra2;clear lon2d_merra2;
    lattemp=ncread(thisfilename,'lat');lontemp=ncread(thisfilename,'lon');
    lat=double(lattemp);lon=double(lontemp);
    for row=1:size(lat,1)
        for col=1:size(lon,1)
            lat2d_merra2(row,col)=lat(row);lon2d_merra2(row,col)=lon(col);
        end
    end
    westernhem=lon2d_merra2>=180;lon2d_merra2(westernhem)=lon2d_merra2(westernhem)-360;lat2d_merra2=flipud(lat2d_merra2);


    latsz_era5=size(lat2d_era5,1);lonsz_era5=size(lat2d_era5,2);
    latsz_merra2=size(lat2d_merra2,1);lonsz_merra2=size(lat2d_merra2,2);
    latsz_jra55=size(lat2d_jra55,1);lonsz_jra55=size(lat2d_jra55,2);
end


if reloadthings==1
    f=load(strcat(saveloc,'canesm5output.mat'));tw_annmax_canesm5=f.tw_annmax_canesm5;
        canesm5_enssz=50;
    f=load(strcat(saveloc,'miroc6output.mat'));tw_annmax_miroc6=f.tw_annmax_miroc6;
        miroc6_enssz=10;
    f=load(strcat(saveloc,'mpigeoutput.mat'));tw_annmax_mpige=f.tw_annmax_mpige;
        mpige_enssz=50; %originally 30
    f=load(strcat(saveloc,'era5output.mat'));tw_annmax_era5=f.tw_annmax_era5;
    f=load(strcat(saveloc,'era5output_regions.mat'));
        tw_regdecstdanoms_era5=f.tw_regdecstdanoms_era5;tw_regdecstdanoms_adj_era5=f.tw_regdecstdanoms_adj_era5;
        t_regdecstdanoms_era5=f.t_regdecstdanoms_era5;t_regdecstdanoms_adj_era5=f.t_regdecstdanoms_adj_era5;
        q_regdecstdanoms_era5=f.q_regdecstdanoms_era5;q_regdecstdanoms_adj_era5=f.q_regdecstdanoms_adj_era5;
        tw_regdecvals_era5=f.tw_regdecvals_era5;tw_regdecvals_adj_era5=f.tw_regdecvals_adj_era5;
        t_regdecvals_era5=f.t_regdecvals_era5;t_regdecvals_adj_era5=f.t_regdecvals_adj_era5;
        q_regdecvals_era5=f.q_regdecvals_era5;
        tw_doyofregannmax_era5=f.tw_doyofregannmax_era5;
        yearofmax_tw_era5=f.yearofmax_tw_era5;yearofmax_tw_adj_era5=f.yearofmax_tw_adj_era5;
        yearofmax_t_era5=f.yearofmax_t_era5;yearofmax_t_adj_era5=f.yearofmax_t_adj_era5;
        yearofmax_q_era5=f.yearofmax_q_era5;yearofmax_q_adj_era5=f.yearofmax_q_adj_era5;
        doyofmax_tw_era5=f.doyofmax_tw_era5;doyofmax_tw_adj_era5=f.doyofmax_tw_adj_era5;
        doyofmax_t_era5=f.doyofmax_t_era5;doyofmax_t_adj_era5=f.doyofmax_t_adj_era5;
        doyofmax_q_era5=f.doyofmax_q_era5;doyofmax_q_adj_era5=f.doyofmax_q_adj_era5;
        maxstdanom_tw_era5=f.maxstdanom_tw_era5;maxstdanom_tw_adj_era5=f.maxstdanom_tw_adj_era5;
        maxstdanom_t_era5=f.maxstdanom_t_era5;maxstdanom_t_adj_era5=f.maxstdanom_t_adj_era5;
        maxstdanom_q_era5=f.maxstdanom_q_era5;maxstdanom_q_adj_era5=f.maxstdanom_q_adj_era5;
        peakmonth_tw_era5=f.peakmonth_tw_era5;peakmonth_t_era5=f.peakmonth_t_era5;peakmonth_q_era5=f.peakmonth_q_era5;
    f=load(strcat(saveloc,'merra2output_regions.mat'));
        tw_regdecstdanoms_merra2=f.tw_regdecstdanoms_merra2;tw_regdecstdanoms_adj_merra2=f.tw_regdecstdanoms_adj_merra2;
        t_regdecstdanoms_merra2=f.t_regdecstdanoms_merra2;t_regdecstdanoms_adj_merra2=f.t_regdecstdanoms_adj_merra2;
        q_regdecstdanoms_merra2=f.q_regdecstdanoms_merra2;q_regdecstdanoms_adj_merra2=f.q_regdecstdanoms_adj_merra2;
        tw_regdecvals_merra2=f.tw_regdecvals_merra2;tw_regdecvals_adj_merra2=f.tw_regdecvals_adj_merra2;
        t_regdecvals_merra2=f.t_regdecvals_merra2;t_regdecvals_adj_merra2=f.t_regdecvals_adj_merra2;
        q_regdecvals_merra2=f.q_regdecvals_merra2;
        tw_doyofregannmax_merra2=f.tw_doyofregannmax_merra2;
        yearofmax_tw_merra2=f.yearofmax_tw_merra2;yearofmax_tw_adj_merra2=f.yearofmax_tw_adj_merra2;
        yearofmax_t_merra2=f.yearofmax_t_merra2;yearofmax_t_adj_merra2=f.yearofmax_t_adj_merra2;
        yearofmax_q_merra2=f.yearofmax_q_merra2;yearofmax_q_adj_merra2=f.yearofmax_q_adj_merra2;
        doyofmax_tw_merra2=f.doyofmax_tw_merra2;doyofmax_tw_adj_merra2=f.doyofmax_tw_adj_merra2;
        doyofmax_t_merra2=f.doyofmax_t_merra2;doyofmax_t_adj_merra2=f.doyofmax_t_adj_merra2;
        doyofmax_q_merra2=f.doyofmax_q_merra2;doyofmax_q_adj_merra2=f.doyofmax_q_adj_merra2;
        maxstdanom_tw_merra2=f.maxstdanom_tw_merra2;maxstdanom_tw_adj_merra2=f.maxstdanom_tw_adj_merra2;
        maxstdanom_t_merra2=f.maxstdanom_t_merra2;maxstdanom_t_adj_merra2=f.maxstdanom_t_adj_merra2;
        maxstdanom_q_merra2=f.maxstdanom_q_merra2;maxstdanom_q_adj_merra2=f.maxstdanom_q_adj_merra2;
        peakmonth_tw_merra2=f.peakmonth_tw_merra2;peakmonth_t_merra2=f.peakmonth_t_merra2;peakmonth_q_merra2=f.peakmonth_q_merra2;
    f=load(strcat(saveloc,'mylatestjra55output_regions.mat'));
        tw_doyofregannmax_jra55=f.tw_doyofregannmax_jra55;
        tw_regdecvals_jra55=f.tw_regdecvals_jra55;t_regdecvals_jra55=f.t_regdecvals_jra55;q_regdecvals_jra55=f.q_regdecvals_jra55;
    f=load(strcat(saveloc,'jra55output_regions.mat'));
        tw_regdecstdanoms_adj_jra55=f.tw_regdecstdanoms_adj_jra55;
        t_regdecstdanoms_adj_jra55=f.t_regdecstdanoms_adj_jra55;
        q_regdecstdanoms_adj_jra55=f.q_regdecstdanoms_adj_jra55;
        tw_regdecvals_adj_jra55=f.tw_regdecvals_adj_jra55;
        t_regdecvals_adj_jra55=f.t_regdecvals_adj_jra55;
        yearofmax_tw_adj_jra55=f.yearofmax_tw_adj_jra55;
        yearofmax_t_adj_jra55=f.yearofmax_t_adj_jra55;
        yearofmax_q_adj_jra55=f.yearofmax_q_adj_jra55;
        doyofmax_tw_adj_jra55=f.doyofmax_tw_adj_jra55;
        doyofmax_t_adj_jra55=f.doyofmax_t_adj_jra55;
        doyofmax_q_adj_jra55=f.doyofmax_q_adj_jra55;
        maxstdanom_tw_adj_jra55=f.maxstdanom_tw_adj_jra55;
        maxstdanom_t_adj_jra55=f.maxstdanom_t_adj_jra55;
        maxstdanom_q_adj_jra55=f.maxstdanom_q_adj_jra55;
        peakmonth_tw_jra55=f.peakmonth_tw_jra55;peakmonth_t_jra55=f.peakmonth_t_jra55;peakmonth_q_jra55=f.peakmonth_q_jra55;
    f=load(strcat(saveloc,'regannmaxes'));
        tw_regannmax_adj_era5=f.tw_regannmax_adj_era5;t_regannmax_adj_era5=f.t_regannmax_adj_era5;
        tw_regannmax_era5=f.tw_regannmax_era5;t_regannmax_era5=f.t_regannmax_era5;
        tw_regannmax_adj_merra2=f.tw_regannmax_adj_merra2;t_regannmax_adj_merra2=f.t_regannmax_adj_merra2;
        tw_regannmax_merra2=f.tw_regannmax_merra2;t_regannmax_merra2=f.t_regannmax_merra2;
        tw_regannmax_adj_jra55=f.tw_regannmax_adj_jra55;t_regannmax_adj_jra55=f.t_regannmax_adj_jra55;
        tw_regannmax_jra55=f.tw_regannmax_jra55;t_regannmax_jra55=f.t_regannmax_jra55;

    f=load(strcat(saveloc,'longestduration_models.mat'));
        longestdurtoplot_canesm52=f.longestdurtoplot_canesm52;longestdurtoplot_canesm52_byensmem=f.longestdurtoplot_canesm52_byensmem;
        longestdurtoplot_miroc62=f.longestdurtoplot_miroc62;longestdurtoplot_miroc62_byensmem=f.longestdurtoplot_miroc62_byensmem;
        longestdurtoplot_mpige2=f.longestdurtoplot_mpige2;longestdurtoplot_mpige2_byensmem=f.longestdurtoplot_mpige2_byensmem;
    f=load(strcat(saveloc,'gmst_models.mat'));
        gmst_canesm5_difffrom1=f.gmst_canesm5_difffrom1;gmst_miroc6_difffrom1=f.gmst_miroc6_difffrom1;
        gmst_mpige_difffrom1=f.gmst_mpige_difffrom1;

    f=load(strcat(saveloc,'tw_valsandstdanoms_adj_canesm5'));
        tw_vals_adj_canesm5=f.tw_vals_adj_canesm5;tw_stdanoms_adj_canesm5=f.tw_stdanoms_adj_canesm5;
    f=load(strcat(saveloc,'tw_valsandstdanoms_adj_miroc6'));
        tw_vals_adj_miroc6=f.tw_vals_adj_miroc6;tw_stdanoms_adj_miroc6=f.tw_stdanoms_adj_miroc6;
    f=load(strcat(saveloc,'tw_valsandstdanoms_adj_mpige'));
        tw_vals_adj_mpige=f.tw_vals_adj_mpige;tw_stdanoms_adj_mpige=f.tw_stdanoms_adj_mpige;

    impossible_pct_canesm5=csvread(strcat(saveloc,'impossiblepct_canesm5.csv'));impossible_pct_3yr_canesm5=csvread(strcat(saveloc,'impossiblepct_3yr_canesm5.csv'));
    impossible_pct_miroc6=csvread(strcat(saveloc,'impossiblepct_miroc6.csv'));impossible_pct_3yr_miroc6=csvread(strcat(saveloc,'impossiblepct_3yr_miroc6.csv'));
    impossible_pct_mpige=csvread(strcat(saveloc,'impossiblepct_mpige.csv'));impossible_pct_3yr_mpige=csvread(strcat(saveloc,'impossiblepct_3yr_mpige.csv'));
    ret_per_full_canesm5=csvread(strcat(saveloc,'retperfull_canesm5.csv'));ret_per_full_3yr_canesm5=csvread(strcat(saveloc,'retperfull_3yr_canesm5.csv'));
    ret_per_full_miroc6=csvread(strcat(saveloc,'retperfull_miroc6.csv'));ret_per_full_3yr_miroc6=csvread(strcat(saveloc,'retperfull_3yr_miroc6.csv'));
    ret_per_full_mpige=csvread(strcat(saveloc,'retperfull_mpige.csv'));ret_per_full_3yr_mpige=csvread(strcat(saveloc,'retperfull_3yr_mpige.csv'));

    mask_regs_t=csvread(strcat(saveloc,'mask_regs_t_',compareagainst,'.csv'));
    mask_regs_tw=csvread(strcat(saveloc,'mask_regs_tw_',compareagainst,'.csv'));

    %Reload return periods (calculated in Python)
    currecretpers_tw=csvread(strcat(saveloc,'currentrecordreturnperiods_tw_era5.csv'));
    currecretpers_t=csvread(strcat(saveloc,'currentrecordreturnperiods_t_era5.csv'));
    currecretpers_withtopevents_tw=csvread(strcat(saveloc,'currentrecordreturnperiods_withtopevents_tw_era5.csv')); %includes top events
    currecretpers_withtopevents_t=csvread(strcat(saveloc,'currentrecordreturnperiods_withtopevents_t_era5.csv'));

    currecretpers_tw_jra55=csvread(strcat(saveloc,'currentrecordreturnperiods_tw_jra55.csv'));
    %currecretpers_t_jra55=csvread(strcat(saveloc,'currentrecordreturnperiods_t_jra55.csv'));
    currecretpers_withtopevents_tw_jra55=csvread(strcat(saveloc,'currentrecordreturnperiods_withtopevents_tw_jra55.csv')); %includes top events
    %currecretpers_withtopevents_t_jra55=csvread(strcat(saveloc,'currentrecordreturnperiods_withtopevents_t_jra55.csv'));


    twretpers_log=log10(currecretpers_tw);
    tretpers_log=log10(currecretpers_t);

    %Annual observed GMST from NASA
    gmst_obs=[0.06;0.03;0.05;-0.20;-0.11;-0.06;-0.02;-0.08;0.05;0.03;-0.08;0.01;0.16;-0.07;-0.01;
        -0.10;0.18;0.07;0.16;0.26;0.32;0.14;0.31;0.16;0.12;0.18;0.32;0.39;0.27;0.45;
        0.41; 0.22; 0.23; 0.31; 0.45; 0.33; 0.46; 0.61; 0.38; 0.39; 0.54; 0.63; 
        0.62; 0.53; 0.68; 0.64; 0.66; 0.54; 0.66; 0.72; 0.61; 0.65; 0.68; 0.74; 0.90; 
        1.01; 0.92; 0.85; 0.98; 1.01; 0.85; 0.89; 1.17]; %1961-2023, as also in Python script; computed versus 1951-1980 mean
        %Note that 1961-1990 mean is 0.1C

    %Define El Nino and La Nina events
    ensoindex=load(strcat(icloud,'General_Academics/Research/KeyFiles/indicesmonthlynino34.txt'))';
    ensoindex1d=reshape(ensoindex(2:end,:),[74*12 1]); %1950-2023
    ensoindex1d=ensoindex1d(133:888); %1961-2023
    %Get JFM-mean index for each year
    for y=1:numyr
        jfmmeanensoidx(y)=mean(ensoindex1d(y*12-11:y*12-9));
    end
end



%Use same regions as in Thompson et al. 2023
if defineregions==1
    regidxs=ncread(strcat(saveloc,'region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc'),'region');
    reglat1d=ncread(strcat(saveloc,'region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc'),'lat');
    reglon1d=ncread(strcat(saveloc,'region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc'),'lon');
    [lat2d_reg,lon2d_reg]=latlon_2dfrom1d(reglat1d,reglon1d);
    lat2d_reg=flipud(lat2d_reg');lon2d_reg=flipud(lon2d_reg');
    regidxs=flipud(regidxs');
    numreg=double(max(max(regidxs))+1);

    %Interpolate regions to match size of various datasets
    regs_era5sz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_era5,lat2d_era5)+1; %so ordinates start at 1 rather than 0
    regs_merra2sz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_merra2,lat2d_merra2)+1; %so ordinates start at 1 rather than 0
    regs_jra55sz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_jra55,lat2d_jra55)+1; %so ordinates start at 1 rather than 0
    regs_canesm5sz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_canesm5,lat2d_canesm5)+1; %so ordinates start at 1 rather than 0
    regs_miroc6sz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_miroc6,lat2d_miroc6)+1; %so ordinates start at 1 rather than 0
    regs_mpigesz=interp2(lon2d_reg,lat2d_reg,regidxs,lon2d_mpige,lat2d_mpige)+1; %so ordinates start at 1 rather than 0

    clear lat2d_reg;clear lon2d_reg;clear regidxs;

    %To resolve the handful of odd cases where a point straddles a region
        %boundary after interpolation, simply assign it to the region to the north or east
    inputs={regs_era5sz;regs_merra2sz;regs_canesm5sz;regs_miroc6sz;regs_mpigesz};
    for loop=1:5
        input=inputs{loop};

        arrnew=NaN.*ones(size(input));
        tofix=rem(input,1)~=0;
        for i=2:size(tofix,1)
            for j=1:size(tofix,2)-1
                if tofix(i,j)==1
                    if ~isnan(input(i-1,j)) && rem(input(i-1,j),1)==0 %check if can replace with point to north
                        arrnew(i,j)=input(i-1,j);
                    elseif ~isnan(input(i,j+1)) && rem(input(i,j+1),1)==0 %if needed, check if can replace with point to east
                        arrnew(i,j)=input(i,j+1);
                    else %point is in ocean, or maybe on land but nothing will be done
                    end
                end
            end
        end
        for i=1:size(tofix,1)
            for j=1:size(tofix,2)
                if tofix(i,j)~=1
                    arrnew(i,j)=input(i,j);
                end
            end
        end
        if loop==1
            regs_era5sz=arrnew;
        elseif loop==2
            regs_merra2sz=arrnew;
        elseif loop==3
            regs_canesm5sz=arrnew;
        elseif loop==4
            regs_miroc6sz=arrnew;
        elseif loop==5
            regs_mpigesz=arrnew;
        end
    end

    %Also, for era5sz, fix a handful of artifacts at region boundaries,
    %where point straddles region but result is an integer so the above
    %loop doesn't pick it up
    regs_era5sz(98:119,157)=1;
    regs_era5sz(436:446,1237)=167;
    regs_era5sz(466:488,1237)=164;
    regs_era5sz(466:468,1285)=170;


    %Reminder: to quickly see where a region is located:
    %imsq(regs_era5sz==130)

    %Get regions' center points (useful for identifying them later on,
    %including in extension Python script)
    clear regcenterlats;clear regcenterlons;
    for reg=1:double(max(max(regs_era5sz)))
        regplotted=regs_era5sz==reg;

        rowsincl=[];colsincl=[];
        for i=1:latsz_era5
            for j=1:lonsz_era5
                if regplotted(i,j)==1
                    rowsincl=[rowsincl;i];
                    colsincl=[colsincl;j];
                end
            end
        end
        centerrow=round(median(sort(rowsincl)));centercol=round(median(sort(colsincl)));

        regcenterlats(reg)=double(lat2d_era5(centerrow,centercol));
        regcenterlons(reg)=double(lon2d_era5(centerrow,centercol));
    end
    A=[regcenterlats' regcenterlons'];
    writematrix(A,strcat(saveloc,'region_centerpoints.csv'));
end


%Process CanESM5
if processcanesm5==1
    z1=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_canesm5,lat2d_canesm5);
    psfc_tmp=pressurefromelev(z1);
    %psfc_canesm5=repmat(psfc_tmp,[1 1 365*numyr]);

    if calcmaxs_canesm5==1
        tw_annmax_canesm5=NaN.*ones(canesm5_enssz,64,128,numyr);

        psfc_canesm5=repmat(psfc_tmp,[1 1 365*numyrtodo]);

        for phys=1:2
            for i=1:canesm5_enssz/2
                idx=phys*25-25+i; %ens member
    
                tasmax=NaN.*ones(64,128,numyrtodo*365);hussmean=NaN.*ones(64,128,numyrtodo*365);
    
                %Historical
                %day 51101 -- 1990/01/01
                %day 40516 -- 1961/01/01
                if dohist==1
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'tasmax',[1 1 40516],[Inf Inf Inf])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
    
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'huss',[1 1 40516],[Inf Inf Inf]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
                end
        
    
                %SSP5-8.5
                %Read only 2015-2023
                if dofut==1
                %tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                %    'tasmax',[1 1 1],[Inf Inf 365*ny_fut])-273.15; %C
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'tasmax',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
        
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'huss',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
                end


                %Calculate wet-bulb temperature using daily-max T and daily-mean q
                tw_canesm5_thisensmem=calcwbt_daviesjones(tasmax,psfc_canesm5.*100,hussmean);
        
                %Annual-max Tw
                firstyrhere=ny_hist+futfirstyr; %default is 1
                for y=1:numyr
                %for y=ny_hist+futfirstyr:ny_hist+futlastyr
                    thisrelyr=y-firstyrhere+1;
                    ystart=365*thisrelyr-364;yend=365*thisrelyr;
                    tw_annmax_canesm5(idx,:,:,y)=max(tw_canesm5_thisensmem(:,:,ystart:yend),[],3);
                end
    
                clear tasmax;clear hussmean;
                disp(phys);disp(i);disp(clock);
                save(strcat(saveloc,'canesm5output.mat'),'tw_annmax_canesm5','-v7.3');
            end
        end
        save(strcat(saveloc,'canesm5output.mat'),'tw_annmax_canesm5','-v7.3');
    end

    %To update and/or combine with the above!!
    %Get annual mean Tw and T for each year
    if calcmeans_canesm5==1
        tw_annmean_canesm5=NaN.*ones(10,64,128,numyr);t_annmean_canesm5=NaN.*ones(10,64,128,numyr);
        psfc_canesm5=repmat(psfc_tmp,[1 1 365*numyr]);
        for phys=1:2
            for i=1:25
                %Part A: about 40 sec%
                idx=phys*25-25+i; %ens member
    
                tasmax=NaN.*ones(64,128,numyr*365);hussmean=NaN.*ones(64,128,numyr*365);
    
                %Historical
                %day 51101 -- 1990/01/01
                %day 40516 -- 1961/01/01
                if dohist==1
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'tasmax',[1 1 40516],[Inf Inf Inf])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
    
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'huss',[1 1 40516],[Inf Inf Inf]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
                end
        
    
                %SSP5-8.5
                %Read only 2015-2023
                if dofut==1
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'tasmax',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
        
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'huss',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
                end
    
    
                %Calculate wet-bulb temperature using daily-max T and daily-mean q
                %DO PYTHON VERSION, AS MATLAB VERSION IS CRASH-INDUCING
                %tw_canesm5_thisensmem=calcwbt_daviesjones(tasmax,psfc_canesm5.*100,hussmean);
                arrt=tasmax;arrp=psfc_canesm5.*100;arrq=hussmean.*1000;suffixforfilename=strcat('_ensmem',num2str(idx));
                calctw_matlabpythonhelper_part1;

                %End of Part A%


                %%%%%
                %(STOP *Run Python code as instructed in calctw_matlabpythonhelper_part1* STOP)
                %Only thing that should need changing there is suffixforfilename
                %about 4 min
                %%%%%

                
                %Part B: about 30 sec%
                myfname='tw_canesm5';
                calctw_matlabpythonhelper_part2;

                tw_canesm5_thisensmem=tw2m;clear tw2m;
        
                %Annual-mean Tw
                firstyrhere=ny_hist+futfirstyr; %default is 1
                %for y=1:numyr
                for y=ny_hist+futfirstyr:ny_hist+futlastyr
                    thisrelyr=y-firstyrhere+1;
                    ystart=365*thisrelyr-364;yend=365*thisrelyr;
                    tw_annmean_canesm5(idx,:,:,y)=mean(tw_canesm5_thisensmem(:,:,ystart:yend),3);
                    t_annmean_canesm5(idx,:,:,y)=mean(tasmax(:,:,ystart:yend),3);
                end
    
                disp(phys);disp(i);disp(clock);
                save(strcat(saveloc,'canesm5_globalannualmeans.mat'),'tw_annmean_canesm5','t_annmean_canesm5','-v7.3');
    
                clear tasmax;clear hussmean;
                %End of Part B%
            end
            clear tw_canesm5_thisensmem;
        end
    end
    clear psfc_canesm5;
    disp('Finished processing CanESM5!');disp(clock);
end


if processmiroc6==1
    %Historical set-up
    decades={'19600101-19691231';'19700101-19791231';'19800101-19891231';'19900101-19991231';'20000101-20091231';'20100101-20141231'};
    arrayindices_start=[1;365*9+1;365*19+1;365*29+1;365*39+1;365*49+1];arrayindices_stop=[365*9;365*19;365*29;365*39;365*49;365*54]; %where data will go
    fileindices_start=[366;1;1;1;1;1];fileindices_stop=[365*10;365*10;365*10;365*10;365*10;365*5]; %where data comes from

    z2=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_miroc6,lat2d_miroc6);
    psfc_tmp=pressurefromelev(z2);
    psfc_miroc6=repmat(psfc_tmp,[1 1 365*numyrtodo]);

    if calcmaxs_miroc6==1
        tw_annmax_miroc6=NaN.*ones(miroc6_enssz,128,256,numyr);
    
        phys=1;
        for i=1:miroc6_enssz %1:10
            for latloop=1:4 %necessary to avoid crashing Matlab   
                tasmax=NaN.*ones(32,256,numyrtodo*365);hussmean=NaN.*ones(32,256,numyrtodo*365);
    
                if latloop==1
                    l1=1;l2=32;
                elseif latloop==2
                    l1=33;l2=64;
                elseif latloop==3
                    l1=65;l2=96;
                elseif latloop==4
                    l1=97;l2=128;
                end
        
                %Historical
                if dohist==1
                for d=1:length(decades)
                    thisdec=decades{d};
                    arraystart=arrayindices_start(d);arraystop=arrayindices_stop(d);
                    filestart=fileindices_start(d);filestop=fileindices_stop(d);
        
                    tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'tasmax')-273.15; %C
                    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1)/4,size(tmp2,2),size(tmp2,3));
                    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                    tasmax(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
        
                    tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'huss'); %kg/kg
                    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1)/4,size(tmp2,2),size(tmp2,3));
                    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                    hussmean(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
                end
                end
        
        
                %SSP5-8.5
                %Read only 2015-2023
                if dofut==1
                tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'tasmax')-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1)/4,size(tmp2,2),size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                tasmax(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp3;
        
                tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'huss'); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1)/4,size(tmp2,2),size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                hussmean(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp3;
                end
                
                %Calculate wet-bulb temperature using daily-max T and daily-mean q           
                tw_miroc6_thisensmem=calcwbt_daviesjones(tasmax,psfc_miroc6(l1:l2,:,:).*100,hussmean);
                clear tasmax;clear hussmean;
        
                %Annual-max Tw
                firstyrhere=ny_hist+futfirstyr; %default is 1
                %for y=1:numyr
                for y=ny_hist+futfirstyr:ny_hist+futlastyr
                    thisrelyr=y-firstyrhere+1;
                    ystart=365*thisrelyr-364;yend=365*thisrelyr;
                    tw_annmax_miroc6(i,l1:l2,:,y)=max(tw_miroc6_thisensmem(:,:,ystart:yend),[],3);
                end
        
                clear tasmax;clear hussmean;
                save(strcat(saveloc,'miroc6output.mat'),'tw_annmax_miroc6','-v7.3');
            end
            disp(i);disp(clock);clear tw_miroc6_thisensmem;
        end
        save(strcat(saveloc,'miroc6output.mat'),'tw_annmax_miroc6','-v7.3');
    end

    %Get ens-mean global-mean-surface Tw and T for MIROC6 for each year
    if calcmeans_miroc6==1
        tw_annmean_miroc6=NaN.*ones(5,128,256,numyr);t_annmean_miroc6=NaN.*ones(5,128,256,numyr);
        phys=1;
        for i=1:5
            %Part A: about 4 min%
            tasmax=NaN.*ones(128,256,numyr*365);hussmean=NaN.*ones(128,256,numyr*365);
            l1=1;l2=128;
    
            %Historical
            for d=1:length(decades)
                thisdec=decades{d};
                arraystart=arrayindices_start(d);arraystop=arrayindices_stop(d);
                filestart=fileindices_start(d);filestop=fileindices_stop(d);
    
                tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'tasmax')-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1),size(tmp2,2),size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                tasmax(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
    
                tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'huss'); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1),size(tmp2,2),size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                hussmean(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
            end
    
    
            %SSP5-8.5
            %Read only 2015-2023
            tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'tasmax')-273.15; %C
            tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1),size(tmp2,2),size(tmp2,3));
            for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
            tasmax(:,:,ny_hist*365+1:end)=tmp3(:,:,1:365*ny_fut);clear tmp3;
    
            tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'huss'); %kg/kg
            tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2,1),size(tmp2,2),size(tmp2,3));
            for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
            hussmean(:,:,ny_hist*365+1:end)=tmp3(:,:,1:365*ny_fut);clear tmp3;
    
            
            %Calculate wet-bulb temperature using daily-max T and daily-mean q           
            %tw_miroc6_thisensmem=calcwbt_daviesjones(tasmax,psfc_miroc6(l1:l2,:,:).*100,hussmean);

            arrt=tasmax;arrp=psfc_miroc6.*100;arrq=hussmean.*1000;suffixforfilename=strcat('_ensmem',num2str(i));
            calctw_matlabpythonhelper_part1;
            
            clear hussmean;

            %End of Part A%


            %%%%%
            %(STOP *Run Python code as instructed in calctw_matlabpythonhelper_part1* STOP)
            %Between ensemble members, the only thing that should need changing there is suffixforfilename
            %about 17 min
            %%%%%

    
            %Annual-mean Tw and T
            %Part B: about 6 min%
            myfname='tw_miroc6';
            calctw_matlabpythonhelper_part2;

            tw_miroc6_thisensmem=tw2m;clear tw2m;

            for y=1:numyr
                ystart=365*y-364;yend=365*y;
                tw_annmean_miroc6(i,l1:l2,:,y)=mean(tw_miroc6_thisensmem(:,:,ystart:yend),3);
                t_annmean_miroc6(i,l1:l2,:,y)=mean(tasmax(:,:,ystart:yend),3);
            end
    
            clear tasmax;clear hussmean;
            clear tw_miroc6_thisensmem;
            save(strcat(saveloc,'miroc6_globalannualmeans.mat'),'tw_annmean_miroc6','t_annmean_miroc6','-v7.3');
            disp(i);disp(clock);
            %End of Part B%
        end
    end
    clear psfc_miroc6;
    disp('Finished processing MIROC6!');disp(clock);
end


%Process MPI-GE
if processmpige==1
    z3=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_mpige,lat2d_mpige);
    psfc_tmp=pressurefromelev(z3);
    psfc_mpige=repmat(psfc_tmp,[1 1 365*numyrtodo]);
    mlatsz=size(lat2d_mpige,1);mlonsz=size(lat2d_mpige,2);

    if calcmaxs_mpige==1
        tw_annmax_mpige=NaN.*ones(mpige_enssz,96,192,numyr);
        for i=1:mpige_enssz
            tasmax=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);hussmean=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);
    
            %Historical
            %Read 1961-2014
            if dohist==1
                filepath_tas_part1=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');
                filepath_huss_part1=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');

                decades={'19500101-19691231';'19700101-19891231';'19900101-20091231';'20100101-20141231'};
                arrayindices_start=[1;365*9+1;365*29+1;365*49+1];arrayindices_stop=[365*9;365*29;365*49;365*54]; %where data will go
                fileindices_start=[365*11+1;1;1;1];fileindices_stop=[365*20;365*20;365*20;365*5]; %where data comes from

                histpart=1;futpart=0;
                readmodeldata_helper;
            end
    
    
            %SSP5-8.5
            %Read only 2015-2023 for now
            if dofut==1
                filepath_tas=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');
                filepath_huss=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');

                histpart=0;futpart=1;
                readmodeldata_helper;
            end
    
            %Calculate Tw (using daily-max T and daily-mean q) and get annual max (calc time 3 min per ensemble member)
            tw_mpige_thisensmem=NaN.*ones(size(tasmax));
            firstyrhere=1;
            %firstyrhere=ny_hist+futfirstyr; %default is 1
            for y=1:numyr
            %for y=ny_hist+futfirstyr:ny_hist+futlastyr
                thisrelyr=y-firstyrhere+1;
                ystart=365*thisrelyr-364;yend=365*thisrelyr;
                tw_mpige_thisensmem=calcwbt_daviesjones(tasmax(:,:,ystart:yend),psfc_mpige(:,:,ystart:yend).*100,hussmean(:,:,ystart:yend));
                tw_annmax_mpige(i,:,:,y)=max(tw_mpige_thisensmem,[],3);
            end
            clear tasmax;clear hussmean;
            disp(i);disp(clock);clear tw_mpige_thisensmem;

            %save(strcat(saveloc,'mpigeoutput.mat'),'tw_annmax_mpige','-v7.3');
        end
    end

    if calcmeans_mpige==1
        tw_annmean_mpige=NaN.*ones(6,mlatsz,mlonsz,numyr);t_annmean_mpige=NaN.*ones(6,mlatsz,mlonsz,numyr);
        for i=1:6
            %Part A: about 2 min%
            tasmax=NaN.*ones(mlatsz,mlonsz,numyr*365);hussmean=NaN.*ones(mlatsz,mlonsz,numyr*365);
    
            %Historical
            %Read 1961-2014
            if dohist==1
                filepath_tas_part1=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');
                filepath_huss_part1=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');

                decades={'19500101-19691231';'19700101-19891231';'19900101-20091231';'20100101-20141231'};
                arrayindices_start=[1;365*9+1;365*29+1;365*49+1];arrayindices_stop=[365*9;365*29;365*49;365*54]; %where data will go
                fileindices_start=[365*11+1;1;1;1];fileindices_stop=[365*20;365*20;365*20;365*5]; %where data comes from

                histpart=1;futpart=0;
                readmodeldata_helper;
            end
    
    
            %SSP5-8.5
            %Read only 2015-2023 for now
            if dofut==1
                filepath_tas=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');
                filepath_huss=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');

                histpart=0;futpart=1;
                readmodeldata_helper;
            end

            arrt=tasmax;arrp=psfc_mpige.*100;arrq=hussmean.*1000;suffixforfilename=strcat('_ensmem',num2str(i));
            calctw_matlabpythonhelper_part1;
            
            clear hussmean;

            %End of Part A%


            %%%%%
            %(STOP *Run Python code as instructed in calctw_matlabpythonhelper_part1* STOP)
            %Between ensemble members, the only thing that should need changing there is suffixforfilename
            %about 10 min
            %%%%%


            %Part B: about 2 min%
            myfname='tw_mpige';
            calctw_matlabpythonhelper_part2;

            tw_mpige_thisensmem=tw2m;clear tw2m;
    
            %Annual-mean Tw
            for y=1:numyr
                ystart=365*y-364;yend=365*y;
                %tw_mpige_thisensmem=calcwbt_daviesjones(tasmax(:,:,ystart:yend),psfc_mpige(:,:,ystart:yend).*100,hussmean(:,:,ystart:yend));
                
                tw_annmean_mpige(i,:,:,y)=mean(tw_mpige_thisensmem(:,:,ystart:yend),3);
                t_annmean_mpige(i,:,:,y)=mean(tasmax(:,:,ystart:yend),3);
            end
            disp(i);disp(clock);clear tasmax;clear hussmean;clear tw_mpige_thisensmem;
            save(strcat(saveloc,'mpige_globalannualmeans.mat'),'tw_annmean_mpige','t_annmean_mpige','-v7.3');
            %End of Part B%
        end
    end
    clear psfc_mpige;
    disp('Finished processing MPI-GE!');disp(clock);
end


if convertmodeldatatoregsandnormalize==1
    for model=1:3 %CanESM5, MIROC6, MPI-ESM
        if model==1
            enssz=canesm5_enssz;tw_regannmax_canesm5=NaN.*ones(enssz,numreg,numyr);lat2d=lat2d_canesm5;lon2d=lon2d_canesm5;regsz=regs_canesm5sz;
            twannmax=tw_annmax_canesm5;modelname='canesm5';
        elseif model==2
            enssz=miroc6_enssz;tw_regannmax_miroc6=NaN.*ones(enssz,numreg,numyr);lat2d=lat2d_miroc6;lon2d=lon2d_miroc6;regsz=regs_miroc6sz;
            twannmax=tw_annmax_miroc6;modelname='miroc6';
        elseif model==3
            enssz=mpige_enssz;tw_regannmax_mpige=NaN.*ones(enssz,numreg,numyr);lat2d=lat2d_mpige;lon2d=lon2d_mpige;regsz=regs_mpigesz;
            twannmax=tw_annmax_mpige;modelname='mpige';
        end

        %Organize the data we're actually interested in
        tw_regannmaxs=cell(numreg);
        for ensmem=1:enssz
            c=zeros(numreg,1);
            for i=1:size(lat2d,1)
                for j=1:size(lat2d,2)
                    thisreg=round(regsz(i,j));
                    if ~isnan(thisreg)
                        c(thisreg)=c(thisreg)+1;
                        tw_regannmaxs{thisreg}(ensmem,c(thisreg),:)=twannmax(ensmem,i,j,:); %last dim is years
                    end
                end
            end
        end
        clear twannmax;
        %Spatial mean for each region, by ensemble member and year
        if model==1
            for r=1:numreg;tw_regannmax_canesm5(:,r,:)=squeeze(mean(squeeze(tw_regannmaxs{r}),2));end
            writematrix(round2(permute(tw_regannmax_canesm5,[1 3 2]),0.001),strcat(saveloc,'canesm5output_regions.txt'));
                %this ordering of dims is more conducive for analyzing in Python
        elseif model==2
            for r=1:numreg;tw_regannmax_miroc6(:,r,:)=squeeze(mean(squeeze(tw_regannmaxs{r}),2));end
            writematrix(round2(permute(tw_regannmax_miroc6,[1 3 2]),0.001),strcat(saveloc,'miroc6output_regions.txt'));
        elseif model==3
            for r=1:numreg;tw_regannmax_mpige(:,r,:)=squeeze(mean(squeeze(tw_regannmaxs{r}),2));end
            writematrix(round2(permute(tw_regannmax_mpige,[1 3 2]),0.001),strcat(saveloc,'mpigeoutput_regions.txt'));
        end
        clear tw_regannmaxs;
    
        %Adjust each model by its ens-mean GMST (lat-weighted) as in Thompson et al. 2023
        f=load(strcat(saveloc,modelname,'_globalannualmeans.mat'));
        if model==1;t_annmean=f.t_annmean_canesm5;elseif model==2;t_annmean=f.t_annmean_miroc6;elseif model==3;t_annmean=f.t_annmean_mpige;end
        invalid=t_annmean==0;t_annmean(invalid)=NaN;
        latweight=cos(deg2rad(lat2d));latweight=latweight(:,1);
        latbandmeanbyyear=squeeze(mean(squeeze(mean(t_annmean,'omitnan')),2));clear t_annmean;
        clear gmst_model;for y=1:numyr;gmst_model(y)=mean(latweight.*latbandmeanbyyear(:,y));end


        %Normalize model timeseries
        %The Python code already has functions to do this, but doing it
            %here more easily allows for visualization and verification that
            %the code is doing what we expect it to
        %Mean of 1961-1990 should be same as obs 1961-1990 mean, i.e. ~0.1C
        %Furthermore, THIS IS THE CORE OF WHAT IS USED IN REPLICATING T23 FIG 4, adjust
            %reg-ann-max Tw timeseries by corresp year's ens-mean GMST, so
            %that all results are for a GMST of 1.0C
        gmst19611990mean=mean(gmst_model(1:30));
        gmst_model_norm=gmst_model-gmst19611990mean+mean(gmst_obs(1:30));
        gmst_model_difffrom1=1-gmst_model_norm;
    
        if model==1
            %gmst_canesm5=gmst_norm;%writematrix(gmst,strcat(saveloc,'canesm5gmst.txt'));
            tw_regannmax_canesm5_adj=tw_regannmax_canesm5+reshape(repmat(repmat(gmst_model_difffrom1,[numreg 1]),[canesm5_enssz 1 1]),[canesm5_enssz numreg numyr]);
            tw_regoverallmax_canesm5_adj=max(tw_regannmax_canesm5_adj,[],3);
            writematrix(round2(tw_regoverallmax_canesm5_adj,0.001),strcat(saveloc,'tw_regoverallmax_canesm5_adj.txt'));

            gmst_canesm5_difffrom1=gmst_model_difffrom1;
            save(strcat(saveloc,'gmst_models.mat'),'gmst_canesm5_difffrom1','-append');
        elseif model==2
            tw_regannmax_miroc6_adj=tw_regannmax_miroc6+reshape(repmat(repmat(gmst_model_difffrom1,[numreg 1]),[miroc6_enssz 1 1]),[miroc6_enssz numreg numyr]);
            tw_regoverallmax_miroc6_adj=max(tw_regannmax_miroc6_adj,[],3);
            writematrix(round2(tw_regoverallmax_miroc6_adj,0.001),strcat(saveloc,'tw_regoverallmax_miroc6_adj.txt'));

            gmst_miroc6_difffrom1=gmst_model_difffrom1;
            save(strcat(saveloc,'gmst_models.mat'),'gmst_miroc6_difffrom1','-append');
        elseif model==3
            tw_regannmax_mpige_adj=tw_regannmax_mpige+reshape(repmat(repmat(gmst_model_difffrom1,[numreg 1]),[mpige_enssz 1 1]),[mpige_enssz numreg numyr]);
            tw_regoverallmax_mpige_adj=max(tw_regannmax_mpige_adj,[],3);
            writematrix(round2(tw_regoverallmax_mpige_adj,0.001),strcat(saveloc,'tw_regoverallmax_mpige_adj.txt'));

            gmst_mpige_difffrom1=gmst_model_difffrom1;
            save(strcat(saveloc,'gmst_models.mat'),'gmst_mpige_difffrom1','-append');
        end


        %In Python script:
        %-fit a GEV to tw_regannmax_era5_adj, separately for each region
        %-compare this GEV to extrema of tw_regannmax_MODEL_adj
    end
    gmst_obs_difffrom1=1-gmst_obs;
    %tw_regannmax_era5_adj=tw_regannmax_era5+repmat(gmst_obs_difffrom1,[1 numreg])';
    %writematrix(round2(tw_regannmax_era5_adj,0.001),strcat(saveloc,'tw_regannmax_era5_adj.txt'));
end

%Partially reread CanESM5 to save daily data
%Also compute some stats for later use
if rereadcanesm5_savedaily==1
    if canesm5_dopart1==1
        z1=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_canesm5,lat2d_canesm5);
        psfc_tmp=pressurefromelev(z1);
        
        mlatsz=size(lat2d_canesm5,1);mlonsz=size(lat2d_canesm5,2);regsz=regs_canesm5sz;
        tw_regdaily_adj_canesm5=NaN.*ones(canesm5_enssz,numyrtodo,numreg,365);
    
        for phys=1:2
            for idx=1:canesm5_enssz/2
                i=phys*25-25+idx; %ens member
        
                tasmax=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);hussmean=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);
                psfc_canesm5=repmat(psfc_tmp,[1 1 numyrtodo*365]);
        
                %Historical
                %day 51101 -- 1990/01/01
                %day 40516 -- 1961/01/01
                if dohist==1
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_historical_r',num2str(idx),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'tasmax',[1 1 40516],[Inf Inf Inf])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
        
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_historical_r',num2str(idx),'i1p',num2str(phys),'f1_gn_18500101-20141231.nc'),...
                    'huss',[1 1 40516],[Inf Inf Inf]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,1:ny_hist*365)=tmp3;clear tmp2;clear tmp3;
                end
        
        
                %SSP5-8.5
                %Read only 2015-2023
                if dofut==1
                tmp1=ncread(strcat(canesm5dir,'tasmax_day_CanESM5_ssp585_r',num2str(idx),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'tasmax',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo])-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                tasmax(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
        
                tmp1=ncread(strcat(canesm5dir,'huss_day_CanESM5_ssp585_r',num2str(idx),'i1p',num2str(phys),'f1_gn_20150101-21001231.nc'),...
                    'huss',[1 1 365*(futfirstyr-1)+1],[Inf Inf 365*ny_fut_todo]); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
                hussmean(:,:,ny_hist_todo*365+1:end)=tmp3;clear tmp2;clear tmp3;
                end
    
                %Calculate wet-bulb temperature using daily-max T and daily-mean q           
                arrt=tasmax;arrp=psfc_canesm5.*100;arrq=hussmean.*1000;
                suffixforfilename=strcat('_ensmem',num2str(i));clear hussmean;
                calctw_matlabpythonhelper_part1;
            end
        end
        clear psfc_canesm5;
    end

    %%%%%
    %%Now, in a Jupyter notebook, run final cell of speedywetbulb_copy.ipynb [which is in home directory for ease of access]
        %This reads the saved values, computes DJ wet-bulb temperature, and saves to another netcdf file
    %Between ensemble members, the only thing that should need changing there is suffixforfilename
    %%%%%

    %Save daily Tw and split into regions
    if canesm5_dopart2==1
        myfname='tw_canesm5_daily'; %prefix to be used for saving final output
        for i=1:canesm5_enssz
            suffixes={strcat('_ensmem',num2str(i))}; %files to be processed for this ens mem

            %tw_canesm5_thisensmem=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);
            calctw_matlabpythonhelper_part2;

            tw_canesm5_thisensmem=load(strcat(saveloc,myfname,suffixes{1},'.mat')).tw2m;


            for y=1:numyr
                %Separate daily data into regions
                tw_regdata=cell(numreg);
                c=zeros(numreg,1);
                for latc=1:mlatsz
                    for lonc=1:mlonsz
                        thisreg=round(regsz(latc,lonc));
                        if ~isnan(thisreg)
                            c(thisreg)=c(thisreg)+1;
                            tw_regdata{thisreg}(c(thisreg),:)=tw_canesm5_thisensmem(latc,lonc,y*365-364:y*365); %lat | lon | days of year
                        end
                    end
                end
                %Spatial mean for each region and day
                for r=1:numreg;tw_regdailythisyear_canesm5(r,:)=squeeze(mean(squeeze(tw_regdata{r})));end
    
                %Adjust to a GMST of +1C
                adjfactor_thisyr=gmst_canesm5_difffrom1(y);
                tw_regdaily_adj_canesm5(i,y,:,:)=tw_regdailythisyear_canesm5+repmat(adjfactor_thisyr,[numreg 365]);
                clear tw_regdailythisyear_canesm5;
            end

            disp(i);disp(clock);clear tw_canesm5_thisensmem;
        end

        %Convert to std anoms
        tw_vals_adj_canesm5=tw_regdaily_adj_canesm5;clear tw_regdaily_adj_canesm5;
        for i=1:canesm5_enssz
            for reg=1:numreg
                for m=1:12
                    monthmeantw(reg,m)=mean(reshape(squeeze(tw_vals_adj_canesm5(i,firstclimoyr:lastclimoyr,reg,monthstarts(m):monthstops(m))),[monthlens(m)*numclimoyr 1]),'omitnan');
                end
                [~,peakmonth_tw_canesm5(i,reg)]=max(monthmeantw(reg,:));
            end
    
            %Calculate mean warm-season daily-max temperature, wet-bulb temperature, and specific humidity
            peakmonth_tw=squeeze(peakmonth_tw_canesm5(i,:));tw_vals=permute(squeeze(tw_vals_adj_canesm5(i,:,:,:)),[2 3 1]);
        
            varstocalc={'tw'};
            firstclimoyear=firstclimoyr;lastclimoyear=lastclimoyr;numclimoyears=numclimoyr;
            meanwarmseasondailymaxhelper;
        
            meanwsmax_tw_canesm5(i,:)=meanwsmax_tw;stdwsmax_tw_canesm5(i,:)=stdwsmax_tw;
        end
       
        %Calculate daily extreme index for all years, for Tw, using Thompson et al. method
        clear tw_stdanoms_adj_canesm5;
        for i=1:canesm5_enssz
            for reg=1:numreg
                tw_stdanoms_adj_canesm5(i,:,reg,:)=(squeeze(tw_vals_adj_canesm5(i,:,reg,:))'-squeeze(repmat(meanwsmax_tw_canesm5(i,reg),[1 365 numyr])))./...
                    squeeze(repmat(stdwsmax_tw_canesm5(i,reg),[1 365 numyr]));
            end
        end
        save(strcat(saveloc,'tw_valsandstdanoms_adj_canesm5'),'tw_vals_adj_canesm5','tw_stdanoms_adj_canesm5');

        %Calc longest duration above 1, 2, and 3 std anoms
        for loop=1:3
            stdanomcutoff_here=loop;
    
            longestdurtoplot=zeros(canesm5_enssz,mlatsz,mlonsz); %just to initialize
            for i=1:canesm5_enssz
                tmp=0.*regs_canesm5sz;
                for reg=1:numreg
                    longestdursofar=0;thisdur=0;
                    for y=1:numyr
                        for doy=1:365
                            if tw_stdanoms_adj_canesm5(i,doy,reg,y)>=stdanomcutoff_here %extr heat starts or continues
                                thisdur=thisdur+1;
                            else
                                if thisdur>longestdursofar
                                    longestdursofar=thisdur;
                                    longestdur_date(reg,1)=y+firstyear-1;longestdur_date(reg,2)=doy-1; %end date of sequence
                                end
                                thisdur=0;
                            end
                        end
                    end
                    longestdurbyreg(reg)=longestdursofar;
            
                    tmp(regs_canesm5sz==reg)=longestdurbyreg(reg);
                end
                longestdurtoplot(i,:,:)=tmp;
            end
            %Ensemble-mean thereof
            if loop==1
                longestdurtoplot_canesm51=squeeze(mean(longestdurtoplot));longestdurtoplot_canesm51_byensmem=longestdurtoplot;
            elseif loop==2
                longestdurtoplot_canesm52=squeeze(mean(longestdurtoplot));longestdurtoplot_canesm52_byensmem=longestdurtoplot;
            elseif loop==3
                longestdurtoplot_canesm53=squeeze(mean(longestdurtoplot));longestdurtoplot_canesm53_byensmem=longestdurtoplot;
            end
        end
    
        save(strcat(saveloc,'longestduration_models.mat'),...
            'longestdurtoplot_canesm51','longestdurtoplot_canesm52','longestdurtoplot_canesm53',...
            'longestdurtoplot_canesm51_byensmem','longestdurtoplot_canesm52_byensmem','longestdurtoplot_canesm53_byensmem','-append');
    end
end

%Partially reread MIROC6 to save daily data
%Also compute some stats for later use
if rereadmiroc6_savedaily==1
    if miroc6_dopart1==1
        decades={'19600101-19691231';'19700101-19791231';'19800101-19891231';'19900101-19991231';'20000101-20091231';'20100101-20141231'};
        arrayindices_start=[1;365*9+1;365*19+1;365*29+1;365*39+1;365*49+1];arrayindices_stop=[365*9;365*19;365*29;365*39;365*49;365*54]; %where data will go
        fileindices_start=[366;1;1;1;1;1];fileindices_stop=[365*10;365*10;365*10;365*10;365*10;365*5]; %where data comes from
    
        z2=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_miroc6,lat2d_miroc6);
        psfc_tmp=pressurefromelev(z2);
        psfc_miroc6=repmat(psfc_tmp,[1 1 365*numyrtodo]);
        mlatsz=size(lat2d_miroc6,1);mlonsz=size(lat2d_miroc6,2);regsz=regs_miroc6sz;
    
        tw_regdaily_adj_miroc6=NaN.*ones(miroc6_enssz,numyr,numreg,365);
    
        for i=1:miroc6_enssz
            tasmax=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);hussmean=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);
    
            phys=1;
            l1s=[1;33;65;97];l2s=[32;64;96;128];
            for latloop=1:4 %necessary to avoid crashing Matlab   
                tasmax=NaN.*ones(mlatsz/4,mlonsz,numyrtodo*365);hussmean=NaN.*ones(mlatsz/4,mlonsz,numyrtodo*365);
    
                l1=l1s(latloop);l2=l2s(latloop);
        
                %Historical
                if dohist==1
                for d=1:length(decades)
                    thisdec=decades{d};
                    arraystart=arrayindices_start(d);arraystop=arrayindices_stop(d);
                    filestart=fileindices_start(d);filestop=fileindices_stop(d);
        
                    tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'tasmax')-273.15; %C
                    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(mlatsz/4,mlonsz,size(tmp2,3));
                    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                    tasmax(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
        
                    tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_historical_r',num2str(i),'i1p',num2str(phys),'f1_gn_',thisdec,'.nc'),'huss'); %kg/kg
                    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(mlatsz/4,mlonsz,size(tmp2,3));
                    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                    hussmean(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp3;
                end
                end
        
        
                %SSP5-8.5
                %Read only 2015-2023
                if dofut==1
                tmp1=ncread(strcat(miroc6dir,'tasmax_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'tasmax')-273.15; %C
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(mlatsz/4,mlonsz,size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                tasmax(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp3;
        
                tmp1=ncread(strcat(miroc6dir,'huss_day_MIROC6_ssp585_r',num2str(i),'i1p',num2str(phys),'f1_gn_20150101-20241231.nc'),'huss'); %kg/kg
                tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(mlatsz/4,mlonsz,size(tmp2,3));
                for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(l1:l2,:,sz));end;clear tmp2;
                hussmean(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp3;
                end
                
                %Calculate wet-bulb temperature using daily-max T and daily-mean q           
                arrt=tasmax;arrp=psfc_miroc6(l1:l2,:,:).*100;arrq=hussmean.*1000;
                suffixforfilename=strcat('_ensmem',num2str(i),'_loop',num2str(latloop));clear hussmean;
                calctw_matlabpythonhelper_part1;
            end
        end
        clear psfc_miroc6;
    end

    %%%%%
    %%Now, in a Jupyter notebook, run final cell of speedywetbulb_copy.ipynb [which is in home directory for ease of access]
        %This reads the saved values, computes DJ wet-bulb temperature, and saves to another netcdf file
    %Between ensemble members, the only thing that should need changing there is suffixforfilename
    %%%%%

    %Save daily Tw and split into regions
    if miroc6_dopart2==1
        myfname='tw_miroc6_daily'; %prefix to be used for saving final output
        for i=6:miroc6_enssz
            suffixes={strcat('_ensmem',num2str(i),'_loop1');strcat('_ensmem',num2str(i),'_loop2');...
                strcat('_ensmem',num2str(i),'_loop3');strcat('_ensmem',num2str(i),'_loop4')}; %files to be processed for this ens mem

            tw_miroc6_thisensmem=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);
            calctw_matlabpythonhelper_part2;
            for latloop=1:4
                l1=l1s(latloop);l2=l2s(latloop);

                tw_miroc6_thisensmem(l1:l2,:,:)=load(strcat(saveloc,myfname,suffixes{latloop},'.mat')).tw2m;
            end

            for y=1:numyr
                %Separate daily data into regions
                tw_regdata=cell(numreg);
                c=zeros(numreg,1);
                for latc=1:mlatsz
                    for lonc=1:mlonsz
                        thisreg=round(regsz(latc,lonc));
                        if ~isnan(thisreg)
                            c(thisreg)=c(thisreg)+1;
                            tw_regdata{thisreg}(c(thisreg),:)=tw_miroc6_thisensmem(latc,lonc,y*365-364:y*365); %lat | lon | days of year
                        end
                    end
                end
                %Spatial mean for each region and day
                for r=1:numreg;tw_regdailythisyear_miroc6(r,:)=squeeze(mean(squeeze(tw_regdata{r})));end
    
                %Adjust to a GMST of +1C
                adjfactor_thisyr=gmst_miroc6_difffrom1(y);
                tw_regdaily_adj_miroc6(i,y,:,:)=tw_regdailythisyear_miroc6+repmat(adjfactor_thisyr,[numreg 365]);
                clear tw_regdailythisyear_miroc6;
            end

            disp(i);disp(clock);clear tw_miroc6_thisensmem;
        end

        %Convert to std anoms (10 sec)
        tw_vals_adj_miroc6=tw_regdaily_adj_miroc6;clear tw_regdaily_adj_miroc6;
        for i=1:miroc6_enssz
            for reg=1:numreg
                for m=1:12
                    monthmeantw(reg,m)=mean(reshape(squeeze(tw_vals_adj_miroc6(i,firstclimoyr:lastclimoyr,reg,monthstarts(m):monthstops(m))),[monthlens(m)*numclimoyr 1]),'omitnan');
                end
                [~,peakmonth_tw_miroc6(i,reg)]=max(monthmeantw(reg,:));
            end
    
            %Calculate mean warm-season daily-max temperature, wet-bulb temperature, and specific humidity
            peakmonth_tw=squeeze(peakmonth_tw_miroc6(i,:));tw_vals=permute(squeeze(tw_vals_adj_miroc6(i,:,:,:)),[2 3 1]);
        
            varstocalc={'tw'};
            firstclimoyear=firstclimoyr;lastclimoyear=lastclimoyr;numclimoyears=numclimoyr;
            meanwarmseasondailymaxhelper;
        
            meanwsmax_tw_miroc6(i,:)=meanwsmax_tw;stdwsmax_tw_miroc6(i,:)=stdwsmax_tw;
        end
       
        %Calculate daily extreme index for all years, for Tw, using Thompson et al. method
        clear tw_stdanoms_adj_miroc6;
        for i=1:miroc6_enssz
            for reg=1:numreg
                tw_stdanoms_adj_miroc6(i,:,reg,:)=(squeeze(tw_vals_adj_miroc6(i,:,reg,:))'-squeeze(repmat(meanwsmax_tw_miroc6(i,reg),[1 365 numyr])))./...
                    squeeze(repmat(stdwsmax_tw_miroc6(i,reg),[1 365 numyr]));
            end
        end
        save(strcat(saveloc,'tw_valsandstdanoms_adj_miroc6'),'tw_vals_adj_miroc6','tw_stdanoms_adj_miroc6');

        %Calc longest duration above 1, 2, and 3 std anoms
        for loop=1:3
            stdanomcutoff_here=loop;
    
            longestdurtoplot=zeros(miroc6_enssz,mlatsz,mlonsz); %just to initialize
            for i=1:miroc6_enssz
                tmp=0.*regs_miroc6sz;
                for reg=1:numreg
                    longestdursofar=0;thisdur=0;
                    for y=1:numyr
                        for doy=1:365
                            if tw_stdanoms_adj_miroc6(i,doy,reg,y)>=stdanomcutoff_here %extr heat starts or continues
                                thisdur=thisdur+1;
                            else
                                if thisdur>longestdursofar
                                    longestdursofar=thisdur;
                                    longestdur_date(reg,1)=y+firstyear-1;longestdur_date(reg,2)=doy-1; %end date of sequence
                                end
                                thisdur=0;
                            end
                        end
                    end
                    longestdurbyreg(reg)=longestdursofar;
            
                    tmp(regs_miroc6sz==reg)=longestdurbyreg(reg);
                end
                longestdurtoplot(i,:,:)=tmp;
            end
            %Ensemble-mean thereof
            if loop==1
                longestdurtoplot_miroc61=squeeze(mean(longestdurtoplot));longestdurtoplot_miroc61_byensmem=longestdurtoplot;
            elseif loop==2
                longestdurtoplot_miroc62=squeeze(mean(longestdurtoplot));longestdurtoplot_miroc62_byensmem=longestdurtoplot;
            elseif loop==3
                longestdurtoplot_miroc63=squeeze(mean(longestdurtoplot));longestdurtoplot_miroc63_byensmem=longestdurtoplot;
            end
        end
    
        save(strcat(saveloc,'longestduration_models.mat'),...
            'longestdurtoplot_miroc61','longestdurtoplot_miroc62','longestdurtoplot_miroc63',...
            'longestdurtoplot_miroc61_byensmem','longestdurtoplot_miroc62_byensmem','longestdurtoplot_miroc63_byensmem','-append');
    end
end


%Partially reread MPI-GE to save daily data
%Also compute some stats for later use
if rereadmpige_savedaily==1
    z3=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_mpige,lat2d_mpige);
    psfc_tmp=pressurefromelev(z3);
    psfc_mpige=repmat(psfc_tmp,[1 1 365]);
    mlatsz=size(lat2d_mpige,1);mlonsz=size(lat2d_mpige,2);regsz=regs_mpigesz;

    tw_regdaily_adj_mpige=NaN.*ones(nummpigeensmemstodo,numyr,numreg,365);

    for i=1:nummpigeensmemstodo
        tasmax=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);hussmean=NaN.*ones(mlatsz,mlonsz,numyrtodo*365);

        %Historical
        %Read 1961-2014
        if dohist==1
            filepath_tas_part1=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');
            filepath_huss_part1=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_historical_r',num2str(i),'i1p1f1_gn_');

            decades={'19500101-19691231';'19700101-19891231';'19900101-20091231';'20100101-20141231'};
            arrayindices_start=[1;365*9+1;365*29+1;365*49+1];arrayindices_stop=[365*9;365*29;365*49;365*54]; %where data will go
            fileindices_start=[365*11+1;1;1;1];fileindices_stop=[365*20;365*20;365*20;365*5]; %where data comes from

            histpart=1;futpart=0;
            readmodeldata_helper;
        end


        %SSP5-8.5
        %Read only 2015-2023 for now
        if dofut==1
            filepath_tas=strcat(mpigedir,'tasmax_day/tasmax_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');
            filepath_huss=strcat(mpigedir,'huss_day/huss_day_MPI-ESM1-2-LR_ssp585_r',num2str(i),'i1p1f1_gn_20150101-20341231.nc');

            histpart=0;futpart=1;
            readmodeldata_helper;
        end

        %Calculate Tw (using daily-max T and daily-mean q) (calc time 3 min per ensemble member)
        firstyrhere=1;
        for y=1:numyr
            tw_regdailythisyear_mpige=NaN.*ones(numreg,365);

            thisrelyr=y-firstyrhere+1;
            ystart=365*thisrelyr-364;yend=365*thisrelyr;
            tw_mpige_thisensmem=calcwbt_daviesjones(tasmax(:,:,ystart:yend),psfc_mpige.*100,hussmean(:,:,ystart:yend));

            %Separate into regions
            tw_regdata=cell(numreg);
            c=zeros(numreg,1);
            for latc=1:mlatsz
                for lonc=1:mlonsz
                    thisreg=round(regsz(latc,lonc));
                    if ~isnan(thisreg)
                        c(thisreg)=c(thisreg)+1;
                        tw_regdata{thisreg}(c(thisreg),:)=tw_mpige_thisensmem(latc,lonc,:); %lat | lon | days of year
                    end
                end
            end
            %Spatial mean for each region and day
            for r=1:numreg;tw_regdailythisyear_mpige(r,:)=squeeze(mean(squeeze(tw_regdata{r})));end

            %Adjust to a GMST of +1C
            adjfactor_thisyr=gmst_mpige_difffrom1(y);
            tw_regdaily_adj_mpige(i,y,:,:)=tw_regdailythisyear_mpige+repmat(adjfactor_thisyr,[numreg 365]);
            clear tw_regdailythisyear_mpige;
        end
        clear tasmax;clear hussmean;
        disp(i);disp(clock);clear tw_mpige_thisensmem;
    end
    clear psfc_mpige;

    %Convert to std anoms
    tw_vals_adj_mpige=tw_regdaily_adj_mpige;clear tw_regdaily_adj_mpige;
    for i=1:nummpigeensmemstodo
        for reg=1:numreg
            for m=1:12
                monthmeantw(reg,m)=mean(reshape(squeeze(tw_vals_adj_mpige(i,firstclimoyr:lastclimoyr,reg,monthstarts(m):monthstops(m))),[monthlens(m)*numclimoyr 1]),'omitnan');
            end
            [~,peakmonth_tw_mpige(i,reg)]=max(monthmeantw(reg,:));
        end

        %Calculate mean warm-season daily-max temperature, wet-bulb temperature, and specific humidity
        peakmonth_tw=squeeze(peakmonth_tw_mpige(i,:));tw_vals=permute(squeeze(tw_vals_adj_mpige(i,:,:,:)),[2 3 1]);
    
        varstocalc={'tw'};
        firstclimoyear=firstclimoyr;lastclimoyear=lastclimoyr;numclimoyears=numclimoyr;
        meanwarmseasondailymaxhelper;
    
        meanwsmax_tw_mpige(i,:)=meanwsmax_tw;stdwsmax_tw_mpige(i,:)=stdwsmax_tw;
    end
   
    %Calculate daily extreme index for all years, for Tw, using Thompson et al. method
    clear tw_stdanoms_adj_mpige;
    for i=1:nummpigeensmemstodo
        for reg=1:numreg
            tw_stdanoms_adj_mpige(i,:,reg,:)=(squeeze(tw_vals_adj_mpige(i,:,reg,:))'-squeeze(repmat(meanwsmax_tw_mpige(i,reg),[1 365 numyr])))./...
                squeeze(repmat(stdwsmax_tw_mpige(i,reg),[1 365 numyr]));
        end
    end
    save(strcat(saveloc,'tw_valsandstdanoms_adj_mpige'),'tw_vals_adj_mpige','tw_stdanoms_adj_mpige');

    %Calc longest duration above 1, 2, and 3 std anoms
    for loop=1:3
        stdanomcutoff_here=loop;

        longestdurtoplot=zeros(nummpigeensmemstodo,mlatsz,mlonsz); %just to initialize
        for i=1:nummpigeensmemstodo
            tmp=0.*regs_mpigesz;
            for reg=1:numreg
                longestdursofar=0;thisdur=0;
                for y=1:numyr
                    for doy=1:365
                        if tw_stdanoms_adj_mpige(i,doy,reg,y)>=stdanomcutoff_here %extr heat starts or continues
                            thisdur=thisdur+1;
                        else
                            if thisdur>longestdursofar
                                longestdursofar=thisdur;
                                longestdur_date(reg,1)=y+firstyear-1;longestdur_date(reg,2)=doy-1; %end date of sequence
                            end
                            thisdur=0;
                        end
                    end
                end
                longestdurbyreg(reg)=longestdursofar;
        
                tmp(regs_mpigesz==reg)=longestdurbyreg(reg);
            end
            longestdurtoplot(i,:,:)=tmp;
        end
        %Ensemble-mean thereof
        if loop==1
            longestdurtoplot_mpige1=squeeze(mean(longestdurtoplot));longestdurtoplot_mpige1_byensmem=longestdurtoplot;
        elseif loop==2
            longestdurtoplot_mpige2=squeeze(mean(longestdurtoplot));longestdurtoplot_mpige2_byensmem=longestdurtoplot;
        elseif loop==3
            longestdurtoplot_mpige3=squeeze(mean(longestdurtoplot));longestdurtoplot_mpige3_byensmem=longestdurtoplot;
        end
    end

    save(strcat(saveloc,'longestduration_models.mat'),...
        'longestdurtoplot_mpige1','longestdurtoplot_mpige2','longestdurtoplot_mpige3',...
        'longestdurtoplot_mpige1_byensmem','longestdurtoplot_mpige2_byensmem','longestdurtoplot_mpige3_byensmem','-append');
end



if processera5==1
    %Starting from scratch
    %tw_regannmax_era5=NaN.*ones(numreg,numyr); %reg number | year
    %tw_doyofregannmax_era5=NaN.*ones(numreg,numyr); %reg number | year
    %tw_regdecvals_era5=NaN.*ones(numreg,365,numyr);t_regdecvals_era5=NaN.*ones(numreg,365,numyr);q_regdecvals_era5=NaN.*ones(numreg,365,numyr);

    %Just adding a few years
    %yrstoadd=3;
    %tw_regdecvals_era5=cat(3,tw_regdecvals_era5,NaN.*ones(numreg,365,yrstoadd));
    %t_regdecvals_era5=cat(3,t_regdecvals_era5,NaN.*ones(numreg,365,yrstoadd));
    %q_regdecvals_era5=cat(3,q_regdecvals_era5,NaN.*ones(numreg,365,yrstoadd));

    %Static elevation (~surface pressure) array
    z0=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_era5,lat2d_era5);
    z0(:,1437:1440)=z0(:,1433:1436); %it's nearly all ocean at this lon anyway, so it makes little difference
    psfc_tmp=pressurefromelev(z0);
    psfc_era5=repmat(psfc_tmp,[1 1 365]);


    for y=firstyear:lastyear
    %for y=2021:2023
        if y<=1999;era5dir=era5dir_1;else;era5dir=era5dir_2;end
        relyr=mod(y,10);if relyr==0;relyr=10;end
        if y>=2021;numyrsindec=3;else;numyrsindec=10;end %deal with partial (current) decade

        if rem(y,10)==1            
            tw_era5_decade=NaN.*ones(numreg,365,numyrsindec);
            t_era5_decade=NaN.*ones(numreg,365,numyrsindec);
            q_era5_decade=NaN.*ones(numreg,365,numyrsindec);
        end
        

        %Read in temperature (1 min 30 sec)
        if y==2023
            t=squeeze(ncread(strcat(era5dir,'t2m_',num2str(y),'.nc'),'t2m',[1 1 1 1],[Inf Inf 1 Inf]))-273.15; %C
        else
            t=ncread(strcat(era5dir,'t2m_',num2str(y),'.nc'),'t2m')-273.15; %C
        end

        %Get daily-max T
        tdailymax=NaN.*ones(latsz_era5,lonsz_era5,365);
        for day=1:365;tdailymax(:,:,day)=recenter(squeeze(max(t(:,:,day*8-7:day*8),[],3))');end
        clear t;

        %Read in dewpoint temperature (1 min 30 sec)
        if y==2023
            td=squeeze(ncread(strcat(era5dir,'td2m_',num2str(y),'.nc'),'d2m',[1 1 1 1],[Inf Inf 1 Inf]))-273.15; %C
        else
            td=ncread(strcat(era5dir,'td2m_',num2str(y),'.nc'),'d2m')-273.15; %C
        end

        %Get daily-mean Td 
        tddailymean=NaN.*ones(latsz_era5,lonsz_era5,365);
        for day=1:365;tddailymean(:,:,day)=recenter(squeeze(mean(td(:,:,day*8-7:day*8),3))');end
        clear td;


        %To save time:
        %Before computing Tw, get daily regional means of T and Td (45 sec)
        tregmeans=cell(numreg,1);tdregmeans=cell(numreg,1);psfcregmeans=cell(numreg,1);c=zeros(numreg,1); %237 regions
        latweights=cell(numreg,1);
        for i=1:latsz_era5
            for j=1:lonsz_era5
                thisreg=round(regs_era5sz(i,j));
                if ~isnan(thisreg)
                    c(thisreg)=c(thisreg)+1;
                    tregmeans{thisreg}(c(thisreg),:)=tdailymax(i,j,:);
                    tdregmeans{thisreg}(c(thisreg),:)=tddailymean(i,j,:);
                    psfcregmeans{thisreg}(c(thisreg))=psfc_era5(i,j);

                    latweights{thisreg}(c(thisreg))=cos(deg2rad(lat2d_era5(i,j)));
                end
            end
        end
        clear tdailymax;clear tddailymean;

        %Spatial mean by region
        tregmean=NaN.*ones(numreg,365);tdregmean=NaN.*ones(numreg,365);psfcregmean=NaN.*ones(numreg,365);
        for r=1:numreg
            if latwgt_era5==0 %no lat weighting (original)
                tregmean(r,:)=mean(tregmeans{r});tdregmean(r,:)=mean(tdregmeans{r});psfcregmean(r,:)=mean(psfcregmeans{r});
            else %latitudinally weighted
                weights=repmat(latweights{r}',[1 365]);
                for doy=1:365
                    tregmean(r,doy)=sum(weights(:,doy).*tregmeans{r}(:,doy))./sum(weights(:,doy));
                    tdregmean(r,doy)=sum(weights(:,doy).*tdregmeans{r}(:,doy))./sum(weights(:,doy));
                    psfcregmean(r,doy)=sum(weights(:,doy).*psfcregmeans{r})./sum(weights(:,doy));
                end
            end
        end
        clear tregmeans;clear tdregmeans;clear psfcregmeans;

        arrp=psfcregmean.*100;arrt=tregmean;arrq=calcqfromTd_dynamicP(tdregmean,arrp)./1000;
        tw_era5=calcwbt_daviesjones(arrt,arrp,arrq);

        %Annual-max Tw by region, and the date on which it occurred
        tw_regannmax_era5(:,y-(firstyear-1))=squeeze(max(tw_era5,[],2));
        for r=1:numreg
            [~,tw_doyofregannmax_era5(r,y-(firstyear-1))]=max(tw_era5(r,:));
        end


        %Save into decadal array
        tw_era5_decade(:,:,relyr)=tw_era5;
        t_era5_decade(:,:,relyr)=arrt;
        q_era5_decade(:,:,relyr)=arrq;

        if relyr==numyrsindec
            %Save big arrays of actual values, for inspection and actual later calculations
            tw_regdecvals_era5(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=tw_era5_decade;
            t_regdecvals_era5(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=t_era5_decade;
            q_regdecvals_era5(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=q_era5_decade;

            save(strcat(saveloc,'mylatestera5output_regions.mat'),'tw_regannmax_era5','tw_doyofregannmax_era5', ...
                'tw_regdecvals_era5','t_regdecvals_era5','q_regdecvals_era5');
        end
       
        disp(y);disp(clock);
    end

    %Remove clearly bad moisture data in many regions
    q_regdecvals_era5(:,348:365,22)=NaN;tw_regdecvals_era5(:,348:365,22)=NaN;
    q_regdecvals_era5(:,362:365,23)=NaN;tw_regdecvals_era5(:,362:365,23)=NaN;
    q_regdecvals_era5(:,355:365,24)=NaN;tw_regdecvals_era5(:,355:365,24)=NaN;

    %Dec 2023 is also bad for some reason
    t_regdecvals_era5(:,335:365,63)=NaN;q_regdecvals_era5(:,335:365,63)=NaN;tw_regdecvals_era5(:,335:365,63)=NaN;

    save(strcat(saveloc,'mylatestera5output_regions.mat'),'tw_regannmax_era5','tw_doyofregannmax_era5', ...
        'tw_regdecvals_era5','t_regdecvals_era5','q_regdecvals_era5');

    %Write to text files for reading in Python (so GMST-adjusted versions can be calculated)
    clear t_regannmax_era5;clear tw_regannmax_era5;clear t_regannmax_adj_era5;clear tw_regannmax_adj_era5;
    for reg=1:numreg
        t_regannmax_era5(reg,:)=max(squeeze(t_regdecvals_era5(reg,:,:)));
        tw_regannmax_era5(reg,:)=max(squeeze(tw_regdecvals_era5(reg,:,:)));
    end
    writematrix(round2(t_regannmax_era5,0.001),strcat(saveloc,'era5output_t_regions.txt'));
    writematrix(round2(tw_regannmax_era5,0.001),strcat(saveloc,'era5output_tw_regions.txt'));

    save(strcat(saveloc,'regannmaxes'),'t_regannmax_era5','tw_regannmax_era5','-append');


    %Also, save in a text file to facilitate reloading in Python
    writematrix(round2(tw_regannmax_era5,0.001),strcat(saveloc,'era5output_tw_regions.txt'));
    writematrix(round2(tw_doyofregannmax_era5,0.001),strcat(saveloc,'era5doyofmax_tw_regions.txt'));
end



if processjra55==1
    %Starting from scratch
    tw_regannmax_jra55=NaN.*ones(numreg,numyr); %reg number | year
    tw_doyofregannmax_jra55=NaN.*ones(numreg,numyr); %reg number | year
    tw_regdecvals_jra55=NaN.*ones(numreg,365,numyr);t_regdecvals_jra55=NaN.*ones(numreg,365,numyr);q_regdecvals_jra55=NaN.*ones(numreg,365,numyr);

    %Static elevation (~surface pressure) array
    z0=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_jra55,lat2d_jra55);
    z0(:,639:640)=z0(:,637:638); %it's nearly all ocean at this lon anyway, so it makes little difference
    psfc_tmp=pressurefromelev(z0);
    psfc_jra55=repmat(psfc_tmp,[1 1 365]);


    for y=firstyear:lastyear
        relyr=mod(y,10);if relyr==0;relyr=10;end
        if y>=2021;numyrsindec=3;else;numyrsindec=10;end %deal with partial (current) decade

        if rem(y,10)==1            
            tw_jra55_decade=NaN.*ones(numreg,365,numyrsindec);
            t_jra55_decade=NaN.*ones(numreg,365,numyrsindec);
            q_jra55_decade=NaN.*ones(numreg,365,numyrsindec);
        end
        

        %Read in temperature (1 min 30 sec)
        if y<=2013
            t=ncread(strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'010100_',num2str(y),'123118.raymond741169.nc'),'TMP_GDS4_HTGL')-273.15; %C
        else %recent years are contained in monthly files
            jan=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'010100_',num2str(y),'013118.raymond741169.nc');
            if rem(y,4)==0
                feb=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'020100_',num2str(y),'022918.raymond741169.nc');
            else
                feb=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'020100_',num2str(y),'022818.raymond741169.nc');
            end
            mar=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'030100_',num2str(y),'033118.raymond741169.nc');
            apr=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'040100_',num2str(y),'043018.raymond741169.nc');
            may=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'050100_',num2str(y),'053118.raymond741169.nc');
            jun=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'060100_',num2str(y),'063018.raymond741169.nc');
            jul=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'070100_',num2str(y),'073118.raymond741169.nc');
            aug=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'080100_',num2str(y),'083118.raymond741169.nc');
            sep=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'090100_',num2str(y),'093018.raymond741169.nc');
            oct=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'100100_',num2str(y),'103118.raymond741169.nc');
            nov=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'110100_',num2str(y),'113018.raymond741169.nc');
            dec=strcat(jra55dir,'anl_surf.011_tmp.reg_tl319.',num2str(y),'120100_',num2str(y),'123118.raymond741169.nc');

            t=cat(3,ncread(jan,'TMP_GDS4_HTGL'),ncread(feb,'TMP_GDS4_HTGL'),ncread(mar,'TMP_GDS4_HTGL'),ncread(apr,'TMP_GDS4_HTGL'),...
                ncread(may,'TMP_GDS4_HTGL'),ncread(jun,'TMP_GDS4_HTGL'),ncread(jul,'TMP_GDS4_HTGL'),ncread(aug,'TMP_GDS4_HTGL'),...
                ncread(sep,'TMP_GDS4_HTGL'),ncread(oct,'TMP_GDS4_HTGL'),ncread(nov,'TMP_GDS4_HTGL'),ncread(dec,'TMP_GDS4_HTGL'))-273.15; %C
        end

        %Get daily-max T (from 6-hourly data)
        tdailymax=NaN.*ones(latsz_jra55,lonsz_jra55,365);
        for day=1:365;tdailymax(:,:,day)=recenter(squeeze(max(t(:,:,day*4-3:day*4),[],3))');end
        clear t;

        %Read in specific humidity (1 min 30 sec)
        if y<=2013
            q=ncread(strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'010100_',num2str(y),'123118.raymond741169.nc'),'SPFH_GDS4_HTGL'); %kg/kg
        else %recent years are contained in monthly files
            jan=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'010100_',num2str(y),'013118.raymond741169.nc');
            if rem(y,4)==0
                feb=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'020100_',num2str(y),'022918.raymond741169.nc');
            else
                feb=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'020100_',num2str(y),'022818.raymond741169.nc');
            end
            mar=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'030100_',num2str(y),'033118.raymond741169.nc');
            apr=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'040100_',num2str(y),'043018.raymond741169.nc');
            may=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'050100_',num2str(y),'053118.raymond741169.nc');
            jun=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'060100_',num2str(y),'063018.raymond741169.nc');
            jul=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'070100_',num2str(y),'073118.raymond741169.nc');
            aug=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'080100_',num2str(y),'083118.raymond741169.nc');
            sep=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'090100_',num2str(y),'093018.raymond741169.nc');
            oct=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'100100_',num2str(y),'103118.raymond741169.nc');
            nov=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'110100_',num2str(y),'113018.raymond741169.nc');
            dec=strcat(jra55dir,'anl_surf.051_spfh.reg_tl319.',num2str(y),'120100_',num2str(y),'123118.raymond741169.nc');

            q=cat(3,ncread(jan,'SPFH_GDS4_HTGL'),ncread(feb,'SPFH_GDS4_HTGL'),ncread(mar,'SPFH_GDS4_HTGL'),ncread(apr,'SPFH_GDS4_HTGL'),...
                ncread(may,'SPFH_GDS4_HTGL'),ncread(jun,'SPFH_GDS4_HTGL'),ncread(jul,'SPFH_GDS4_HTGL'),ncread(aug,'SPFH_GDS4_HTGL'),...
                ncread(sep,'SPFH_GDS4_HTGL'),ncread(oct,'SPFH_GDS4_HTGL'),ncread(nov,'SPFH_GDS4_HTGL'),ncread(dec,'SPFH_GDS4_HTGL')); %kg/kg
        end

        %Get daily-mean q (from 6-hourly data)
        qdailymean=NaN.*ones(latsz_jra55,lonsz_jra55,365);
        for day=1:365;qdailymean(:,:,day)=recenter(squeeze(mean(q(:,:,day*4-3:day*4),3))');end
        clear q;


        %To save time:
        %Before computing Tw, get daily regional means of T and q (30 sec)
        tregmeans=cell(numreg,1);qregmeans=cell(numreg,1);psfcregmeans=cell(numreg,1);c=zeros(numreg,1); %237 regions
        latweights=cell(numreg,1);
        for i=1:latsz_jra55
            for j=1:lonsz_jra55
                thisreg=round(regs_jra55sz(i,j));
                if ~isnan(thisreg)
                    c(thisreg)=c(thisreg)+1;
                    tregmeans{thisreg}(c(thisreg),:)=tdailymax(i,j,:);
                    qregmeans{thisreg}(c(thisreg),:)=qdailymean(i,j,:);
                    psfcregmeans{thisreg}(c(thisreg))=psfc_jra55(i,j);

                    latweights{thisreg}(c(thisreg))=cos(deg2rad(lat2d_jra55(i,j)));
                end
            end
        end
        %clear tdailymax;clear qdailymean;
        %Spatial mean by region
        tregmean=NaN.*ones(numreg,365);qregmean=NaN.*ones(numreg,365);psfcregmean=NaN.*ones(numreg,365);
        for r=1:numreg
            if latwgt_jra55==0 %no lat weighting (original)
                tregmean(r,:)=mean(tregmeans{r});qregmean(r,:)=mean(qregmeans{r});psfcregmean(r,:)=mean(psfcregmeans{r});
            else %latitudinally weighted
                weights=repmat(latweights{r}',[1 365]);
                for doy=1:365
                    tregmean(r,doy)=sum(weights(:,doy).*tregmeans{r}(:,doy))./sum(weights(:,doy));
                    qregmean(r,doy)=sum(weights(:,doy).*qregmeans{r}(:,doy))./sum(weights(:,doy));
                    psfcregmean(r,doy)=sum(weights(:,doy).*psfcregmeans{r}')./sum(weights(:,doy));
                end
            end
        end
        clear tregmeans;clear qregmeans;

        arrp=psfcregmean.*100;arrt=tregmean;arrq=qregmean;
        tw_jra55=calcwbt_daviesjones(arrt,arrp,arrq);

        %Annual-max Tw by region, and the date on which it occurred
        tw_regannmax_jra55(:,y-(firstyear-1))=squeeze(max(tw_jra55,[],2));
        for r=1:numreg
            [~,tw_doyofregannmax_jra55(r,y-(firstyear-1))]=max(tw_jra55(r,:));
        end


        %Save into decadal array
        tw_jra55_decade(:,:,relyr)=tw_jra55;
        t_jra55_decade(:,:,relyr)=arrt;
        q_jra55_decade(:,:,relyr)=arrq;

        if relyr==numyrsindec
            %Save big arrays of actual values, for inspection and actual later calculations
            tw_regdecvals_jra55(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=tw_jra55_decade;
            t_regdecvals_jra55(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=t_jra55_decade;
            q_regdecvals_jra55(:,:,y-(firstyear-1)-(numyrsindec-1):y-(firstyear-1))=q_jra55_decade;

            save(strcat(saveloc,'mylatestjra55output_regions.mat'),'tw_regannmax_jra55','tw_doyofregannmax_jra55', ...
                'tw_regdecvals_jra55','t_regdecvals_jra55','q_regdecvals_jra55');
        end
       
        disp(y);disp(clock);
    end


    save(strcat(saveloc,'mylatestjra55output_regions.mat'),'tw_regannmax_jra55','tw_doyofregannmax_jra55', ...
        'tw_regdecvals_jra55','t_regdecvals_jra55','q_regdecvals_jra55');

    %Write to text files for reading in Python (so GMST-adjusted versions can be calculated)
    clear t_regannmax_jra55;clear tw_regannmax_jra55;clear t_regannmax_adj_jra55;clear tw_regannmax_adj_jra55;
    for reg=1:numreg
        t_regannmax_jra55(reg,:)=max(squeeze(t_regdecvals_jra55(reg,:,:)));
        tw_regannmax_jra55(reg,:)=max(squeeze(tw_regdecvals_jra55(reg,:,:)));
    end
    writematrix(round2(t_regannmax_jra55,0.001),strcat(saveloc,'jra55output_t_regions.txt'));
    writematrix(round2(tw_regannmax_jra55,0.001),strcat(saveloc,'jra55output_tw_regions.txt'));

    save(strcat(saveloc,'regannmaxes'),'t_regannmax_jra55','tw_regannmax_jra55','-append');


    %Also, save in a text file to facilitate reloading in Python
    writematrix(round2(tw_regannmax_jra55,0.001),strcat(saveloc,'jra55output_tw_regions.txt'));
    writematrix(round2(tw_doyofregannmax_jra55,0.001),strcat(saveloc,'jra55doyofmax_tw_regions.txt'));
end

%As in Thompson et al. 2022, 2023, adjust to a GMST of 1C
if adjustgmst==1
    gmst_obs_difffrom1=1-gmst_obs;

    %GMST adjustments are calculated in Python script (REDO IF NECESSARY); then, import them here
    %Factors vary by region and by year, naturally being most positive at the start of the timeseries and zero or slightly negative at the end
    era5adjfactors_tw=csvread(strcat(saveloc,'era5adjfactors_tw.csv'));
    merra2adjfactors_tw=csvread(strcat(saveloc,'merra2adjfactors_tw.csv'));
    jra55adjfactors_tw=csvread(strcat(saveloc,'jra55adjfactors_tw.csv'));

    canesm5adjfactors_tw=reshape(csvread(strcat(saveloc,'canesm5adjfactors_tw.csv')),[canesm5_enssz numreg numyr]);
    miroc6adjfactors_tw=reshape(csvread(strcat(saveloc,'miroc6adjfactors_tw.csv')),[miroc6_enssz numreg numyr]);
    mpigeadjfactors_tw=reshape(csvread(strcat(saveloc,'mpigeadjfactors_tw.csv')),[mpige_enssz numreg numyr]);

    %Note: T can only be run after having done Python GMST Tw adjustment,
    %then running code that creates e.g. t_regannmax_adj_jra55
    era5adjfactors_t=csvread(strcat(saveloc,'era5adjfactors_t.csv'));
    merra2adjfactors_t=csvread(strcat(saveloc,'merra2adjfactors_t.csv'));
    jra55adjfactors_t=csvread(strcat(saveloc,'jra55adjfactors_t.csv'));

    tw_regdecvals_adj_era5=NaN.*ones(size(tw_regdecvals_era5));t_regdecvals_adj_era5=NaN.*ones(size(t_regdecvals_era5));
    tw_regdecvals_adj_merra2=NaN.*ones(size(tw_regdecvals_merra2));t_regdecvals_adj_merra2=NaN.*ones(size(t_regdecvals_merra2));
    tw_regdecvals_adj_jra55=NaN.*ones(size(tw_regdecvals_jra55));t_regdecvals_adj_jra55=NaN.*ones(size(t_regdecvals_jra55));
    %tw_regdecvals_adj_canesm5=NaN.*ones(size(tw_regdecvals_canesm5));
    %tw_regdecvals_adj_miroc6=NaN.*ones(size(tw_regdecvals_miroc6));
    %tw_regdecvals_adj_mpige=NaN.*ones(size(tw_regdecvals_mpige));
    for r=1:numreg
        for y=1:numyr
            myvals=squeeze(tw_regdecvals_era5(r,:,y));myadjfactor=era5adjfactors_tw(r,y);
                tw_regdecvals_adj_era5(r,:,y)=myvals+myadjfactor;
            myvals=squeeze(t_regdecvals_era5(r,:,y));myadjfactor=era5adjfactors_t(r,y);
                t_regdecvals_adj_era5(r,:,y)=myvals+myadjfactor;

            myvals=squeeze(tw_regdecvals_merra2(r,:,y));myadjfactor=merra2adjfactors_tw(r,y);
                tw_regdecvals_adj_merra2(r,:,y)=myvals+myadjfactor;
            myvals=squeeze(t_regdecvals_merra2(r,:,y));myadjfactor=merra2adjfactors_t(r,y);
                t_regdecvals_adj_merra2(r,:,y)=myvals+myadjfactor;

            myvals=squeeze(tw_regdecvals_jra55(r,:,y));myadjfactor=jra55adjfactors_tw(r,y);
                tw_regdecvals_adj_jra55(r,:,y)=myvals+myadjfactor;
            myvals=squeeze(t_regdecvals_jra55(r,:,y));myadjfactor=jra55adjfactors_t(r,y);
                t_regdecvals_adj_jra55(r,:,y)=myvals+myadjfactor;

            %myvals=squeeze(tw_regdecvals_canesm5(r,:,y));myadjfactor=canesm5adjfactors_tw(r,y);
            %    tw_regdecvals_adj_canesm5(r,:,y)=myvals+myadjfactor;
            %myvals=squeeze(tw_regdecvals_miroc6(r,:,y));myadjfactor=miroc6adjfactors_tw(r,y);
            %    tw_regdecvals_adj_miroc6(r,:,y)=myvals+myadjfactor;
            %myvals=squeeze(tw_regdecvals_mpige(r,:,y));myadjfactor=mpigeadjfactors_tw(r,y);
            %    tw_regdecvals_adj_mpige(r,:,y)=myvals+myadjfactor;
        end
    end

    clear t_regannmax_adj_era5;clear tw_regannmax_adj_era5;
    clear t_regannmax_adj_merra2;clear tw_regannmax_adj_merra2;
    clear t_regannmax_adj_jra55;clear tw_regannmax_adj_jra55;
    clear tw_regannmax_adj_canesm5;clear tw_regannmax_adj_miroc6;clear tw_regannmax_adj_mpige;
    for reg=1:numreg
        t_regannmax_adj_era5(reg,:)=max(squeeze(t_regdecvals_adj_era5(reg,:,:)));
        tw_regannmax_adj_era5(reg,:)=max(squeeze(tw_regdecvals_adj_era5(reg,:,:)));
        t_regannmax_adj_merra2(reg,:)=max(squeeze(t_regdecvals_adj_merra2(reg,:,:)));
        tw_regannmax_adj_merra2(reg,:)=max(squeeze(tw_regdecvals_adj_merra2(reg,:,:)));
        t_regannmax_adj_jra55(reg,:)=max(squeeze(t_regdecvals_adj_jra55(reg,:,:)));
        tw_regannmax_adj_jra55(reg,:)=max(squeeze(tw_regdecvals_adj_jra55(reg,:,:)));

        %tw_regannmax_adj_canesm5(reg,:)=max(squeeze(tw_regdecvals_adj_canesm5(reg,:,:)));
        %tw_regannmax_adj_miroc6(reg,:)=max(squeeze(tw_regdecvals_adj_miroc6(reg,:,:)));
        %tw_regannmax_adj_mpige(reg,:)=max(squeeze(tw_regdecvals_adj_mpige(reg,:,:)));
    end

    save(strcat(saveloc,'regannmaxes'),'t_regannmax_adj_era5','tw_regannmax_adj_era5',...
        't_regannmax_adj_merra2','tw_regannmax_adj_merra2','t_regannmax_adj_jra55','tw_regannmax_adj_jra55','-append');

    writematrix(round2(tw_regannmax_adj_era5,0.001),strcat(saveloc,'era5output_tw_adj_regions.txt'));
    writematrix(round2(tw_regannmax_adj_merra2,0.001),strcat(saveloc,'merra2output_tw_adj_regions.txt'));
    writematrix(round2(tw_regannmax_adj_jra55,0.001),strcat(saveloc,'jra55output_tw_adj_regions.txt'));

    %save(strcat(saveloc,'regannmaxes_models'),'tw_regannmax_adj_canesm5','tw_regannmax_adj_miroc6','tw_regannmax_adj_mpige');

    %writematrix(round2(tw_regannmax_adj_canesm5,0.001),strcat(saveloc,'canesm5output_tw_adj_regions.txt'));
    %writematrix(round2(tw_regannmax_adj_miroc6,0.001),strcat(saveloc,'miroc6output_tw_adj_regions.txt'));
    %writematrix(round2(tw_regannmax_adj_mpige,0.001),strcat(saveloc,'mpigeoutput_tw_adj_regions.txt'));
end


if era5_calcstdanoms==1
    %Compute month with greatest mean temperature and mean wet-bulb temperature
    %**Now, uses data adjusted to a GMST of +1C, as in Thompson et al. papers**
    %Reminder: climatologies are based on 1981-2010, as in Thompson et al. 2022
    choice=2;
    if choice==1 %original -- for SI figs only
        tw_vals_era5=tw_regdecvals_era5;t_vals_era5=t_regdecvals_era5;q_vals_era5=q_regdecvals_era5;
    elseif choice==2 %updated; q is not directly adjusted but this isn't used in any of the final figures or analysis anyway
        tw_vals_era5=tw_regdecvals_adj_era5;t_vals_era5=t_regdecvals_adj_era5;q_vals_era5=q_regdecvals_era5;
    end

    for reg=1:numreg
        for m=1:12
            monthmeantw(reg,m)=mean(reshape(squeeze(tw_vals_era5(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
            monthmeant(reg,m)=mean(reshape(squeeze(t_vals_era5(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
            monthmeanq(reg,m)=mean(reshape(squeeze(q_vals_era5(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
        end
        [~,peakmonth_tw_era5(reg)]=max(monthmeantw(reg,:));
        [~,peakmonth_t_era5(reg)]=max(monthmeant(reg,:));
        [~,peakmonth_q_era5(reg)]=max(monthmeanq(reg,:));
    end

    %Get annual-mean RH for each location
    td_vals_era5=calcTdfromq(q_vals_era5.*1000);
    rh_vals_era5=calcrhfromTandTd(t_vals_era5,td_vals_era5);
    for reg=1:numreg
        annualmeanrh_era5(reg)=mean(reshape(squeeze(rh_vals_era5(reg,1:365,firstclimoyr:lastclimoyr)),[365*numclimoyr 1]),'omitnan');
        annualstdrh_era5(reg)=std(reshape(squeeze(rh_vals_era5(reg,1:365,firstclimoyr:lastclimoyr)),[365*numclimoyr 1]),'omitnan');
    end

    %Calculate mean warm-season daily-max temperature, wet-bulb temperature, and specific humidity
    peakmonth_tw=peakmonth_tw_era5;peakmonth_t=peakmonth_t_era5;peakmonth_q=peakmonth_q_era5;
    tw_vals=tw_vals_era5;t_vals=t_vals_era5;q_vals=q_vals_era5;

    varstocalc={'t';'tw';'q'};
    firstclimoyear=firstclimoyr;lastclimoyear=lastclimoyr;numclimoyears=numclimoyr;
    meanwarmseasondailymaxhelper;

    meanwsmax_tw_era5=meanwsmax_tw;meanwsmax_t_era5=meanwsmax_t;meanwsmax_q_era5=meanwsmax_q;
    stdwsmax_tw_era5=stdwsmax_tw;stdwsmax_t_era5=stdwsmax_t;stdwsmax_q_era5=stdwsmax_q;

    %Also calculate mean warm-season RH
    for reg=1:numreg
        month1=peakmonth_tw_era5(reg)-1;month2=peakmonth_tw_era5(reg)+1;
        if month1==0 %actually Dec
            wsmeanrh_era5(reg)=mean(cat(1,reshape(squeeze(rh_vals_era5(reg,335:365,firstclimoyr:lastclimoyr)),[31*numclimoyr 1]),...
                reshape(squeeze(rh_vals_era5(reg,1:59,firstclimoyr:lastclimoyr)),[59*numclimoyr 1])),'omitnan');
        elseif month2==13 %actually Jan
            wsmeanrh_era5(reg)=mean(cat(1,reshape(squeeze(rh_vals_era5(reg,305:365,firstclimoyr:lastclimoyr)),[61*numclimoyr 1]),...
                reshape(squeeze(rh_vals_era5(reg,1:31,firstclimoyr:lastclimoyr)),[31*numclimoyr 1])),'omitnan');
        else
            wsmeanrh_era5(reg)=mean(reshape(squeeze(rh_vals_era5(reg,monthstarts(month1):monthstops(month2),firstclimoyr:lastclimoyr)),...
                [(monthstops(month2)-monthstarts(month1)+1)*numclimoyr 1]),'omitnan');
        end
    end

    %Finally, calculate daily extreme index for all years, again for T, Tw, and q, using Thompson et al. method
    clear t_stdanoms_era5;clear q_stdanoms_era5;clear tw_stdanoms_era5;
    for reg=1:numreg
        t_stdanoms_era5(reg,:,:)=(squeeze(t_vals_era5(reg,:,:))-squeeze(repmat(meanwsmax_t_era5(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_t_era5(reg),[1 365 numyr]));
        tw_stdanoms_era5(reg,:,:)=(squeeze(tw_vals_era5(reg,:,:))-squeeze(repmat(meanwsmax_tw_era5(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_tw_era5(reg),[1 365 numyr]));
        q_stdanoms_era5(reg,:,:)=(squeeze(q_vals_era5(reg,:,:))-squeeze(repmat(meanwsmax_q_era5(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_q_era5(reg),[1 365 numyr]));

        %Get year and DOY of daily-extreme-index max in each region
        tmp=squeeze(t_stdanoms_era5(reg,:,:));
            [maxstdanom_t_tmp_era5(reg),yearofmax_t_tmp_era5(reg)]=max(max(tmp));[~,doyofmax_t_tmp_era5(reg)]=max(tmp(:,yearofmax_t_tmp_era5(reg)));
        tmp=squeeze(tw_stdanoms_era5(reg,:,:));
            [maxstdanom_tw_tmp_era5(reg),yearofmax_tw_tmp_era5(reg)]=max(max(tmp));[~,doyofmax_tw_tmp_era5(reg)]=max(tmp(:,yearofmax_tw_tmp_era5(reg)));
        tmp=squeeze(q_stdanoms_era5(reg,:,:));
            [maxstdanom_q_tmp_era5(reg),yearofmax_q_tmp_era5(reg)]=max(max(tmp));[~,doyofmax_q_tmp_era5(reg)]=max(tmp(:,yearofmax_q_tmp_era5(reg)));

        yearofmax_t_tmp_era5(reg)=yearofmax_t_tmp_era5(reg)+firstyear-1;
        yearofmax_tw_tmp_era5(reg)=yearofmax_tw_tmp_era5(reg)+firstyear-1;
        yearofmax_q_tmp_era5(reg)=yearofmax_q_tmp_era5(reg)+firstyear-1;
    end

    if choice==1
        tw_regdecstdanoms_era5=tw_stdanoms_era5;t_regdecstdanoms_era5=t_stdanoms_era5;q_regdecstdanoms_era5=q_stdanoms_era5;
        maxstdanom_tw_era5=maxstdanom_tw_tmp_era5;maxstdanom_t_era5=maxstdanom_t_tmp_era5;maxstdanom_q_era5=maxstdanom_q_tmp_era5;
        yearofmax_tw_era5=yearofmax_tw_tmp_era5;yearofmax_t_era5=yearofmax_t_tmp_era5;yearofmax_q_era5=yearofmax_q_tmp_era5;
        doyofmax_tw_era5=doyofmax_tw_tmp_era5;doyofmax_t_era5=doyofmax_t_tmp_era5;doyofmax_q_era5=doyofmax_q_tmp_era5;
        save(strcat(saveloc,'era5output_regions.mat'),'maxstdanom_tw_era5','tw_regdecstdanoms_era5',...
            'maxstdanom_t_era5','t_regdecstdanoms_era5','maxstdanom_q_era5','q_regdecstdanoms_era5',...
            'yearofmax_t_era5','doyofmax_t_era5','yearofmax_tw_era5','doyofmax_tw_era5','yearofmax_q_era5','doyofmax_q_era5',...
            'tw_regdecvals_era5','t_regdecvals_era5','wsmeanrh_era5','peakmonth_tw_era5','peakmonth_t_era5','peakmonth_q_era5','-append');
    elseif choice==2
        tw_regdecstdanoms_adj_era5=tw_stdanoms_era5;t_regdecstdanoms_adj_era5=t_stdanoms_era5;q_regdecstdanoms_adj_era5=q_stdanoms_era5;
        maxstdanom_tw_adj_era5=maxstdanom_tw_tmp_era5;maxstdanom_t_adj_era5=maxstdanom_t_tmp_era5;maxstdanom_q_adj_era5=maxstdanom_q_tmp_era5;
        yearofmax_tw_adj_era5=yearofmax_tw_tmp_era5;yearofmax_t_adj_era5=yearofmax_t_tmp_era5;yearofmax_q_adj_era5=yearofmax_q_tmp_era5;
        doyofmax_tw_adj_era5=doyofmax_tw_tmp_era5;doyofmax_t_adj_era5=doyofmax_t_tmp_era5;doyofmax_q_adj_era5=doyofmax_q_tmp_era5;
        save(strcat(saveloc,'era5output_regions.mat'),'maxstdanom_tw_adj_era5','tw_regdecstdanoms_adj_era5',...
            'maxstdanom_t_adj_era5','t_regdecstdanoms_adj_era5','maxstdanom_q_adj_era5','q_regdecstdanoms_adj_era5',...
            'yearofmax_t_adj_era5','doyofmax_t_adj_era5','yearofmax_tw_adj_era5','doyofmax_tw_adj_era5',...
            'yearofmax_q_adj_era5','doyofmax_q_adj_era5','tw_regdecvals_adj_era5','t_regdecvals_adj_era5',...
            'wsmeanrh_era5','peakmonth_tw_era5','peakmonth_t_era5','peakmonth_q_era5','-append');
    end
end





if jra55_calcstdanoms==1
    %Compute month with greatest mean temperature and mean wet-bulb temperature
    %Climatologies are based on 1981-2010, as in Thompson et al. 2022
    choice=2;
    if choice==1 %original
        tw_vals_jra55=tw_regdecvals_jra55;t_vals_jra55=t_regdecvals_jra55;q_vals_jra55=q_regdecvals_jra55;
    elseif choice==2 %updated; q is not adjusted but this isn't used in any of the final figures or analysis anyway
        tw_vals_jra55=tw_regdecvals_adj_jra55;t_vals_jra55=t_regdecvals_adj_jra55;q_vals_jra55=q_regdecvals_jra55;
    end

    for reg=1:numreg
        for m=1:12
            monthmeantw(reg,m)=mean(reshape(squeeze(tw_vals_jra55(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
            monthmeant(reg,m)=mean(reshape(squeeze(t_vals_jra55(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
            monthmeanq(reg,m)=mean(reshape(squeeze(q_vals_jra55(reg,monthstarts(m):monthstops(m),firstclimoyr:lastclimoyr)),[monthlens(m)*numclimoyr 1]),'omitnan');
        end
        [~,peakmonth_tw_jra55(reg)]=max(monthmeantw(reg,:));
        [~,peakmonth_t_jra55(reg)]=max(monthmeant(reg,:));
        [~,peakmonth_q_jra55(reg)]=max(monthmeanq(reg,:));
    end

    %Get mean RH for each place
    td_vals_jra55=calcTdfromq(q_vals_jra55.*1000);
    rh_vals_jra55=calcrhfromTandTd(t_vals_jra55,td_vals_jra55);
    for reg=1:numreg
        annualmeanrh_jra55(reg)=mean(reshape(squeeze(rh_vals_jra55(reg,1:365,firstclimoyr:lastclimoyr)),[365*numclimoyr 1]),'omitnan');
        annualstdrh_jra55(reg)=std(reshape(squeeze(rh_vals_jra55(reg,1:365,firstclimoyr:lastclimoyr)),[365*numclimoyr 1]),'omitnan');
    end

    %Calculate mean warm-season daily-max temperature, wet-bulb temperature, and specific humidity
    peakmonth_tw=peakmonth_tw_jra55;peakmonth_t=peakmonth_t_jra55;peakmonth_q=peakmonth_q_jra55;
    tw_vals=tw_vals_jra55;t_vals=t_vals_jra55;q_vals=q_vals_jra55;

    varstocalc={'t';'tw';'q'};
    firstclimoyear=firstclimoyr;lastclimoyear=lastclimoyr;numclimoyears=numclimoyr;
    meanwarmseasondailymaxhelper;

    meanwsmax_tw_jra55=meanwsmax_tw;meanwsmax_t_jra55=meanwsmax_t;meanwsmax_q_jra55=meanwsmax_q;
    stdwsmax_tw_jra55=stdwsmax_tw;stdwsmax_t_jra55=stdwsmax_t;stdwsmax_q_jra55=stdwsmax_q;


    %Finally, calculate daily extreme index for all years, again for T, Tw, and q, using Thompson et al. method
    clear t_stdanoms_jra55;clear q_stdanoms_jra55;clear tw_stdanoms_jra55;
    for reg=1:numreg
        t_stdanoms_jra55(reg,:,:)=(squeeze(t_vals_jra55(reg,:,:))-squeeze(repmat(meanwsmax_t_jra55(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_t_jra55(reg),[1 365 numyr]));
        tw_stdanoms_jra55(reg,:,:)=(squeeze(tw_vals_jra55(reg,:,:))-squeeze(repmat(meanwsmax_tw_jra55(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_tw_jra55(reg),[1 365 numyr]));
        q_stdanoms_jra55(reg,:,:)=(squeeze(q_vals_jra55(reg,:,:))-squeeze(repmat(meanwsmax_q_jra55(reg),[1 365 numyr])))./...
            squeeze(repmat(stdwsmax_q_jra55(reg),[1 365 numyr]));

        %Get year and DOY of DEI max in each region
        tmp=squeeze(t_stdanoms_jra55(reg,:,:));
            [maxstdanom_t_tmp_jra55(reg),yearofmax_t_tmp_jra55(reg)]=max(max(tmp));[~,doyofmax_t_tmp_jra55(reg)]=max(tmp(:,yearofmax_t_tmp_jra55(reg)));
        tmp=squeeze(tw_stdanoms_jra55(reg,:,:));
            [maxstdanom_tw_tmp_jra55(reg),yearofmax_tw_tmp_jra55(reg)]=max(max(tmp));[~,doyofmax_tw_tmp_jra55(reg)]=max(tmp(:,yearofmax_tw_tmp_jra55(reg)));
        tmp=squeeze(q_stdanoms_jra55(reg,:,:));
            [maxstdanom_q_tmp_jra55(reg),yearofmax_q_tmp_jra55(reg)]=max(max(tmp));[~,doyofmax_q_tmp_jra55(reg)]=max(tmp(:,yearofmax_q_tmp_jra55(reg)));

        yearofmax_t_tmp_jra55(reg)=yearofmax_t_tmp_jra55(reg)+firstyear-1;
        yearofmax_tw_tmp_jra55(reg)=yearofmax_tw_tmp_jra55(reg)+firstyear-1;
        yearofmax_q_tmp_jra55(reg)=yearofmax_q_tmp_jra55(reg)+firstyear-1;
    end
    
    if choice==1
        tw_regdecstdanoms_jra55=tw_stdanoms_jra55;t_regdecstdanoms_jra55=t_stdanoms_jra55;q_regdecstdanoms_jra55=q_stdanoms_jra55;
        maxstdanom_tw_jra55=maxstdanom_tw_tmp_jra55;maxstdanom_t_jra55=maxstdanom_t_tmp_jra55;maxstdanom_q_jra55=maxstdanom_q_tmp_jra55;
        yearofmax_tw_jra55=yearofmax_tw_tmp_jra55;yearofmax_t_jra55=yearofmax_t_tmp_jra55;yearofmax_q_jra55=yearofmax_q_tmp_jra55;
        doyofmax_tw_jra55=doyofmax_tw_tmp_jra55;doyofmax_t_jra55=doyofmax_t_tmp_jra55;doyofmax_q_jra55=doyofmax_q_tmp_jra55;
        save(strcat(saveloc,'jra55output_regions.mat'),'maxstdanom_tw_jra55','tw_regdecstdanoms_jra55',...
            'maxstdanom_t_jra55','t_regdecstdanoms_jra55','maxstdanom_q_jra55','q_regdecstdanoms_jra55',...
            'yearofmax_t_jra55','doyofmax_t_jra55','yearofmax_tw_jra55','doyofmax_tw_jra55','yearofmax_q_jra55','doyofmax_q_jra55',...
            'tw_regdecvals_jra55','t_regdecvals_jra55','peakmonth_tw_jra55','peakmonth_t_jra55','peakmonth_q_jra55','-append');
    elseif choice==2
        tw_regdecstdanoms_adj_jra55=tw_stdanoms_jra55;t_regdecstdanoms_adj_jra55=t_stdanoms_jra55;q_regdecstdanoms_adj_jra55=q_stdanoms_jra55;
        maxstdanom_tw_adj_jra55=maxstdanom_tw_tmp_jra55;maxstdanom_t_adj_jra55=maxstdanom_t_tmp_jra55;maxstdanom_q_adj_jra55=maxstdanom_q_tmp_jra55;
        yearofmax_tw_adj_jra55=yearofmax_tw_tmp_jra55;yearofmax_t_adj_jra55=yearofmax_t_tmp_jra55;yearofmax_q_adj_jra55=yearofmax_q_tmp_jra55;
        doyofmax_tw_adj_jra55=doyofmax_tw_tmp_jra55;doyofmax_t_adj_jra55=doyofmax_t_tmp_jra55;doyofmax_q_adj_jra55=doyofmax_q_tmp_jra55;
        save(strcat(saveloc,'jra55output_regions.mat'),'maxstdanom_tw_adj_jra55','tw_regdecstdanoms_adj_jra55',...
            'maxstdanom_t_adj_jra55','t_regdecstdanoms_adj_jra55','maxstdanom_q_adj_jra55','q_regdecstdanoms_adj_jra55',...
            'yearofmax_t_adj_jra55','doyofmax_t_adj_jra55','yearofmax_tw_adj_jra55','doyofmax_tw_adj_jra55',...
            'yearofmax_q_adj_jra55','doyofmax_q_adj_jra55','tw_regdecvals_adj_jra55','t_regdecvals_adj_jra55',...
            'peakmonth_tw_jra55','peakmonth_t_jra55','peakmonth_q_jra55','-append');
    end
end





%Finds regions that have the same top event across:
%1. shifted climatological period (sliding 30-year windows -- considering only those which include the year of the overall top event)
        %default (as also used in T22) is 1981-2010
%2. ERA5 and JRA55 (formerly MERRA2)
%3. (bonus) Tw and T
if completelyconsistenttopevents==1
    %1. Sliding windows (runtime 1 min)
    for loop=1:2 %ERA5, comparison
        if loop==1
            tw_vals=tw_regdecvals_adj_era5;t_vals=t_regdecvals_adj_era5;
            peakmonth_tw=peakmonth_tw_era5;peakmonth_t=peakmonth_t_era5;
        elseif loop==2
            if strcmp(compareagainst,'merra2')
                tw_vals=tw_regdecvals_adj_merra2;t_vals=t_regdecvals_adj_merra2;
                peakmonth_tw=peakmonth_tw_merra2;peakmonth_t=peakmonth_t_merra2;
            elseif strcmp(compareagainst,'jra55')
                tw_vals=tw_regdecvals_adj_jra55;t_vals=t_regdecvals_adj_jra55;
                peakmonth_tw=peakmonth_tw_jra55;peakmonth_t=peakmonth_t_jra55;
            end
        end

        varstocalc={'t';'tw'};
        clear meanwsmax_tw_sliding;
        for starty=1:34 %starting years varying from 1961 to 1994
            endy=starty+29;
            firstclimoyear=starty;lastclimoyear=endy;numclimoyears=30;
            meanwarmseasondailymaxhelper;
    
            meanwsmax_tw_sliding(starty,:)=meanwsmax_tw;meanwsmax_t_sliding(starty,:)=meanwsmax_t;
            stdwsmax_tw_sliding(starty,:)=stdwsmax_tw;stdwsmax_t_sliding(starty,:)=stdwsmax_t;
        end
    
    
        %Then, calculate daily extreme index for all years using Thompson et
        %al. method, with various sliding climatological periods
        clear maxstdanom_t_sliding;clear maxstdanom_tw_sliding;
        clear yearofmax_t_sliding;clear yearofmax_tw_sliding;
        clear doyofmax_t_sliding;clear doyofmax_tw_sliding;
        for starty=1:34 %1961-1994
            clear t_stdanoms_sliding;clear tw_stdanoms_sliding;
            for reg=1:numreg
                t_stdanoms_sliding(reg,:,:)=(squeeze(t_vals(reg,:,:))-squeeze(repmat(meanwsmax_t_sliding(starty,reg),[1 365 numyr])))./...
                    squeeze(repmat(stdwsmax_t_sliding(starty,reg),[1 365 numyr]));
                tw_stdanoms_sliding(reg,:,:)=(squeeze(tw_vals(reg,:,:))-squeeze(repmat(meanwsmax_tw_sliding(starty,reg),[1 365 numyr])))./...
                    squeeze(repmat(stdwsmax_tw_sliding(starty,reg),[1 365 numyr]));
        
                %Get year and DOY of daily-extreme-index max in each region
                tmp=squeeze(t_stdanoms_sliding(reg,:,:));
                    [maxstdanom_t_sliding(starty,reg),yearofmax_t_sliding(starty,reg)]=max(max(tmp));
                    [~,doyofmax_t_sliding(starty,reg)]=max(tmp(:,yearofmax_t_sliding(starty,reg)));
                tmp=squeeze(tw_stdanoms_sliding(reg,:,:));
                    [maxstdanom_tw_sliding(starty,reg),yearofmax_tw_sliding(starty,reg)]=max(max(tmp));
                    [~,doyofmax_tw_sliding(starty,reg)]=max(tmp(:,yearofmax_tw_sliding(starty,reg)));
        
                yearofmax_t_sliding(starty,reg)=yearofmax_t_sliding(starty,reg)+firstyear-1;
                yearofmax_tw_sliding(starty,reg)=yearofmax_tw_sliding(starty,reg)+firstyear-1;
            end
        end

        if loop==1
            maxstdanom_t_era5_sliding=maxstdanom_t_sliding;yearofmax_t_era5_sliding=yearofmax_t_sliding;
            doyofmax_t_era5_sliding=doyofmax_t_sliding;
            maxstdanom_tw_era5_sliding=maxstdanom_tw_sliding;yearofmax_tw_era5_sliding=yearofmax_tw_sliding;
            doyofmax_tw_era5_sliding=doyofmax_tw_sliding;
        elseif loop==2
            if strcmp(compareagainst,'merra2')
                maxstdanom_t_merra2_sliding=maxstdanom_t_sliding;yearofmax_t_merra2_sliding=yearofmax_t_sliding;
                doyofmax_t_merra2_sliding=doyofmax_t_sliding;
                maxstdanom_tw_merra2_sliding=maxstdanom_tw_sliding;yearofmax_tw_merra2_sliding=yearofmax_tw_sliding;
                doyofmax_tw_merra2_sliding=doyofmax_tw_sliding;
            elseif strcmp(compareagainst,'jra55')
                maxstdanom_t_jra55_sliding=maxstdanom_t_sliding;yearofmax_t_jra55_sliding=yearofmax_t_sliding;
                doyofmax_t_jra55_sliding=doyofmax_t_sliding;
                maxstdanom_tw_jra55_sliding=maxstdanom_tw_sliding;yearofmax_tw_jra55_sliding=yearofmax_tw_sliding;
                doyofmax_tw_jra55_sliding=doyofmax_tw_sliding;
            end
        end
    end

    if strcmp(compareagainst,'merra2')
        tw_regdecstdanoms_adj_other=tw_regdecstdanoms_adj_merra2;t_regdecstdanoms_adj_other=t_regdecstdanoms_adj_merra2;
        maxstdanom_tw_other_sliding=maxstdanom_tw_merra2_sliding;maxstdanom_t_other_sliding=maxstdanom_t_merra2_sliding;
        yearofmax_tw_other_sliding=yearofmax_tw_merra2_sliding;yearofmax_t_other_sliding=yearofmax_t_merra2_sliding;
        doyofmax_tw_other_sliding=doyofmax_tw_merra2_sliding;doyofmax_t_other_sliding=doyofmax_t_merra2_sliding;
    elseif strcmp(compareagainst,'jra55')
        tw_regdecstdanoms_adj_other=tw_regdecstdanoms_adj_jra55;t_regdecstdanoms_adj_other=t_regdecstdanoms_adj_jra55;
        maxstdanom_tw_other_sliding=maxstdanom_tw_jra55_sliding;maxstdanom_t_other_sliding=maxstdanom_t_jra55_sliding;
        yearofmax_tw_other_sliding=yearofmax_tw_jra55_sliding;yearofmax_t_other_sliding=yearofmax_t_jra55_sliding;
        doyofmax_tw_other_sliding=doyofmax_tw_jra55_sliding;doyofmax_t_other_sliding=doyofmax_t_jra55_sliding;
    end
    
    passmatrix=NaN.*ones(numreg,11);
    yoflargesttw_other=NaN.*ones(numreg,5);doyoflargesttw_other=NaN.*ones(numreg,5);
    for reg=1:numreg
        alltw_era5=squeeze(tw_regdecstdanoms_adj_era5(reg,:,:));
        [~,yofoverallmaxtw_era5(reg)]=max(max(alltw_era5));[~,doyofoverallmaxtw_era5(reg)]=max(alltw_era5(:,yofoverallmaxtw_era5(reg)));

        alltw_other=squeeze(tw_regdecstdanoms_adj_other(reg,:,:));
        [~,yofoverallmaxtw_other(reg)]=max(max(alltw_other));[~,doyofoverallmaxtw_other(reg)]=max(alltw_other(:,yofoverallmaxtw_other(reg)));

        %Get other-dataset top 5 days for Tw
        tmp=reshape(alltw_other,[365*numyr 1]);myranks=tiedrank(tmp);
        idxoflargest=find(myranks==max(myranks));
            yoflargesttw_other(reg,1)=round2(idxoflargest/365,1,'floor')+1;doyoflargesttw_other(reg,1)=round(rem(idxoflargest/365,1)*365);
        idxof2ndlargest=find(myranks==max(myranks)-1);
            yoflargesttw_other(reg,2)=round2(idxof2ndlargest/365,1,'floor')+1;doyoflargesttw_other(reg,2)=round(rem(idxof2ndlargest/365,1)*365);
        idxof3rdlargest=find(myranks==max(myranks)-2);
            yoflargesttw_other(reg,3)=round2(idxof3rdlargest/365,1,'floor')+1;doyoflargesttw_other(reg,3)=round(rem(idxof3rdlargest/365,1)*365);
        idxof4thlargest=find(myranks==max(myranks)-3);
            yoflargesttw_other(reg,4)=round2(idxof4thlargest/365,1,'floor')+1;doyoflargesttw_other(reg,4)=round(rem(idxof4thlargest/365,1)*365);
        idxof5thlargest=find(myranks==max(myranks)-4);
            yoflargesttw_other(reg,5)=round2(idxof5thlargest/365,1,'floor')+1;doyoflargesttw_other(reg,5)=round(rem(idxof5thlargest/365,1)*365);

        allt_era5=squeeze(t_regdecstdanoms_adj_era5(reg,:,:));
        [~,yofoverallmaxt_era5(reg)]=max(max(allt_era5));[~,doyofoverallmaxt_era5(reg)]=max(allt_era5(:,yofoverallmaxt_era5(reg)));

        allt_other=squeeze(t_regdecstdanoms_adj_other(reg,:,:));
        [~,yofoverallmaxt_other(reg)]=max(max(allt_other));[~,doyofoverallmaxt_other(reg)]=max(allt_other(:,yofoverallmaxt_other(reg)));

        %Get other-dataset top 5 days for T
        tmp=reshape(allt_other,[365*numyr 1]);myranks=tiedrank(tmp);
        idxoflargest=find(myranks==max(myranks));
            yoflargestt_other(reg,1)=round2(idxoflargest/365,1,'floor')+1;doyoflargestt_other(reg,1)=round(rem(idxoflargest/365,1)*365);
        idxof2ndlargest=find(myranks==max(myranks)-1);
            yoflargestt_other(reg,2)=round2(idxof2ndlargest/365,1,'floor')+1;doyoflargestt_other(reg,2)=round(rem(idxof2ndlargest/365,1)*365);
        idxof3rdlargest=find(myranks==max(myranks)-2);
            yoflargestt_other(reg,3)=round2(idxof3rdlargest/365,1,'floor')+1;doyoflargestt_other(reg,3)=round(rem(idxof3rdlargest/365,1)*365);
        idxof4thlargest=find(myranks==max(myranks)-3);
            yoflargestt_other(reg,4)=round2(idxof4thlargest/365,1,'floor')+1;doyoflargestt_other(reg,4)=round(rem(idxof4thlargest/365,1)*365);
        idxof5thlargest=find(myranks==max(myranks)-4);
            yoflargestt_other(reg,5)=round2(idxof5thlargest/365,1,'floor')+1;doyoflargestt_other(reg,5)=round(rem(idxof5thlargest/365,1)*365);

        %A. do sliding climatologies make any difference?
        if sum(yearofmax_tw_era5_sliding(:,reg))==(yofoverallmaxtw_era5(reg)+firstyear-1)*34 &&...
            sum(doyofmax_tw_era5_sliding(:,reg))==(doyofoverallmaxtw_era5(reg))*34 %ERA5 Tw max date and year are robust to shifting
            passmatrix(reg,1)=1;
        end
        if sum(yearofmax_tw_other_sliding(:,reg))==(yofoverallmaxtw_other(reg)+firstyear-1)*34 &&...
            sum(doyofmax_tw_other_sliding(:,reg))==(doyofoverallmaxtw_other(reg))*34 %Other-dataset Tw max date and year are robust to shifting
            passmatrix(reg,2)=1;
        end
        if sum(yearofmax_t_era5_sliding(:,reg))==(yofoverallmaxt_era5(reg)+firstyear-1)*34 &&...
            sum(doyofmax_t_era5_sliding(:,reg))==(doyofoverallmaxt_era5(reg))*34 %ERA5 T max date and year are robust to shifting
            passmatrix(reg,3)=1;
        end
        if sum(yearofmax_t_other_sliding(:,reg))==(yofoverallmaxt_other(reg)+firstyear-1)*34 &&...
            sum(doyofmax_t_other_sliding(:,reg))==(doyofoverallmaxt_other(reg))*34 %Other-dataset T max date and year are robust to shifting
            passmatrix(reg,4)=1;
        end

        %For Tw, do ERA5 and other-dataset max dates exactly match?
        if yofoverallmaxtw_era5(reg)==yofoverallmaxtw_other(reg) && doyofoverallmaxtw_era5(reg)==doyofoverallmaxtw_other(reg);passmatrix(reg,5)=1;end
        %For T, do ERA5 and other-dataset max dates exactly match?
        if yofoverallmaxt_era5(reg)==yofoverallmaxt_other(reg) && doyofoverallmaxt_era5(reg)==doyofoverallmaxt_other(reg);passmatrix(reg,6)=1;end
        %Finally, do ERA5 and other-dataset dates exactly match for both Tw and T (note that T and Tw dates can themselves be different)?
        if passmatrix(reg,5)==1 && passmatrix(reg,6)==1;passmatrix(reg,7)=1;end


        %More relaxed -- allowing ERA5 top day to simply be in the
        %other-dataset top 5 days (and allowing 1 day of offset)
        %Tw
        if (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,1) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,1))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,2) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,2))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,3) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,3))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,4) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,4))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,5) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,5))<=1)
            passmatrix(reg,8)=1;
        end
        %T
        if (yofoverallmaxt_era5(reg)==yoflargestt_other(reg,1) && abs(doyofoverallmaxt_era5(reg)-doyoflargestt_other(reg,1))<=1) ||...
                (yofoverallmaxt_era5(reg)==yoflargestt_other(reg,2) && abs(doyofoverallmaxt_era5(reg)-doyoflargestt_other(reg,2))<=1) ||...
                (yofoverallmaxt_era5(reg)==yoflargestt_other(reg,3) && abs(doyofoverallmaxt_era5(reg)-doyoflargestt_other(reg,3))<=1) ||...
                (yofoverallmaxt_era5(reg)==yoflargestt_other(reg,4) && abs(doyofoverallmaxt_era5(reg)-doyoflargestt_other(reg,4))<=1) ||...
                (yofoverallmaxt_era5(reg)==yoflargestt_other(reg,5) && abs(doyofoverallmaxt_era5(reg)-doyoflargestt_other(reg,5))<=1)
            passmatrix(reg,9)=1;
        end
        %Do ERA5 and other-dataset dates match for both Tw and T (again, with this looser definition)?
        if passmatrix(reg,8)==1 && passmatrix(reg,9)==1;passmatrix(reg,10)=1;end


        

        %Finally -- is ERA5 top Tw day among the top 5 other-dataset days in both
        %T and Tw (i.e. is the same event top by both metrics)?
        if ((yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,1) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,1))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,2) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,2))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,3) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,3))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,4) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,4))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargesttw_other(reg,5) && abs(doyofoverallmaxtw_era5(reg)-doyoflargesttw_other(reg,5))<=1)) &&...
                ((yofoverallmaxtw_era5(reg)==yoflargestt_other(reg,1) && abs(doyofoverallmaxtw_era5(reg)-doyoflargestt_other(reg,1))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargestt_other(reg,2) && abs(doyofoverallmaxtw_era5(reg)-doyoflargestt_other(reg,2))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargestt_other(reg,3) && abs(doyofoverallmaxtw_era5(reg)-doyoflargestt_other(reg,3))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargestt_other(reg,4) && abs(doyofoverallmaxtw_era5(reg)-doyoflargestt_other(reg,4))<=1) ||...
                (yofoverallmaxtw_era5(reg)==yoflargestt_other(reg,5) && abs(doyofoverallmaxtw_era5(reg)-doyoflargestt_other(reg,5))<=1))
            passmatrix(reg,11)=1;
        end

        %Where date of top Tw exactly matches between ERA5 and other dataset *AND* is >3 std anoms in both datasets
        toptwexact=NaN.*ones(numreg,1);
        for i=1:237;if passmatrix(i,5)==1;if tw_regdecstdanoms_adj_era5(i,doyofoverallmaxtw_era5(i),yofoverallmaxtw_era5(i))>=3;...
                        toptwexact(i)=1;end;end;end
        %Ditto for T
        toptexact=NaN.*ones(numreg,1);
        for i=1:237;if passmatrix(i,6)==1;if t_regdecstdanoms_adj_era5(i,doyofoverallmaxt_era5(i),yofoverallmaxt_era5(i))>=3;...
                        toptexact(i)=1;end;end;end
    end


    %Create map
    consistenttw=regs_era5sz; %to initialize
    consistentt=regs_era5sz; %to initialize
    sameeventtandtw=regs_era5sz; %to initialize
    for reg=1:numreg
        consistenttw(regs_era5sz==reg)=passmatrix(reg,8);
        consistentt(regs_era5sz==reg)=passmatrix(reg,9);
        sameeventtandtw(regs_era5sz==reg)=passmatrix(reg,11);
    end

    figure(777);clf;hold on;figname='completeconsistency';curpart=1;highqualityfiguresetup;

    data={double(lat2d_era5);double(lon2d_era5);consistenttw};cmap=flipud(colormaps('orangewhite','2','not'));
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0.5;'underlaycaxismax';1.5;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'omitfirstsubplotcolorbar';1;'nonewfig';1;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[0.1 0.7 0.8 0.29]);
    t=text(-0.02,0.49,'a)','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);

    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);consistentt};cmap=flipud(colormaps('orangewhite','2','not'));
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0.5;'underlaycaxismax';1.5;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'omitfirstsubplotcolorbar';1;'nonewfig';1;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.37 0.8 0.29]);
    t=text(-0.02,0.49,'b)','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);

    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);sameeventtandtw};cmap=flipud(colormaps('orangewhite','2','not'));
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0.5;'underlaycaxismax';1.5;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'omitfirstsubplotcolorbar';1;'nonewfig';1;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.04 0.8 0.29]);
    t=text(-0.02,0.49,'c)','units','normalized');set(t,'fontweight','bold','fontname','arial','fontsize',12);

    curpart=2;highqualityfiguresetup;
end




%Was previously based on ERA5/JRA55 mean; now based just on ERA5 for consistency with figures
if preparetable2==1
    retpers=currecretpers_tw;
    tmparr=[[1:numreg]' retpers];
    tmparr=sortrows(tmparr,2,'ascend');

    %vals=(tw_regdecvals_adj_era5+tw_regdecvals_adj_jra55)./2;
    vals=tw_regdecvals_adj_era5;

    for i=1:size(tmparr,1)
        if tmparr(i,2)<=100 %<100-year return periods using GMST-adjusted ERA5/JRA55-mean data
            %Get date and value of record
            thisreg=tmparr(i,1);
            [tmparr(i,3),tmparr(i,4)]=max(max(vals(thisreg,:,:)));
            [~,tmparr(i,5)]=max(vals(thisreg,:,tmparr(i,4)));
        end
    end
    tmparr(:,4)=tmparr(:,4)+firstyear-1;
end


if morestats==1
    %Stats for regions' extreme events
    %Recorded in resultscomparison.xlsx
    c1=0;c2=0;c3=0;c4=0;c5=0;c6=0;c7=0;c8=0;c9=0;clear A;clear B;clear C;clear D;clear E;clear F;clear G;clear H;clear I;
    clear tdatehere;clear twdatehere;clear twdatehere_unadj;clear twdatehere_adj;clear tdatehere_adj;
    clear Fdate;clear Gdate;clear Hdate;clear Idate;
    for reg=1:numreg
        %For results_comparison Tw table
        regdata=squeeze(tw_regdecstdanoms_adj_era5(reg,:,:));
        [twmax,yearofmax]=max(max(regdata));
        [~,doyofmax]=max(regdata(:,yearofmax));
        if twmax>=4
            c1=c1+1;
            reghere=reg;
            yearhere=yearofmax+firstyear-1;
            doyhere=doyofmax;
            twdatehere{c1}=DOYtoDate(doyhere,yearhere);
            maghere=round2(twmax,0.01);
            valhere=round2(tw_regdecvals_adj_era5(reg,doyofmax,yearofmax),0.01);
            jra55valhere=round2(tw_regdecvals_adj_jra55(reg,doyofmax,yearofmax),0.01);
            merra2valhere=round2(tw_regdecvals_adj_merra2(reg,doyofmax,yearofmax),0.01);

            A(c1,1)=reghere;
            A(c1,2)=yearhere;
            A(c1,3)=doyhere;
            A(c1,4)=maghere;
            A(c1,5)=valhere;
            A(c1,6)=jra55valhere;
            A(c1,7)=merra2valhere;
        end 

        %For results_comparison T table
        regdata=squeeze(t_regdecstdanoms_adj_era5(reg,:,:));
        [tmax,yearofmax]=max(max(regdata));
        [~,doyofmax]=max(regdata(:,yearofmax));
        if tmax>=4
            c2=c2+1;
            reghere=reg;
            yearhere=yearofmax+firstyear-1;
            doyhere=doyofmax;
            tdatehere{c2}=DOYtoDate(doyhere,yearhere);
            maghere=round2(tmax,0.01);
            valhere=round2(t_regdecvals_adj_era5(reg,doyofmax,yearofmax),0.01);
            jra55valhere=round2(t_regdecvals_adj_jra55(reg,doyofmax,yearofmax),0.01);
            merra2valhere=round2(t_regdecvals_adj_merra2(reg,doyofmax,yearofmax),0.01);

            B(c2,1)=reghere;
            B(c2,2)=yearhere;
            B(c2,3)=doyhere;
            B(c2,4)=maghere;
            B(c2,5)=valhere;
            B(c2,6)=jra55valhere;
            B(c2,7)=merra2valhere;
        end 

        %For Table S3
        regdata_era5=squeeze(tw_regdecstdanoms_era5(reg,:,:));
        [twmax_era5,yearofmax_era5]=max(max(regdata_era5));
        [~,doyofmax_era5]=max(regdata_era5(:,yearofmax_era5));

        regdata_jra55=squeeze(tw_regdecstdanoms_jra55(reg,:,:));
        [twmax_jra55,yearofmax_jra55]=max(max(regdata_jra55));
        [~,doyofmax_jra55]=max(regdata_jra55(:,yearofmax_jra55));

        if (twmax_era5+twmax_jra55)/2>=3
            c3=c3+1;
            reghere=reg;
            yearhere=yearofmax_era5+firstyear-1;
            doyhere=doyofmax_era5;
            twdatehere_unadj{c3}=DOYtoDate(doyhere,yearhere);
            maghere=round2(twmax_era5,0.01);
            valhere=round2(tw_regdecvals_era5(reg,doyofmax_era5,yearofmax_era5),0.01);
            jra55stdanomhere=round2(tw_regdecstdanoms_jra55(reg,doyofmax_jra55,yearofmax_jra55),0.01);
            jra55valhere=round2(tw_regdecvals_jra55(reg,doyofmax_jra55,yearofmax_jra55),0.01);
            %merra2valhere=round2(tw_regdecvals_merra2(reg,doyofmax,yearofmax),0.01);

            C(c3,1)=reghere;
            C(c3,2)=yearhere;
            C(c3,3)=doyhere;
            C(c3,4)=maghere;
            C(c3,5)=valhere;
            C(c3,6)=jra55stdanomhere;
            C(c3,7)=jra55valhere;
            %C(c3,8)=merra2valhere;
        end 

        %Most extreme Tw events -- GMST-adjusted
        %Criterion: Std anom >= 3.5 in ERA5/JRA55 mean
        %Table S1
        anoms_era5=tw_regdecstdanoms_adj_era5;anoms_jra55=tw_regdecstdanoms_adj_jra55;
        vals_era5=tw_regdecvals_adj_era5;vals_jra55=tw_regdecvals_adj_jra55;

        anoms=anoms_era5;

        regdata=squeeze(anoms(reg,:,:));
        [twmax,yearofmax]=max(max(regdata));
        [~,doyofmax]=max(regdata(:,yearofmax));

        if twmax>=3.5
            c4=c4+1;
            reghere=reg;
            yearhere=yearofmax+firstyear-1;
            doyhere=doyofmax;
            twdatehere_adj{c4}=DOYtoDate(doyhere,yearhere);
            maghere=round2(twmax,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax,yearofmax),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax,yearofmax),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax,yearofmax),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax,yearofmax),0.1);

            D(c4,1)=reghere;
            D(c4,2)=yearhere;
            D(c4,3)=DOYtoMonth(doyhere,yearhere);
            D(c4,4)=DOYtoDOM(doyhere,yearhere);
            D(c4,5)=era5maghere;
            D(c4,6)=era5valhere;
            %D(c4,8)=jra55maghere;
            D(c4,7)=jra55valhere;
        end 

        %Most extreme T events -- GMST-adjusted
        %Criterion: Std anom >= 3.5 in ERA5/JRA55 mean
        anoms_era5=t_regdecstdanoms_adj_era5;anoms_jra55=t_regdecstdanoms_adj_jra55;
        vals_era5=t_regdecvals_adj_era5;vals_jra55=t_regdecvals_adj_jra55;

        combined=(anoms_era5+anoms_jra55)./2;

        regdata_combo=squeeze(combined(reg,:,:));
        [tmax_combo,yearofmax_combo]=max(max(regdata_combo));
        [~,doyofmax_combo]=max(regdata_combo(:,yearofmax_combo));

        if tmax_combo>=3.5
            c5=c5+1;
            reghere=reg;
            yearhere=yearofmax_combo+firstyear-1;
            doyhere=doyofmax_combo;
            tdatehere_adj{c5}=DOYtoDate(doyhere,yearhere);
            maghere=round2(tmax_combo,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax_combo,yearofmax_combo),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax_combo,yearofmax_combo),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax_combo,yearofmax_combo),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax_combo,yearofmax_combo),0.1);

            E(c5,1)=reghere;
            E(c5,2)=yearhere;
            E(c5,3)=doyhere;
            E(c5,4)=maghere;
            E(c5,5)=era5maghere;
            E(c5,6)=era5valhere;
            E(c5,7)=jra55maghere;
            E(c5,8)=jra55valhere;
        end 

        %Most extreme Tw events -- NOT GMST-adjusted
        %Criterion: Std anom >= 3.5 in ERA5/JRA55 mean
        %Table S2
        anoms_era5=tw_regdecstdanoms_era5;anoms_jra55=tw_regdecstdanoms_jra55;
        vals_era5=tw_regdecvals_era5;vals_jra55=tw_regdecvals_jra55;
        combined=anoms_era5;
        regdata_combo=squeeze(combined(reg,:,:));
        [twmax_combo,yearofmax_combo]=max(max(regdata_combo));
        [~,doyofmax_combo]=max(regdata_combo(:,yearofmax_combo));

        if twmax_combo>=4
            c6=c6+1;
            reghere=reg;
            yearhere=yearofmax_combo+firstyear-1;
            doyhere=doyofmax_combo;
            Fdate{c6}=DOYtoDate(doyhere,yearhere);
            maghere=round2(twmax_combo,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax_combo,yearofmax_combo),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax_combo,yearofmax_combo),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax_combo,yearofmax_combo),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax_combo,yearofmax_combo),0.1);

            F(c6,1)=reghere;
            F(c6,2)=yearhere;
            F(c6,3)=DOYtoMonth(doyhere,yearhere);
            F(c6,4)=DOYtoDOM(doyhere,yearhere);
            F(c6,5)=era5maghere;
            F(c6,6)=era5valhere;
            F(c6,7)=jra55valhere;
        end 

        %Most extreme T events -- NOT GMST-adjusted
        %Criterion: Std anom >= 3.5 in ERA5/JRA55 mean
        anoms_era5=t_regdecstdanoms_era5;anoms_jra55=t_regdecstdanoms_jra55;
        vals_era5=t_regdecvals_era5;vals_jra55=t_regdecvals_jra55;
        combined=(anoms_era5+anoms_jra55)./2;
        regdata_combo=squeeze(combined(reg,:,:));
        [tmax_combo,yearofmax_combo]=max(max(regdata_combo));
        [~,doyofmax_combo]=max(regdata_combo(:,yearofmax_combo));

        if tmax_combo>=3.5
            c7=c7+1;
            reghere=reg;
            yearhere=yearofmax_combo+firstyear-1;
            doyhere=doyofmax_combo;
            Gdate{c7}=DOYtoDate(doyhere,yearhere);
            maghere=round2(tmax_combo,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax_combo,yearofmax_combo),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax_combo,yearofmax_combo),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax_combo,yearofmax_combo),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax_combo,yearofmax_combo),0.1);

            G(c7,1)=reghere;
            G(c7,2)=yearhere;
            G(c7,3)=doyhere;
            G(c7,4)=maghere;
            G(c7,5)=era5maghere;
            G(c7,6)=era5valhere;
            G(c7,7)=jra55maghere;
            G(c7,8)=jra55valhere;
        end 

        %Most extreme Tw events -- NOT GMST-adjusted
        %Criterion: Tw value in ERA5/JRA55 mean >=27C
        %Table S3
        anoms_era5=tw_regdecstdanoms_era5;anoms_jra55=tw_regdecstdanoms_jra55;
        vals_era5=tw_regdecvals_era5;vals_jra55=tw_regdecvals_jra55;
        anoms=vals_era5;
        regdata_combo=squeeze(anoms(reg,:,:));
        [twmax_combo,yearofmax_combo]=max(max(regdata_combo));
        [~,doyofmax_combo]=max(regdata_combo(:,yearofmax_combo));

        if twmax_combo>=27
            c8=c8+1;
            reghere=reg;
            yearhere=yearofmax_combo+firstyear-1;
            doyhere=doyofmax_combo;
            Hdate{c8}=DOYtoDate(doyhere,yearhere);
            valhere=round2(twmax_combo,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax_combo,yearofmax_combo),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax_combo,yearofmax_combo),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax_combo,yearofmax_combo),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax_combo,yearofmax_combo),0.1);

            H(c8,1)=reghere;
            H(c8,2)=yearhere;
            H(c8,3)=DOYtoMonth(doyhere,yearhere);
            H(c8,4)=DOYtoDOM(doyhere,yearhere);
            H(c8,5)=era5valhere;
            H(c8,6)=era5maghere;
            H(c8,7)=jra55valhere;
        end 

        %Most extreme T events -- NOT GMST-adjusted
        %Criterion: T value in ERA5/JRA55 mean >=45C
        anoms_era5=t_regdecstdanoms_era5;anoms_jra55=t_regdecstdanoms_jra55;
        vals_era5=t_regdecvals_era5;vals_jra55=t_regdecvals_jra55;
        combined=(vals_era5+vals_jra55)./2;
        regdata_combo=squeeze(combined(reg,:,:));
        [tmax_combo,yearofmax_combo]=max(max(regdata_combo));
        [~,doyofmax_combo]=max(regdata_combo(:,yearofmax_combo));

        if tmax_combo>=45
            c9=c9+1;
            reghere=reg;
            yearhere=yearofmax_combo+firstyear-1;
            doyhere=doyofmax_combo;
            Idate{c9}=DOYtoDate(doyhere,yearhere);
            valhere=round2(tmax_combo,0.01);

            era5maghere=round2(anoms_era5(reg,doyofmax_combo,yearofmax_combo),0.01);
            era5valhere=round2(vals_era5(reg,doyofmax_combo,yearofmax_combo),0.1);
            jra55maghere=round2(anoms_jra55(reg,doyofmax_combo,yearofmax_combo),0.01);
            jra55valhere=round2(vals_jra55(reg,doyofmax_combo,yearofmax_combo),0.1);

            I(c9,1)=reghere;
            I(c9,2)=yearhere;
            I(c9,3)=doyhere;
            I(c9,4)=valhere;
            I(c9,5)=era5maghere;
            I(c9,6)=era5valhere;
            I(c9,7)=jra55maghere;
            I(c9,8)=jra55valhere;
        end 
    end

    %Most extreme Tw events -- NOT GMST-adjusted
    %Criterion: Top 10 Tw value when averaged across ERA5 and JRA55
    %In this listing, regions can repeat
    anoms_era5=tw_regdecstdanoms_era5;anoms_jra55=tw_regdecstdanoms_jra55;
    vals_era5=tw_regdecvals_era5;vals_jra55=tw_regdecvals_jra55;
    combined=(vals_era5+vals_jra55)./2;
    combined1d=reshape(combined,[size(combined,1)*size(combined,2)*size(combined,3) 1]);
    [vals,inds]=maxk(combined1d,10);
    reglist=zeros(10,1);doylist=zeros(10,1);yearlist=zeros(10,1);
    for i=1:10
        [reglist(i),doylist(i),yearlist(i)]=ind2sub(size(combined),inds(i));
    end
    yearlist=yearlist+firstyear-1;

    %Not so interesting -- all Central East China!
    

     %Most extreme T events -- NOT GMST-adjusted
    %Criterion: Top 10 T value when averaged across ERA5 and JRA55
    %In this listing, regions can repeat
    anoms_era5=t_regdecstdanoms_era5;anoms_jra55=t_regdecstdanoms_jra55;
    vals_era5=t_regdecvals_era5;vals_jra55=t_regdecvals_jra55;
    combined=(vals_era5+vals_jra55)./2;
    combined1d=reshape(combined,[size(combined,1)*size(combined,2)*size(combined,3) 1]);
    [vals,inds]=maxk(combined1d,10);
    reglist=zeros(10,1);doylist=zeros(10,1);yearlist=zeros(10,1);
    for i=1:10
        [reglist(i),doylist(i),yearlist(i)]=ind2sub(size(combined),inds(i));
    end
    yearlist=yearlist+firstyear-1;

    %For Tw -- not so interesting -- all Central East China!
    %For T -- all E Saudi Arabia

    C=sortrows(C,4,'descend');D=sortrows(D,5,'descend');E=sortrows(E,4,'descend');
    F=sortrows(F,5,'descend');G=sortrows(G,4,'descend');H=sortrows(H,5,'descend');I=sortrows(I,4,'descend');
end

if assessdiurnaleffects==1
    y=1961;
    if y<=1999;era5dir=era5dir_1;else;era5dir=era5dir_2;end

    lon2d_reduc=lon2d_era5(2:2:end,2:2:end);lat2d_reduc=lat2d_era5(2:2:end,2:2:end);
    regs_reduc=regs_era5sz(2:2:end,2:2:end);
    
    z0=interp2(lon2d_192x288,lat2d_192x288,z_192x288,lon2d_reduc,lat2d_reduc);
    z0(:,719:720)=z0(:,717:718); %it's nearly all ocean at this lon anyway, so it makes little difference
    psfc_tmp=pressurefromelev(z0);
    psfc_era5_3hr=repmat(psfc_tmp,[1 1 365*8]);


    %Read in T and Td (1 min 30 sec)
    t=ncread(strcat(era5dir,'t2m_',num2str(y),'.nc'),'t2m')-273.15; %C
    td=ncread(strcat(era5dir,'td2m_',num2str(y),'.nc'),'d2m')-273.15; %C
    t=t(2:2:end,2:2:end,:);td=td(2:2:end,2:2:end,:);

    %Flip and recenter
    trev=NaN.*ones(size(t,2),size(t,1),size(t,3));
    for timestep=1:365*8;trev(:,:,timestep)=recenter(squeeze(t(:,:,timestep))');end
    tdailymax=NaN.*ones(size(t,2),size(t,1),365);
    for day=1:365;tdailymax(:,:,day)=recenter(squeeze(max(t(:,:,day*8-7:day*8),[],3))');end;clear t;

    tdrev=NaN.*ones(size(td,2),size(td,1),size(td,3));
    for timestep=1:365*8;tdrev(:,:,timestep)=recenter(squeeze(td(:,:,timestep))');end
    tddailymean=NaN.*ones(size(td,2),size(td,1),365);
    for day=1:365;tddailymean(:,:,day)=recenter(squeeze(mean(td(:,:,day*8-7:day*8),3))');end;clear td;

    %Prepare to do lat weighting and regionalization
    tregmeans=cell(numreg,1);tdregmeans=cell(numreg,1);psfcregmeans=cell(numreg,1);c=zeros(numreg,1); %237 regions
    tregmeans_daily=cell(numreg,1);tdregmeans_daily=cell(numreg,1);psfcregmeans_daily=cell(numreg,1);
    latweights=cell(numreg,1);
    for i=1:size(lat2d_reduc,1)
        for j=1:size(lat2d_reduc,2)
            thisreg=round(regs_reduc(i,j));
            if ~isnan(thisreg)
                c(thisreg)=c(thisreg)+1;
                tregmeans{thisreg}(c(thisreg),:)=trev(i,j,:);
                tdregmeans{thisreg}(c(thisreg),:)=tdrev(i,j,:);
                psfcregmeans{thisreg}(c(thisreg),:)=psfc_era5_3hr(i,j,:);

                tregmeans_daily{thisreg}(c(thisreg),:)=tdailymax(i,j,:);
                tdregmeans_daily{thisreg}(c(thisreg),:)=tddailymean(i,j,:);
                psfcregmeans_daily{thisreg}(c(thisreg),:)=psfc_era5_3hr(i,j,1:365);

                latweights{thisreg}(c(thisreg))=cos(deg2rad(lat2d_reduc(i,j)));
            end
        end
    end
    clear trev;clear tdrev;clear psfcera5_3hr;


    %Spatial mean by region
    tregmean=NaN.*ones(numreg,365*8);tdregmean=NaN.*ones(numreg,365*8);psfcregmean=NaN.*ones(numreg,365*8);
    tregmean_daily=NaN.*ones(numreg,365);tdregmean_daily=NaN.*ones(numreg,365);psfcregmean_daily=NaN.*ones(numreg,365);
    for r=1:numreg
        if latwgt_era5==0 %no lat weighting (original)
        else %latitudinally weighted
            weights=repmat(latweights{r}',[1 365*8]);weights_daily=repmat(latweights{r}',[1 365]);
            for timestep=1:365*8
                tregmean(r,timestep)=sum(weights(:,timestep).*tregmeans{r}(:,timestep))./sum(weights(:,timestep));
                tdregmean(r,timestep)=sum(weights(:,timestep).*tdregmeans{r}(:,timestep))./sum(weights(:,timestep));
                psfcregmean(r,timestep)=sum(weights(:,timestep).*psfcregmeans{r}(:,timestep))./sum(weights(:,timestep));
            end
            for timestep=1:365
                tregmean_daily(r,timestep)=sum(weights_daily(:,timestep).*tregmeans_daily{r}(:,timestep))./sum(weights_daily(:,timestep));
                tdregmean_daily(r,timestep)=sum(weights_daily(:,timestep).*tdregmeans_daily{r}(:,timestep))./sum(weights_daily(:,timestep));
                psfcregmean_daily(r,timestep)=sum(weights_daily(:,timestep).*psfcregmeans_daily{r}(:,timestep))./sum(weights_daily(:,timestep));
            end
        end
    end
    clear tregmeans;clear tdregmeans;clear psfcregmeans;

    arrp=psfcregmean.*100;arrt=tregmean;arrq=calcqfromTd_dynamicP(tdregmean,arrp)./1000;
    tw_era5=calcwbt_daviesjones(arrt,arrp,arrq);
    clear arrp;clear arrt;clear arrq;

    arrp=psfcregmean_daily.*100;arrt=tregmean_daily;arrq=calcqfromTd_dynamicP(tdregmean_daily,arrp)./1000;
    tw_era5_daily=calcwbt_daviesjones(arrt,arrp,arrq);
    clear arrp;clear arrt;clear arrq;

    save(strcat(saveloc,'era5tw_regions_hourlycalc.mat'),'tw_era5','tw_era5_daily');

    

    %Compare against normal method
    %f=load(strcat(saveloc,'mylatestera5output_regions.mat'));tw_regdecvals_era5=f.tw_regdecvals_era5;
    %origway=squeeze(tw_regdecvals_era5(:,:,1));
    origway=tw_era5_daily;


    %Get daily maxes from more-precise way
    morepreciseway=NaN.*ones(size(origway));
    for endhr=8:8:365*8
        morepreciseway(:,endhr/8)=max(tw_era5(:,endhr-7:endhr),[],2);
    end

    %Difference
    mydiff=origway-morepreciseway;

    %Warm-season difference in the median and in the 95th pctile
    for reg=1:numreg
        if peakmonth_tw_era5(reg)==12
            myvals1=origway(reg,335:365)';myvals2=origway(reg,1:59)';origvals=[myvals1;myvals2];
            myvals1=morepreciseway(reg,335:365)';myvals2=morepreciseway(reg,1:59)';moreprecisevals=[myvals1;myvals2];
        elseif peakmonth_tw_era5(reg)==1
            myvals1=origway(reg,305:365)';myvals2=origway(reg,1:31)';origvals=[myvals1;myvals2];
            myvals1=morepreciseway(reg,305:365)';myvals2=morepreciseway(reg,1:31)';moreprecisevals=[myvals1;myvals2];
        else
            month1=peakmonth_tw_era5(reg)-1;month2=peakmonth_tw_era5(reg)+1;
            origvals=origway(reg,monthstarts(month1):monthstops(month2))';
            moreprecisevals=morepreciseway(reg,monthstarts(month1):monthstops(month2))';
        end
        wsmedian_orig=quantile(origvals,0.5);wsp95_orig=quantile(origvals,0.95);
        wsmedian_moreprecise=quantile(moreprecisevals,0.5);wsp95_moreprecise=quantile(moreprecisevals,0.95);

        twdiff_wsmedian(reg)=wsmedian_orig-wsmedian_moreprecise;
        twdiff_wsp95(reg)=wsp95_orig-wsp95_moreprecise;
    end

    %Plot map
    mediandiff=regs_era5sz; %to initialize
    p95diff=regs_era5sz; %to initialize
    for reg=1:numreg
        mediandiff(regs_era5sz==reg)=twdiff_wsmedian(reg);
        p95diff(regs_era5sz==reg)=twdiff_wsp95(reg);
    end

    figure(200);clf;hold on;curpart=1;highqualityfiguresetup;figname='figsXX';
    data={double(lat2d_era5);double(lon2d_era5);mediandiff};cmap=colormaps('t','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';-1;'underlaycaxismax';1;'underlaystepsize';0.2;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'omitfirstsubplotcolorbar';0;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[0.1 0.5 0.8 0.45]);

    subplot(10,10,100);
    data={double(lat2d_era5);double(lon2d_era5);p95diff};cmap=colormaps('t','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';-1;'underlaycaxismax';1;'underlaystepsize';0.2;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'omitfirstsubplotcolorbar';0;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.02 0.8 0.45]);

    curpart=2;highqualityfiguresetup;
end


%Combines older Fig 2 ERA5 Tw return periods with older Fig 5 GCM outside-GEV-fit fractions
if createfig2==1
    %GCM info
    outsidegev_canesm5=impossible_pct_canesm5;
    outsidegev_miroc6=impossible_pct_miroc6;
    outsidegev_mpige=impossible_pct_mpige;
    outsidegev_mean=(outsidegev_canesm5+outsidegev_miroc6+outsidegev_mpige)./3;
    
    %ERA5 info
    xarr=[100./currecretpers_tw(1:170);100./currecretpers_tw(192:end)]; %exclude Antarctica

    yarr=[outsidegev_mean(1:170);outsidegev_mean(192:end)];
    %scatter(xarr,yarr);

    %For actual probabilities (panel a): cutoffs mark 0.1% and 1%, and 25% and 40% of ensemble members

    
    x90=quantile(xarr,0.9);x67=quantile(xarr,0.666);
    y90=quantile(yarr,0.9);y67=quantile(yarr,0.666);

    xcateg_a=zeros(numreg_exclant,1);ycateg_a=zeros(numreg_exclant,1);
    xcateg_b=zeros(numreg_exclant,1);ycateg_b=zeros(numreg_exclant,1);
    for reg=1:numreg_exclant
        %For panel a: cutoffs are at fixed thresholds
        if xarr(reg)>=1
            xcateg_a(reg)=0.5; %i.e. will be #3 (top)
        elseif xarr(reg)>=0.1
            xcateg_a(reg)=1.5;
        else
            xcateg_a(reg)=2.5;
        end

        if yarr(reg)>=40
            ycateg_a(reg)=0.5; %i.e. will be #3 (top)
        elseif yarr(reg)>=25
            ycateg_a(reg)=1.5;
        else
            ycateg_a(reg)=2.5;
        end


        %For regional relative ranking (panel b): cutoffs marking top 10% and 33% of regions for each measure
        if xarr(reg)>=x90
            xcateg_b(reg)=0.5; %i.e. will be #3 (top)
        elseif xarr(reg)>=x67
            xcateg_b(reg)=1.5;
        else
            xcateg_b(reg)=2.5;
        end

        if yarr(reg)>=y90
            ycateg_b(reg)=0.5; %i.e. will be #3 (top)
        elseif yarr(reg)>=y67
            ycateg_b(reg)=1.5;
        else
            ycateg_b(reg)=2.5;
        end
    end

    topleftcolor=colors('medium green');bottomrightcolor=colors('medium red');toprightcolor=colors('purple');
    incrcontrast=1;palefactor=0.4;numcolors=3;
    colorshading=colorfield2d(topleftcolor,bottomrightcolor,toprightcolor,numcolors,incrcontrast,palefactor);


    %Set up to plot
    axmin=0;axmax=3;
    regcolor_a=ones(numreg_exclant,3);regcolor_b=ones(numreg_exclant,3);
    for reg=1:numreg_exclant
        %For absolute thresholds
        xrelpos_a=round(numcolors*(xcateg_a(reg)-axmin)/(axmax-axmin));
        yrelpos_a=round(numcolors*(ycateg_a(reg)-axmin)/(axmax-axmin));
        if xrelpos_a>numcolors;xrelpos_a=numcolors;end;if yrelpos_a>numcolors;yrelpos_a=numcolors;end
        if xrelpos_a<1;xrelpos_a=1;end;if yrelpos_a<1;yrelpos_a=1;end
        
        if ~isnan(xrelpos_a) && ~isnan(yrelpos_a)
            regcolor_a(reg,:)=colorshading(yrelpos_a,(numcolors+1)-xrelpos_a,:); %x and y axes were initially flipped the wrong way
        end

        %For regional relative ranking
        xrelpos_b=round(numcolors*(xcateg_b(reg)-axmin)/(axmax-axmin));
        yrelpos_b=round(numcolors*(ycateg_b(reg)-axmin)/(axmax-axmin));
        if xrelpos_b>numcolors;xrelpos_b=numcolors;end;if yrelpos_b>numcolors;yrelpos_b=numcolors;end
        if xrelpos_b<1;xrelpos_b=1;end;if yrelpos_b<1;yrelpos_b=1;end
        
        if ~isnan(xrelpos_b) && ~isnan(yrelpos_b)
            regcolor_b(reg,:)=colorshading(yrelpos_b,(numcolors+1)-xrelpos_b,:);
        end
    end
    

    figure(900);clf;hold on;figname='fig2_latest';curpart=1;highqualityfiguresetup;
    lefts=[0.27;0.27];bottoms=[0.51;0.01];wwidth=0.70;hheight=0.48;
    numpanels=2; %orig was 1 (current panel b only)
    if numpanels==1
        i1=-0.17;i2=-0.123;i3=-0.06;i4=0.66;i5=0.31;i6=-0.05;i7=0.05;i8=0.41;i9=0.76;i10=-0.05;i11=-0.047;i12=-0.17;bottompos=0.02;leftpos=0;
    elseif numpanels==2
        i1=-0.14;i2=0.03;i3=-0.04;i4=[0.765;0.765];i5=[0.385;0.385];i6=[0.005;0.005];i7=[0.03;0.05];
        i8=[0.41;0.412];i9=[0.794;0.774];i10=-0.05;i11=0.06;i12=-0.16;bottompos=0.03;leftpos=0.03;
        legendbottom=[0.63;0.13];
        subpaneltitles={'Record-Breaking Events: Likelihoods';'Record-Breaking Events: Regional Rankings'};
    end

    for loop=1:2
        if loop==1
            mycmap=regcolor_a;caxmin=0.5;caxmax=numreg_exclant+0.5;step=1;
            xlabel1='>1%';xlabel2='>0.1%';xlabel3='<0.1%';ylabel1='>40%';ylabel2='>25%';ylabel3='<25%';
        elseif loop==2
            mycmap=regcolor_b;caxmin=0.5;caxmax=numreg_exclant+0.5;step=1;
            xlabel1='>p90';xlabel2='>p67';xlabel3='<p67';ylabel1='>p90';ylabel2='>p67';ylabel3='<p67';
        end
        if loop==2;axes('position',[0.01 0.01 0.01 0.01]);end
        data={double(lat2d_era5);double(lon2d_era5);regs_era5sz};cmap=mycmap;
        vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxmin;'underlaycaxismax';caxmax;'underlaystepsize';step;'underlaycolormap';cmap;
            'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
            'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'omitfirstsubplotcolorbar';1;...
            'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.65};
        datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
        set(gca,'position',[lefts(loop) bottoms(loop) wwidth hheight]);
        title(subpaneltitles{loop},'fontsize',14,'fontweight','bold','fontname','arial');

        if numpanels==2
            t=text(0.02,0.48,splabels{loop},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');
        end
    
        %Make the legend
        if numpanels==1
            axpos=[0.13 0.10 0.1 0.14];
        elseif numpanels==2
            axpos=[0.09 legendbottom(loop) 0.16 0.182];
        end
        axes('position',[0.01 0.01 0.01 0.01]);image(gca,colorshading);
        set(gca,'position',axpos,'xtick',[],'xticklabel',{},'ytick',[],'yticklabel',{});
    
        t=text(i1,i2,{'   Ens Members with ',' Outside-GEV Events'},'units','normalized');
            set(t,'fontsize',10,'fontweight','bold','fontname','arial','rotation',90);
        toppos=0.48;yposdist=toppos-bottompos;
        rightpos=0.48;xposdist=rightpos-leftpos;
    
        t=text(i3,bottompos+i4(loop)*yposdist,ylabel1,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial','rotation',90);
        t=text(i3,bottompos+i5(loop)*yposdist,ylabel2,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial','rotation',90);
        t=text(i3,bottompos+i6(loop)*yposdist,ylabel3,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial','rotation',90);
    
        t=text(leftpos+i7(loop)*xposdist,i10,xlabel3,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial');
        t=text(leftpos+i8(loop)*xposdist,i10,xlabel2,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial');
        t=text(leftpos+i9(loop)*xposdist,i10,xlabel1,'units','normalized');set(t,'fontsize',9,'fontweight','bold','fontname','arial');
    
        t=text(i11,i12,{'ERA5 Inferred Annual ','   Exceedance Prob.'},'units','normalized');
            set(t,'fontsize',10,'fontweight','bold','fontname','arial');
    end

    set(gcf,'color','w');curpart=2;thisheight=9;highqualityfiguresetup;
end


if createfig1==1
    %Panel a: ERA5 max observed Tw std anoms
    stdanomarr=tw_regdecstdanoms_adj_era5;
    stdanomarr_jra55=tw_regdecstdanoms_adj_jra55;

    stdanomstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        stdanomofrec_era5(reg)=max(max(stdanomarr(reg,:,:)));
        stdanomofrec_jra55(reg)=max(max(stdanomarr_jra55(reg,:,:)));

        stdanomstoplot(regs_era5sz==reg)=stdanomofrec_era5(reg);

        %Implement region-exclusion policy: if mag of record differs by >=10% between ERA5 and JRA55
        if 100*((stdanomofrec_era5(reg)-stdanomofrec_jra55(reg))./stdanomofrec_era5(reg))>=10
            stdanomstoplot(regs_era5sz==reg)=NaN;
        end
    end

    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);stdanomstoplot};cmap=colormaps('whitelightreddarkred','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';1.5;'underlaycaxismax';5.0;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.54 0.8 0.42]);
    c=colorbar;c.Location='east outside';c.FontWeight='bold';c.Label.String='\sigma';
    c.Ticks=[2:5];c.TickLabels=[2:5];c.FontSize=cbfontsz+1;c.Label.FontSize=17;
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');
    title('Largest Extreme','fontsize',14,'fontweight','bold','fontname','arial');


    %Panel b: ratio of ERA5 return periods from GEV fits with and without observed maxes
    myratio=currecretpers_tw./currecretpers_withtopevents_tw;
    ratiostoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        ratiostoplot(regs_era5sz==reg)=log10(myratio(reg));
    end
    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);ratiostoplot};cmap=colormaps('lightyellowgreenblue','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0;'underlaycaxismax';2.7;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'colorbarticks';[0;0.48;1;1.48;2;2.48];'colorbarticklabels';{'1';'3';'10';'30';'100';'300'};...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.04 0.8 0.42]);
    c=colorbar;c.Location='east outside';c.FontWeight='bold';c.Label.String='Years';
    c.Ticks=[0;0.48;1;1.48;2;2.48];c.TickLabels={'1';'3';'10';'30';'100';'300'};c.FontSize=cbfontsz+1;
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');
    title('Return-Period Signature of Largest Extreme','fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig1_latest';curpart=2;highqualityfiguresetup;
end

if fig1_siversions==1
    %First for ERA5 T, then for JRA55 Tw

    %Panel a: ERA5 max observed T std anoms
    stdanomarr=t_regdecstdanoms_adj_era5;
    stdanomarr_jra55=t_regdecstdanoms_adj_jra55;

    stdanomstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        stdanomofrec_era5(reg)=max(max(stdanomarr(reg,:,:)));
        stdanomofrec_jra55(reg)=max(max(stdanomarr_jra55(reg,:,:)));

        stdanomstoplot(regs_era5sz==reg)=stdanomofrec_era5(reg);

        %Implement region-exclusion policy: if mag of record differs by >=10% between ERA5 and JRA55
        if 100*((stdanomofrec_era5(reg)-stdanomofrec_jra55(reg))./stdanomofrec_era5(reg))>=10
            stdanomstoplot(regs_era5sz==reg)=NaN;
        end
    end

    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);stdanomstoplot};cmap=colormaps('whitelightreddarkred','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';1.5;'underlaycaxismax';5.0;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.55 0.8 0.42]);
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');


    %Panel b: ratio of ERA5 return periods from GEV fits with and without observed maxes
    myratio=currecretpers_t./currecretpers_withtopevents_t;
    ratiostoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        ratiostoplot(regs_era5sz==reg)=log10(myratio(reg));
    end
    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);ratiostoplot};cmap=colormaps('lightyellowgreenblue','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0;'underlaycaxismax';2.7;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'colorbarticks';[0;0.48;1;1.48;2;2.48];'colorbarticklabels';{'1';'3';'10';'30';'100';'300'};...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.05 0.8 0.42]);
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig1_butwitht';curpart=2;highqualityfiguresetup;




    %Now again for Tw in JRA55
    %Panel a: JRA55 max observed Tw std anoms
    stdanomarr=tw_regdecstdanoms_adj_jra55;
    stdanomarr_era5=tw_regdecstdanoms_adj_era5;

    stdanomstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        stdanomofrec_jra55(reg)=max(max(stdanomarr(reg,:,:)));
        stdanomofrec_era5(reg)=max(max(stdanomarr_era5(reg,:,:)));

        stdanomstoplot(regs_era5sz==reg)=stdanomofrec_jra55(reg);

        %Implement region-exclusion policy: if mag of record differs by >=10% between ERA5 and JRA55
        if 100*((stdanomofrec_jra55(reg)-stdanomofrec_era5(reg))./stdanomofrec_jra55(reg))>=10
            stdanomstoplot(regs_era5sz==reg)=NaN;
        end
    end

    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);stdanomstoplot};cmap=colormaps('whitelightreddarkred','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';1.5;'underlaycaxismax';5.0;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.55 0.8 0.42]);
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');


    %Panel b: ratio of JRA55 return periods from GEV fits with and without observed maxes
    myratio=currecretpers_tw_jra55./currecretpers_withtopevents_tw_jra55;
    ratiostoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        ratiostoplot(regs_era5sz==reg)=log10(myratio(reg));
    end
    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);ratiostoplot};cmap=colormaps('lightyellowgreenblue','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0;'underlaycaxismax';2.7;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;'colorbarticks';[0;0.48;1;1.48;2;2.48];'colorbarticklabels';{'1';'3';'10';'30';'100';'300'};...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.05 0.8 0.42]);
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig1_butwithjra55';curpart=2;highqualityfiguresetup;
end


if createfig3==1
    variab='tw';
    yeararr=eval(['yearofmax_' variab '_adj_era5;']);
    yearstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr(reg);
        yearstoplot(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank(y-firstyear+1)=sum(yeararr==y);
    end

    %Panel a
    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.54 0.8 0.42]);
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');
    title('Year of Largest Extreme','fontsize',14,'fontweight','bold','fontname','arial');

    %Panel b
    yeararr_unadj=eval(['yearofmax_' variab '_era5;']);
    yearstoplot_unadj=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr_unadj(reg);
        yearstoplot_unadj(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank_unadj=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank_unadj(y-firstyear+1)=sum(yeararr_unadj==y);
    end

    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot_unadj};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.04 0.8 0.42]);
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');
    title('Year of Largest Extreme, Not GMST-Adjusted','fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig3_latest';curpart=2;highqualityfiguresetup;
end


if fig3_siversions==1
    %First for ERA5 T, then for JRA55 Tw

    yeararr=yearofmax_t_adj_era5;
    yearstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr(reg);
        yearstoplot(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank(y-firstyear+1)=sum(yeararr==y);
    end

    %Panel a
    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.55 0.8 0.42]);
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    %Panel b
    yeararr_unadj=yearofmax_t_era5;
    yearstoplot_unadj=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr_unadj(reg);
        yearstoplot_unadj(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank_unadj=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank_unadj(y-firstyear+1)=sum(yeararr_unadj==y);
    end

    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot_unadj};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.07 0.8 0.42]);
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig3_butfort';curpart=2;highqualityfiguresetup;


    %JRA55 Tw
    yeararr=yearofmax_tw_adj_jra55;
    yearstoplot=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr(reg);
        yearstoplot(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank(y-firstyear+1)=sum(yeararr==y);
    end

    %Panel a
    figure(888);clf;hold on;curpart=1;highqualityfiguresetup;
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.55 0.8 0.42]);
    t=text(-0.02,0.49,splabels{1},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    %Panel b
    exist yearofmax_tw_jra55;if ans==0;tmp=load(strcat(saveloc,'jra55output_regions.mat'));yearofmax_tw_jra55=tmp.yearofmax_tw_jra55;end
    yeararr_unadj=yearofmax_tw_jra55;
    yearstoplot_unadj=regs_era5sz; %to initialize
    for reg=1:numreg
        yearofrec_era5(reg)=yeararr_unadj(reg);
        yearstoplot_unadj(regs_era5sz==reg)=yearofrec_era5(reg);
    end
    yearrank_unadj=zeros(numyr,1);
    for y=firstyear:lastyear
        yearrank_unadj(y-firstyear+1)=sum(yeararr_unadj==y);
    end

    subplot(100,100,10000);
    data={double(lat2d_era5);double(lon2d_era5);yearstoplot_unadj};cmap=colormaps('classy rainbow','more','pale');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';firstyear;'underlaycaxismax';lastyear;'underlaystepsize';1;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gcf,'color','w');set(gca,'position',[0.1 0.07 0.8 0.42]);
    t=text(-0.02,0.49,splabels{2},'units','normalized');set(t,'fontsize',14,'fontweight','bold','fontname','arial');

    figname='fig3_butforjra55';curpart=2;highqualityfiguresetup;
end


if fig4setup_clustering==1
    %Calculate percentage of all p95-Tw days that occur in the top 25% of years from ERA5 and models
    arr=tw_regdecstdanoms_adj_era5;
    topyearspercentage_era5=zeros(numreg,1);
    cbyyr=zeros(numreg,numyr);
    for reg=1:numreg
        p95thresh=quantile(reshape(tw_regdecstdanoms_adj_era5(reg,:,:),[365*numyr 1]),0.95);
        for y=1:numyr;cbyyr(reg,y)=sum(tw_regdecstdanoms_adj_era5(reg,:,y)>=p95thresh);end
        %Percent in top 25% of years for this reg
        c_sorted=sort(cbyyr(reg,:),'descend');
        topyearspercentage_era5(reg)=100*sum(c_sorted(1:round(numyr/4)))/sum(c_sorted);
    end

    for ensmem=1:canesm5_enssz
        for reg=1:numreg
            p95thresh=quantile(reshape(tw_stdanoms_adj_canesm5(ensmem,:,reg,:),[365*numyr 1]),0.95);
            clear cbyyr;for y=1:numyr;cbyyr(y)=sum(tw_stdanoms_adj_canesm5(ensmem,:,reg,y)>=p95thresh);end
            %Percent in top 25% of years for this reg
            c_sorted=sort(cbyyr,'descend');
            topyearspercentage_canesm5(ensmem,reg)=100*sum(c_sorted(1:round(numyr/4)))/sum(c_sorted);
        end
    end

    for ensmem=1:miroc6_enssz
        for reg=1:numreg
            p95thresh=quantile(reshape(tw_stdanoms_adj_miroc6(ensmem,:,reg,:),[365*numyr 1]),0.95);
            clear cbyyr;for y=1:numyr;cbyyr(y)=sum(tw_stdanoms_adj_miroc6(ensmem,:,reg,y)>=p95thresh);end
            %Percent in top 25% of years for this reg
            c_sorted=sort(cbyyr,'descend');
            topyearspercentage_miroc6(ensmem,reg)=100*sum(c_sorted(1:round(numyr/4)))/sum(c_sorted);
        end
    end

    for ensmem=1:30
        for reg=1:numreg
            p95thresh=quantile(reshape(tw_stdanoms_adj_mpige(ensmem,:,reg,:),[365*numyr 1]),0.95);
            clear cbyyr;for y=1:numyr;cbyyr(y)=sum(tw_stdanoms_adj_mpige(ensmem,:,reg,y)>=p95thresh);end
            %Percent in top 25% of years for this reg
            c_sorted=sort(cbyyr,'descend');
            topyearspercentage_mpige(ensmem,reg)=100*sum(c_sorted(1:round(numyr/4)))/sum(c_sorted);
        end
    end

    %Make 2D array versions
    topyearspercentagetoplot_era5=0.*regs_era5sz;
    for reg=1:numreg;topyearspercentagetoplot_era5(regs_era5sz==reg)=topyearspercentage_era5(reg);end

    topyearspercentagetoplot_canesm5=NaN.*ones(canesm5_enssz,size(regs_canesm5sz,1),size(regs_canesm5sz,2));
    for ensmem=1:canesm5_enssz
        myarr=squeeze(topyearspercentage_canesm5(ensmem,:));
        for reg=1:numreg;topyearspercentagetoplot_canesm5(ensmem,regs_canesm5sz==reg)=myarr(reg);end
    end
    topyearspercentagetoplot_miroc6=NaN.*ones(miroc6_enssz,size(regs_miroc6sz,1),size(regs_miroc6sz,2));
    for ensmem=1:miroc6_enssz
        myarr=squeeze(topyearspercentage_miroc6(ensmem,:));
        for reg=1:numreg;topyearspercentagetoplot_miroc6(ensmem,regs_miroc6sz==reg)=myarr(reg);end
    end
    topyearspercentagetoplot_mpige=NaN.*ones(30,size(regs_mpigesz,1),size(regs_mpigesz,2));
    for ensmem=1:30
        myarr=squeeze(topyearspercentage_mpige(ensmem,:));
        for reg=1:numreg;topyearspercentagetoplot_mpige(ensmem,regs_mpigesz==reg)=myarr(reg);end
    end

    %Interpolate ERA5 and other models to MPI res for purposes of comparison
    topyearspctage_era5_mpigeres=interp2(lon2d_era5,lat2d_era5,topyearspercentagetoplot_era5,lon2d_mpige,lat2d_mpige);
    topyearspctage_canesm5_mpigeres=NaN.*ones(canesm5_enssz,96,192);
    for ensmem=1:canesm5_enssz
        topyearspctage_canesm5_mpigeres(ensmem,:,:)=interp2(lon2d_canesm5,lat2d_canesm5,squeeze(topyearspercentagetoplot_canesm5(ensmem,:,:)),lon2d_mpige,lat2d_mpige);
    end
    topyearspctage_miroc6_mpigeres=NaN.*ones(miroc6_enssz,96,192);
    for ensmem=1:miroc6_enssz
        topyearspctage_miroc6_mpigeres(ensmem,:,:)=interp2(lon2d_miroc6,lat2d_miroc6,squeeze(topyearspercentagetoplot_miroc6(ensmem,:,:)),lon2d_mpige,lat2d_mpige);
    end
    topyearspctage_mpige=topyearspercentagetoplot_mpige;
    topyearspercentage_xensmean=(squeeze(mean(topyearspctage_canesm5_mpigeres))+squeeze(mean(topyearspctage_miroc6_mpigeres))+squeeze(mean(topyearspctage_mpige)))./3;


    %Add model/ERA5 comparison
    clustering_modelera5ratio=NaN.*ones(size(lat2d_mpige));
    for i=1:size(lat2d_mpige,1)
        for j=1:size(lat2d_mpige,2)
            %if ~isnan(signifdiff(i,j)) %ensemble significance test result
                clustering_modelera5ratio(i,j)=topyearspercentage_xensmean(i,j)./topyearspctage_era5_mpigeres(i,j);
            %end
        end
    end
end


if fig4setup_duration==1
    %Calculate longest-duration >=2-std-anom events in ERA5
    longestdurtoplot_era5=0.*regs_era5sz; %just to initialize
    for reg=1:numreg
        longestdursofar=0;thisdur=0;
        for y=1:numyr
            for doy=1:365
                if tw_regdecstdanoms_adj_era5(reg,doy,y)>=2 %extr heat starts or continues
                    thisdur=thisdur+1;
                else
                    if thisdur>longestdursofar
                        longestdursofar=thisdur;
                        longestdur_date(reg,1)=y+firstyear-1;longestdur_date(reg,2)=doy-1; %end date of sequence
                    end
                    thisdur=0;
                end
            end
        end
        longestdurbyreg(reg)=longestdursofar;

        longestdurtoplot_era5(regs_era5sz==reg)=longestdurbyreg(reg);
    end
    %Interpolate ERA5 and other models (2-std-anom versions) to MPI res for purposes of comparison
    longestdurtoplot_era5_mpigeres=interp2(lon2d_era5,lat2d_era5,longestdurtoplot_era5,lon2d_mpige,lat2d_mpige);
    longestdurtoplot_canesm5_mpigeres=interp2(lon2d_canesm5,lat2d_canesm5,longestdurtoplot_canesm52,lon2d_mpige,lat2d_mpige);
    longestdurtoplot_miroc6_mpigeres=interp2(lon2d_miroc6,lat2d_miroc6,longestdurtoplot_miroc62,lon2d_mpige,lat2d_mpige);

    for i=1:canesm5_enssz;ldp_canesm5_byensmem_mpigeres(i,:,:)=interp2(lon2d_canesm5,lat2d_canesm5,...
            squeeze(longestdurtoplot_canesm52_byensmem(i,:,:)),lon2d_mpige,lat2d_mpige);end
    for i=1:miroc6_enssz;ldp_miroc6_byensmem_mpigeres(i,:,:)=interp2(lon2d_miroc6,lat2d_miroc6,...
            squeeze(longestdurtoplot_miroc62_byensmem(i,:,:)),lon2d_mpige,lat2d_mpige);end

    %Analyze longest-duration events averaged from each ensemble
    longestdurtoplot_xensmean=(longestdurtoplot_canesm5_mpigeres+longestdurtoplot_miroc6_mpigeres+longestdurtoplot_mpige2)./3;

    canesm5_dur_p10=squeeze(quantile(ldp_canesm5_byensmem_mpigeres,0.10));canesm5_dur_p90=squeeze(quantile(ldp_canesm5_byensmem_mpigeres,0.90));
    miroc6_dur_p10=squeeze(quantile(ldp_miroc6_byensmem_mpigeres,0.10));miroc6_dur_p90=squeeze(quantile(ldp_miroc6_byensmem_mpigeres,0.90));
    mpige_dur_p10=squeeze(quantile(longestdurtoplot_mpige2_byensmem,0.10));mpige_dur_p90=squeeze(quantile(longestdurtoplot_mpige2_byensmem,0.90));
    canesm5_dur_p25=squeeze(quantile(ldp_canesm5_byensmem_mpigeres,0.25));canesm5_dur_p75=squeeze(quantile(ldp_canesm5_byensmem_mpigeres,0.75));
    miroc6_dur_p25=squeeze(quantile(ldp_miroc6_byensmem_mpigeres,0.25));miroc6_dur_p75=squeeze(quantile(ldp_miroc6_byensmem_mpigeres,0.75));
    mpige_dur_p25=squeeze(quantile(longestdurtoplot_mpige2_byensmem,0.25));mpige_dur_p75=squeeze(quantile(longestdurtoplot_mpige2_byensmem,0.75));

    %Only will plot if 90% of model ens mems are above ERA5, or if 90% are
    %below, for two of the three ensembles -- DECIDED TO OMIT THIS AS UNNECESSARILY STRICT
    %A minimum ERA5 duration of 5 days is needed to make this comparison meaningful
    signifdiff=NaN.*ones(size(lat2d_mpige,1),size(lat2d_mpige,2));
    for i=1:size(lat2d_mpige,1)
        for j=1:size(lat2d_mpige,2)
            if longestdurtoplot_era5_mpigeres(i,j)>=5
                orig=0;
                if orig==1 %10/90
                    if (canesm5_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && miroc6_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (canesm5_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (miroc6_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p10(i,j)>=longestdurtoplot_era5_mpigeres(i,j))
                        signifdiff(i,j)=1;
                    elseif (canesm5_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && miroc6_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (canesm5_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (miroc6_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p90(i,j)<=longestdurtoplot_era5_mpigeres(i,j))
                        signifdiff(i,j)=2;
                    end
                else %25/75
                    if (canesm5_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && miroc6_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (canesm5_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (miroc6_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p25(i,j)>=longestdurtoplot_era5_mpigeres(i,j))
                        signifdiff(i,j)=1;
                    elseif (canesm5_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && miroc6_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (canesm5_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j)) ||...
                            (miroc6_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j) && mpige_dur_p75(i,j)<=longestdurtoplot_era5_mpigeres(i,j))
                        signifdiff(i,j)=2;
                    end
                end
            end
        end
    end

    %Add model/ERA5 comparison
    duration_modelera5pctdiff=NaN.*ones(size(lat2d_mpige));duration_modelera5ratio=NaN.*ones(size(lat2d_mpige));
    for i=1:size(lat2d_mpige,1)
        for j=1:size(lat2d_mpige,2)
            if longestdurtoplot_era5_mpigeres(i,j)>=5 %otherwise comparison is not meaningful
                %if ~isnan(signifdiff(i,j)) %ensemble significance test result
                    duration_modelera5pctdiff(i,j)=100.*(longestdurtoplot_xensmean(i,j)-longestdurtoplot_era5_mpigeres(i,j))./longestdurtoplot_era5_mpigeres(i,j);
                    duration_modelera5ratio(i,j)=longestdurtoplot_xensmean(i,j)./longestdurtoplot_era5_mpigeres(i,j);
                %end
            end
        end
    end
end

if fig4setup_intensity==1
    %ERA5
    stdanomarr=tw_regdecstdanoms_adj_era5;
    stdanomarr_jra55=tw_regdecstdanoms_adj_jra55;

    mostintense_era5=regs_era5sz; %to initialize
    for reg=1:numreg
        stdanomofrec_era5(reg)=max(max(stdanomarr(reg,:,:)));
        stdanomofrec_jra55(reg)=max(max(stdanomarr_jra55(reg,:,:)));

        mostintense_era5(regs_era5sz==reg)=stdanomofrec_era5(reg);

        %Implement region-exclusion policy: if mag of record differs by >=10% between ERA5 and JRA55
        if 100*((stdanomofrec_era5(reg)-stdanomofrec_jra55(reg))./stdanomofrec_era5(reg))>=10
            mostintense_era5(regs_era5sz==reg)=NaN;
        end
    end

    %Models
    stdanomofrec=NaN.*ones(canesm5_enssz,numreg);
    for ensmem=1:canesm5_enssz
        for reg=1:numreg;stdanomofrec(ensmem,reg)=max(max(tw_stdanoms_adj_canesm5(ensmem,:,reg,:)));end
    end
    canesm5mean=mean(stdanomofrec);

    stdanomofrec=NaN.*ones(miroc6_enssz,numreg);
    for ensmem=1:miroc6_enssz
        for reg=1:numreg;stdanomofrec(ensmem,reg)=max(max(tw_stdanoms_adj_miroc6(ensmem,:,reg,:)));end
    end
    miroc6mean=mean(stdanomofrec);

    stdanomofrec=NaN.*ones(30,numreg);
    for ensmem=1:30
        for reg=1:numreg;stdanomofrec(ensmem,reg)=max(max(tw_stdanoms_adj_mpige(ensmem,:,reg,:)));end
    end
    mpigemean=mean(stdanomofrec);

    xensmean=(canesm5mean+miroc6mean+mpigemean)./3;


    mostintense_xensmean=regs_mpigesz;
    for reg=1:numreg
        mostintense_xensmean(regs_mpigesz==reg)=xensmean(reg);
    end

    %Model/ERA5 comparison
    %Interpolate ERA5 to MPI res for purposes of comparison
    mostintense_era5_mpigeres=interp2(lon2d_era5,lat2d_era5,mostintense_era5,lon2d_mpige,lat2d_mpige);
    mostintense_modelera5ratio=mostintense_xensmean./mostintense_era5_mpigeres;
end


if createfig4==1
    %Left column: clustering
    %Middle column: duration
    %Right column: intensity
    figure(500);clf;hold on;curpart=1;highqualityfiguresetup;
    ori='vert'; %'horiz' or 'vert'
    if strcmp(ori,'horiz')
        lefts=[0.02;0.02;0.35;0.35;0.68;0.68];bottoms=[0.6;0.29;0.6;0.29;0.6;0.29];wwidth=0.31;hheight=0.38; %2 rows, 3 cols
        cbloc='southoutside';cbleftadj=0.03;cbbotadj=0.05;cbw=wwidth-0.06;cbh=0.015;
        cbw_ratios=cbw/2;cbh_ratios=cbh;cbw_addladj=cbw_ratios;cbh_addladj=0;
    else
        lefts=[0.005;0.51;0.005;0.51;0.005;0.51];bottoms=[0.67;0.67;0.34;0.34;0.01;0.01];wwidth=0.435;hheight=0.29; %3 rows, 2 cols
        cbloc='eastoutside';cbleftadj=wwidth+0.01;cbbotadj=0.03;cbw=0.015;cbh=hheight-0.06;
        cbw_ratios=cbw;cbh_ratios=cbh/2;cbw_addladj=0;cbh_addladj=cbh_ratios;
    end
    cbfontsz=10;

    %1. Clustering
    %a. ERA5
    k=1;
    data={double(lat2d_era5);double(lon2d_era5);topyearspercentagetoplot_era5};cmap=colormaps('whitelightreddarkred','more','not');
    cmax=25;cstep=2; %colorbar choices based on ERA5
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';35;'underlaycaxismax';75;'underlaystepsize';5;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'omitfirstsubplotcolorbar';1;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    c1=colorbar;c1.Location=cbloc;c1.FontWeight='bold';c1.Ticks=[35:10:75];c1.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw cbh];
    c1.Label.String='Percent';c1.FontSize=cbfontsz;
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');
    title('ERA5','fontsize',12,'fontweight','bold','fontname','arial');
    t=text(2.6,1.5,'Clustering');set(t,'fontsize',15,'fontweight','bold','fontname','arial');

    %b. Model/ERA5 comparison
    subplot(100,100,10000);k=2;
    %Set up colormap
    thisarr=clustering_modelera5ratio;
    maxfactor=1.5;cmap=flipud(colormaps('bluewhiteorange','more','not'));cmap=cmap(20:130,:);
    [cmapbelow1,cmapabove1,cbticks_below1,cbticks_above1,cbticklabels_below1,cbticklabels_above1]=fractionalcolormap(maxfactor,cmap);
    caxismina=1/maxfactor;caxismaxa=1;caxisstep=0.1;
    arr_below1=thisarr.*(thisarr<caxismaxa);
    invalid=arr_below1==0;arr_below1(invalid)=NaN;
    hAxesa=axes;cmapa=colormap(hAxesa,gray);
    %Plot 1
    data={double(lat2d_mpige);double(lon2d_mpige);arr_below1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxismina;'underlaycaxismax';caxismaxa;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapbelow1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapbelow1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6;'omitfirstsubplotcolorbar';1};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c2a=colorbar;c2a.Location=cbloc;c2a.FontWeight='bold';c2a.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw_ratios cbh_ratios];
    c2a.Ticks=cbticks_below1;c2a.TickLabels=cbticklabels_below1;c2a.Limits=[1/maxfactor 1];c2a.FontSize=cbfontsz;
    %Set up second colormap and then plot again
    caxisminb=1;caxismaxb=maxfactor;caxisstep=0.2;
    arr_above1=thisarr.*(thisarr>caxismaxa);
    invalid=arr_above1==0;arr_above1(invalid)=NaN;
    hAxesb=axes;cmapb=colormap(hAxesb,gray);
    %Plot 2
    data={double(lat2d_mpige);double(lon2d_mpige);arr_above1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxisminb;'underlaycaxismax';caxismaxb;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapabove1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapabove1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6;'omitfirstsubplotcolorbar';1};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c2b=colorbar;c2b.Location=cbloc;c2b.FontWeight='bold';c2b.Position=[lefts(k)+cbleftadj+cbw_addladj bottoms(k)+cbbotadj+cbh_addladj cbw_ratios cbh_ratios];
    c2b.Ticks=cbticks_above1;c2b.TickLabels=cbticklabels_above1;c2b.Limits=[1 maxfactor];c2b.FontSize=cbfontsz;
    %c2a.Label.String='                     Ratio';
    curpos=c2a.Label.Position;c2a.Label.Position=[curpos(1)-0.24 curpos(2)];
    title('LE/ERA5 Ratio','fontsize',12,'fontweight','bold','fontname','arial');
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');



    %2. Duration 
    %a. ERA5
    subplot(100,100,10000);k=3;
    data={double(lat2d_era5);double(lon2d_era5);longestdurtoplot_era5};cmap=colormaps('whitelightreddarkred','more','not');
    cmax=25;cstep=2; %colorbar choices based on ERA5
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';0;'underlaycaxismax';cmax;'underlaystepsize';cstep;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'omitfirstsubplotcolorbar';1;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c3=colorbar;c3.Location=cbloc;c3.FontWeight='bold';c3.Ticks=[0:5:25];c3.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw cbh];
    c3.Label.String='Days';c3.FontSize=cbfontsz;
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');
    title('ERA5','fontsize',12,'fontweight','bold','fontname','arial');
    t=text(2.6,1.5,'Duration');set(t,'fontsize',15,'fontweight','bold','fontname','arial');

    
    %b. Model/ERA5 comparison
    subplot(100,100,10000);k=4;
    %Set up colormap
    thisarr=duration_modelera5ratio;
    maxfactor=3;cmap=flipud(colormaps('bluewhiteorange','more','not'));cmap=cmap(20:130,:);
    [cmapbelow1,cmapabove1,cbticks_below1,cbticks_above1,cbticklabels_below1,cbticklabels_above1]=fractionalcolormap(maxfactor,cmap);
    caxismina=1/maxfactor;caxismaxa=1;caxisstep=0.1;
    arr_below1=thisarr.*(thisarr<caxismaxa);
    invalid=arr_below1==0;arr_below1(invalid)=NaN;
    hAxesa=axes;cmapa=colormap(hAxesa,gray);
    %Plot 1
    data={double(lat2d_mpige);double(lon2d_mpige);arr_below1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxismina;'underlaycaxismax';caxismaxa;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapbelow1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapbelow1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c4a=colorbar;c4a.Location=cbloc;c4a.FontWeight='bold';c4a.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw_ratios cbh_ratios];
    c4a.Ticks=cbticks_below1;c4a.TickLabels=cbticklabels_below1;c4a.Limits=[1/maxfactor 1];c4a.FontSize=cbfontsz;
    %Set up second colormap and then plot again
    caxisminb=1;caxismaxb=maxfactor;caxisstep=1/maxfactor;
    arr_above1=thisarr.*(thisarr>caxismaxa);
    invalid=arr_above1==0;arr_above1(invalid)=NaN;
    hAxesb=axes;cmapb=colormap(hAxesb,gray);
    %Plot 2
    data={double(lat2d_mpige);double(lon2d_mpige);arr_above1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxisminb;'underlaycaxismax';caxismaxb;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapabove1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapabove1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c4b=colorbar;c4b.Location=cbloc;c4b.FontWeight='bold';c4b.Position=[lefts(k)+cbleftadj+cbw_addladj bottoms(k)+cbbotadj+cbh_addladj cbw_ratios cbh_ratios];
    c4b.Ticks=cbticks_above1;c4b.TickLabels=cbticklabels_above1;c4b.Limits=[1 maxfactor];c4b.FontSize=cbfontsz;
    %c4a.Label.String='                     Ratio';
    curpos=c4a.Label.Position;c4a.Label.Position=[curpos(1)-0.24 curpos(2)];
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');
    title('LE/ERA5 Ratio','fontsize',12,'fontweight','bold','fontname','arial');


    %3. Intensity
    %a. ERA5
    subplot(100,100,10000);k=5;
    data={double(lat2d_era5);double(lon2d_era5);mostintense_era5};cmap=colormaps('whitelightreddarkred','more','not');
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';1.5;'underlaycaxismax';5.0;'underlaystepsize';0.25;'underlaycolormap';cmap;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmap;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c5=colorbar;c5.Location=cbloc;c5.FontWeight='bold';c5.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw cbh];c5.Ticks=[2:5];c5.FontSize=cbfontsz;
    c5.Label.String='\sigma';c5.Label.FontSize=17;
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');
    title('ERA5','fontsize',12,'fontweight','bold','fontname','arial');
    t=text(2.6,1.5,'Intensity');set(t,'fontsize',15,'fontweight','bold','fontname','arial');

    %b. Model/ERA5 comparison
    subplot(100,100,10000);k=6;
    %Set up colormap
    thisarr=mostintense_modelera5ratio;
    maxfactor=1.5;cmap=flipud(colormaps('bluewhiteorange','more','not'));cmap=cmap(20:130,:);
    [cmapbelow1,cmapabove1,cbticks_below1,cbticks_above1,cbticklabels_below1,cbticklabels_above1]=fractionalcolormap(maxfactor,cmap);
    caxismina=1/maxfactor;caxismaxa=1;caxisstep=0.1;
    arr_below1=thisarr.*(thisarr<caxismaxa);
    invalid=arr_below1==0;arr_below1(invalid)=NaN;
    hAxesa=axes;cmapa=colormap(hAxesa,gray);
    %Plot 1
    data={double(lat2d_mpige);double(lon2d_mpige);arr_below1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxismina;'underlaycaxismax';caxismaxa;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapbelow1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapbelow1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6;'omitfirstsubplotcolorbar';1};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c6a=colorbar;c6a.Location=cbloc;c6a.FontWeight='bold';c6a.Position=[lefts(k)+cbleftadj bottoms(k)+cbbotadj cbw_ratios cbh_ratios];
    c6a.Ticks=cbticks_below1;c6a.TickLabels=cbticklabels_below1;c6a.Limits=[1/maxfactor 1];c6a.FontSize=cbfontsz;
    %Set up second colormap and then plot again
    caxisminb=1;caxismaxb=maxfactor;caxisstep=0.2;
    arr_above1=thisarr.*(thisarr>caxismaxa);
    invalid=arr_above1==0;arr_above1(invalid)=NaN;
    hAxesb=axes;cmapb=colormap(hAxesb,gray);
    %Plot 2
    data={double(lat2d_mpige);double(lon2d_mpige);arr_above1};
    vararginnew={'mapproj';'robinson';'datatounderlay';data;'underlaycaxismin';caxisminb;'underlaycaxismax';caxismaxb;...
        'underlaystepsize';caxisstep;'underlaycolormap';cmapabove1;
        'contour_underlay';0;'contourunderlayfill';1;'contourunderlaycolors';cmapabove1;'centeredon';0;...
        'overlaynow';0;'nonewfig';1;'nanstransparent';1;'colorbarfontsize';12;...
        'stateboundaries';0;'countryboundarycolor';colors('gray');'countryborderlinewidth';0.6;'omitfirstsubplotcolorbar';1};
    datatype='custom';region={-180;90;180;-60};plotModelData(data,region,vararginnew,datatype);
    set(gca,'position',[lefts(k) bottoms(k) wwidth hheight]);
    c6b=colorbar;c6b.Location=cbloc;c6b.FontWeight='bold';c6b.Position=[lefts(k)+cbleftadj+cbw_addladj bottoms(k)+cbbotadj+cbh_addladj cbw_ratios cbh_ratios];
    c6b.Ticks=cbticks_above1;c6b.TickLabels=cbticklabels_above1;c6b.Limits=[1 maxfactor];c6b.FontSize=cbfontsz;
    %c6a.Label.String='                     Ratio';
    curpos=c6a.Label.Position;c6a.Label.Position=[curpos(1)-0.24 curpos(2)];
    t=text(0.01,0.47,splabels{k},'units','normalized');set(t,'fontsize',12,'fontweight','bold','fontname','arial');
    title('LE/ERA5 Ratio','fontsize',12,'fontweight','bold','fontname','arial');

    set(gcf,'color','w');figname='fig4_latest';curpart=2;thisheight=8;highqualityfiguresetup;
end


if findmostsingularevents==1
    for reg=1:numreg
        [largestanom(reg),yofmax(reg)]=max(max(squeeze(tw_regdecstdanoms_adj_era5(reg,:,:))));
        [~,dateofmax(reg)]=max(squeeze(tw_regdecstdanoms_adj_era5(reg,:,yofmax(reg))));

        withsurroundingsmasked=tw_regdecstdanoms_adj_era5;
        if dateofmax(reg)>=8 && dateofmax(reg)<=358
            withsurroundingsmasked(reg,dateofmax(reg)-7:dateofmax(reg)+7,yofmax(reg))=NaN;
        elseif dateofmax(reg)<=7
            withsurroundingsmasked(reg,365+dateofmax(reg)-7:365,yofmax(reg)-1)=NaN;
            withsurroundingsmasked(reg,1:dateofmax(reg)+7,yofmax(reg))=NaN;
        elseif dateofmax(reg)>=359
            withsurroundingsmasked(reg,dateofmax(reg)-7:365,yofmax(reg))=NaN;
            withsurroundingsmasked(reg,1:dateofmax(reg)+7-365,yofmax(reg)+1)=NaN;
        end
        [secondlargestanom(reg),yof2ndmax(reg)]=max(max(squeeze(withsurroundingsmasked(reg,:,:))));
        [~,dateof2ndmax(reg)]=max(squeeze(withsurroundingsmasked(reg,:,yof2ndmax(reg))));

        singularness(reg)=largestanom(reg)-secondlargestanom(reg);

        largestval(reg)=tw_regdecvals_adj_era5(reg,dateofmax(reg),yofmax(reg));
        secondlargestval(reg)=tw_regdecvals_adj_era5(reg,dateof2ndmax(reg),yof2ndmax(reg));
    end

    %Tabulate for regions with singularness >=0.5 std anoms
    c=0;clear myarr;
    for reg=1:numreg
        if singularness(reg)>=0.5
            c=c+1;
            myarr(c,1)=reg;
            myarr(c,2)=singularness(reg);
            myarr(c,3)=largestval(reg);
            myarr(c,4)=secondlargestval(reg);
            myarr(c,5)=yofmax(reg)+firstyear-1;
            myarr(c,6)=DOYtoMonth(dateofmax(reg),myarr(c,3));
            myarr(c,7)=DOYtoDOM(dateofmax(reg),myarr(c,3));
        end
    end
    myarr=sortrows(myarr,2,'descend');
end
