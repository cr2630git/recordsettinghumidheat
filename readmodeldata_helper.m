%Helper for reading model data
%Written for "Record-setting humid heat" manuscript
%Colin Raymond, Apr 2024

if histpart==1
    for d=1:length(decades)
        thisdec=decades{d};
        arraystart=arrayindices_start(d);arraystop=arrayindices_stop(d);
        filestart=fileindices_start(d);filestop=fileindices_stop(d);
        filepath_tas=strcat(filepath_tas_part1,thisdec,'.nc');
        filepath_huss=strcat(filepath_huss_part1,thisdec,'.nc');
    
        tmp1=ncread(filepath_tas,'tasmax')-273.15; %C
        tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
        for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
        tasmax(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp2;clear tmp3;
    
        tmp1=ncread(filepath_huss,'huss'); %kg/kg
        tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
        for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
        hussmean(:,:,arraystart:arraystop)=tmp3(:,:,filestart:filestop);clear tmp2;clear tmp3;
    end
end

if futpart==1
    nytoread=ny_fut;

    tmp1=ncread(filepath_tas,'tasmax',[1 1 1],[Inf Inf 365*nytoread])-273.15; %C
    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
    tasmax(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp2;clear tmp3;


    tmp1=ncread(filepath_huss,'huss',[1 1 1],[Inf Inf 365*nytoread]); %kg/kg
    tmp2=flip(permute(tmp1,[2 1 3]),1);clear tmp1;tmp3=zeros(size(tmp2));
    for sz=1:size(tmp2,3);tmp3(:,:,sz)=recenter(tmp2(:,:,sz));end
    hussmean(:,:,ny_hist_todo*365+1:end)=tmp3(:,:,365*(futfirstyr-1)+1:365*futlastyr);clear tmp2;clear tmp3;
end