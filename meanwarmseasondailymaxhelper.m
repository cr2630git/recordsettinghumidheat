%Helper for calculating mean warm-season daily-max variables
%Written for "Record-setting humid heat" manuscript
%Colin Raymond, Apr 2024

%Note: regdecvals must be input with dims region | day of year | year

for reghere=1:numreg
    for loophere=1:length(varstocalc)
        thisvar=varstocalc{loophere};

        if strcmp(thisvar,'t')
            peakmonth=peakmonth_t(reghere);regdecvals=t_vals;
        elseif strcmp(thisvar,'tw')
            peakmonth=peakmonth_tw(reghere);regdecvals=tw_vals;
        elseif strcmp(thisvar,'q')
            peakmonth=peakmonth_q(reghere);regdecvals=q_vals;
        end

        if peakmonth>=2 && peakmonth<=11
            wsdoystart=monthstarts(peakmonth-1);wsdoystop=monthstops(peakmonth+1);
            wslen=monthlens(peakmonth-1)+monthlens(peakmonth)+monthlens(peakmonth+1);

            meanwsmax=mean(reshape(squeeze(regdecvals(reghere,wsdoystart:wsdoystop,firstclimoyear:lastclimoyear)),[wslen*numclimoyears 1]),'omitnan');
            stdwsmax=std(reshape(squeeze(regdecvals(reghere,wsdoystart:wsdoystop,firstclimoyear:lastclimoyear)),[wslen*numclimoyears 1]),'omitnan');
        elseif peakmonth==1
            wsdoystart=monthstarts(12);wsdoystop=monthstops(2);
            wslen=monthlens(12)+monthlens(1)+monthlens(2);

            data1=reshape(squeeze(regdecvals(reghere,wsdoystart:365,firstclimoyear:lastclimoyear)),[monthlens(12)*numclimoyears 1]);
            data2=reshape(squeeze(regdecvals(reghere,1:wsdoystop,firstclimoyear:lastclimoyear)),[(monthlens(1)+monthlens(2))*numclimoyears 1]);
            meanwsmax=mean([data1;data2],'omitnan');
            stdwsmax=std([data1;data2],'omitnan');
        elseif peakmonth==12
            wsdoystart=monthstarts(11);wsdoystop=monthstops(1);
            wslen=monthlens(11)+monthlens(12)+monthlens(1);

            data1=reshape(squeeze(regdecvals(reghere,wsdoystart:365,firstclimoyear:lastclimoyear)),[(monthlens(11)+monthlens(12))*numclimoyears 1]);
            data2=reshape(squeeze(regdecvals(reghere,1:wsdoystop,firstclimoyear:lastclimoyear)),[monthlens(1)*numclimoyears 1]);
            meanwsmax=mean([data1;data2],'omitnan');
            stdwsmax=std([data1;data2],'omitnan');
        end
        if strcmp(thisvar,'t')
            meanwsmax_t(reghere)=meanwsmax;stdwsmax_t(reghere)=stdwsmax;
        elseif strcmp(thisvar,'tw')
            meanwsmax_tw(reghere)=meanwsmax;stdwsmax_tw(reghere)=stdwsmax;
        elseif strcmp(thisvar,'q')
            meanwsmax_q(reghere)=meanwsmax;stdwsmax_q(reghere)=stdwsmax;
        end
    end
end

