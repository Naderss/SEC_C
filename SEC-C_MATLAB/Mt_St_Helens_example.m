%%

%A toy example of performing template matching that includes, retrieving, 
%prepossessing, performing template matching using SEC-C and
%postprocessing results for Mt St Helens seismicity. Here we only do
%one day of template matching (2018-01-03) and only for one template. 
%User can simply loop over the days and templates and also include more
%stations (here we used 4 stations, one component).

%*************************************************************************
%This code uses irisFetch to interact with DMCs, so you should have
%irisfetch.m in your MATLAB path.
%You can find irisFetch from: https://ds.iris.edu/ds/nodes/dmc/software/downloads/irisFetch.m/
%You also should have SEC_C_fully_normalized.m in your path as well.
%**************************************************************************

% Preprossessings inspired by: https://blogs.mathworks.com/loren/2015/03/03/direct-access-to-seismological-data-using-matlab/

%irisFetch requires a Java JAR file (IRIS-WS.jar) to communicate with the 
%FDSN web services. Therefore, before using irisFetch, retrieve the 
%latest IRIS-WS java JAR file and make it available to MATLAB.
%Obtain the most recent IRIS-WS-2.0.x.jar from:
%http://ds.iris.edu/ds/nodes/dmc/software/downloads/IRIS-WS/

%e.g.
%javaaddpath('IRIS-WS-2.0.17.jar');


%Nader Shakibay Senobari, summer 2018




%/////////////////////////////////////////////////////
%////////////////////////////////////////////////////
%Retrive Catalog events 

%lets set some vars for Mt St Helen area and seimicity
minlat = 46.0189; maxlat = 46.418;
minlon = -122.455; maxlon = -121.8549;
starttime = '2018-01-03 00:00:00';
endtime = '2018-01-04 00:00:00';
maxmag = 5.0; minmag = 0.0;
  
  %retrieving event informations
  
  cata_ev=irisFetch.Events('boxcoordinates',[minlat,maxlat,minlon,maxlon],...
      'maximumMagnitude',maxmag,'minimumMagnitude',minmag,...
      'startTime',starttime,'endTime',endtime,...
      'BASEURL','http://service.iris.edu/fdsnws/event/1/');
  
  
  %retrieving station iformations (only PB network
  PB_sta = irisFetch.Stations('station','PB','*','*','EHZ',...
      'boxcoordinates',[minlat,maxlat,minlon,maxlon]);
  
  
  plot([cata_ev.PreferredLongitude],[cata_ev.PreferredLatitude],'r.',...
      'MarkerSize',14);hold on;
  plot([PB_sta.Longitude],[PB_sta.Latitude],'b^','MarkerFaceColor','b');
  
  %Mt St Helens
  %plot(-122.19086,46.19978,'g*');
  
  
  
  xlim([minlon,maxlon]);
  ylim([minlat,maxlat]);
  
  

  
 %% 
  %retrieving data
    for i=1:length(PB_sta)
      data_tr(i) = irisFetch.Traces(PB_sta(i).NetworkCode,PB_sta(i).StationCode,'*','EHZ','2018-01-03 00:00:00',...
          '2018-01-04 00:00:00','verbose');
    end
    %%
  

 %now lets filter the data between 1 to 10 Hz
  bandfilt_freq1 = 1;
  bandfilt_freq2 = 10;
  bandfilt_order = 4;
  
  for i=1:length(data_tr)
      data = (data_tr(i).data - mean(data_tr(i).data)) ./ data_tr(i).sensitivity;
      data=detrend(data);
      wn1 = bandfilt_freq1/data_tr(i).sampleRate;
      wn2 = bandfilt_freq2/data_tr(i).sampleRate;
      [f1,f2] = butter(bandfilt_order,[wn1 wn2],'stop');
      data = filter(f1,f2,data);
      data_tr(i).data_fil=data;
  end
  %%
  %Now lets extract some templates based on a catalog event
  %you can loop over all events
  
  ev_st=datenum(cata_ev(28).PreferredTime); %event start time
  day_st=data_tr(1).startTime;   %data start time
  
  %lets find the index for event waveforms
  ind_sample=etime(datevec(ev_st),datevec(day_st))*data_tr(1).sampleRate;
  ind_sample=round(ind_sample);
  
  %now lets pick the P arrival and make templates;
  for i=1:length(data_tr)
        L=figure;
        plot(data_tr(i).data_fil(ind_sample-200:ind_sample+1000));
        set(L, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.1, 1, 1]);

        [x_ind,y_ind]=ginput;
        x_ind=round(x_ind);
        
        tem(i).st_ind=ind_sample-300+x_ind;
        tem(i).tr=data_tr(i).data_fil(tem(i).st_ind:tem(i).st_ind+6*data_tr(i).sampleRate);
        
        close;
  end
  
 %% 
  %now perform the matched filter using SEC-C
 
  %note that in order to get SEC_C to work you should prepare the input
  %data carefully. Here we have data and templates for one components. so
  %wee need to make a 3 dimensional matrix for SEC_C's input as below:

templates(:,1,:)=[tem(:).tr];
clear data
data(:,1,:)=[data_tr(:).data_fil];

% we have 4 stations and we want to assign them equial weights:
weights=[1/4,1/4,1/4,1/4]';

% and moveouts based on pick times, actually moveouts can be assigned as pick
% times
moveouts=[tem(:).st_ind]';

%now perform matched filtering using SEC_C_fully_normalized

CC_sum=SEC_C_fully_normalized(data,templates,2^13,moveouts,weights,3);

%now plot CC_sum and 9*MAD (i.e. detection threshold)

plot(CC_sum);hold on;
plot([1,length(CC_sum)],9*[mad(CC_sum,1),mad(CC_sum,1)],'r');
  
%detecting seismic events using 9*MAD threshold and assuming 2 seconds
%seperation between detections
%%
[pks,locs]=findpeaks(CC_sum,'SortStr','descend','MinPeakDistance',2*data_tr(1).sampleRate,...
    'MinPeakHeight',9*mad(CC_sum,1));
%%

%locs are the time of detections (in sample) and pks are the values of CC
%for these detections. Both are useful information for further processes 
%such as locating earthquakes or relocating them. 
  
  
  
      
  
  
  
