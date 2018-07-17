%SEC_C: Super Efficient Cross-correlation MATlAB code for seismic data 

%SEC_C is designed to handle multiple stations, multiple components and 
%multiple templates matched filtering of time series data specifically for
%seismic applications in an efficient time with an efficient memory usage. 
%We have adopted the input parameter style as Beauce et al., 2017 matched filter 
%code (https://github.com/beridel/fast_matched_filter). For testing SEC_C 
%user also can use the Beauce et al., 2017 test code. 

%Normalization part inspired by Mass algorithm 
%(http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html)

%Tested with MATLAB R2017a and later versions; the movsum function and vectorized
%column dot product do not work with R2014a and earlier versions.

%Written by: Nader Shakibay Senobari, Spring 2018 
%Citation:  Shakibay Senobari et al., 2018, submitted to SRL. 


function [ccc_sum]=SEC_C(data,templates,k,moveouts,weights,chan_scale2one)

%Input parameters 

%data: 3D matrix [number of data samples, number of components, number  
%of stations]

%templates: 4D matrix [number of samples, number of components, number  
%of stations, number of templates]

%k: A tunning parameter, we recommend assigning a power of two for k 
%(e.g. 2^12 for one day of 20 Hz data or 2^13 for one day of 50-100 Hz data) 
%for efficient performance. However, in all cases we advise running some 
%test cases to find values of k to optimize the run time. For more 
%information we refer you to the citation above.

%moveouts: 2D matrix [number of station, number of templates]

%weights: 2D matrix [number of stations, number of templates]

%chan_scale2one: =1 if you want to CCC_sum to be scaled to 1 (divide
%weights by the number of channels), =0 if you don't want to.

%note: moveouts and weights can be extended to include components (i.e.
%become a 3D matrix). In this case user can easily change line 94 to
%account for that.

%output: 
%CCC_sum: 2D matrix [length of data-length of template+1, number of
%templates] 

%User also can save individual CCC for each channel-template, see below

%///////////////////////////////SEC_C/////////////////////////////////////
%let's gather some information about data and templates such as the length of
%templates (m), number of components (n_c), number of stations (n_s),
%number of tempplates (n_t) and length of data. 
[m, n_c, n_s, n_t] = size(templates);  
l_data=length(data(:,1,1));

ccc_sum(l_data-m+1,n_t)=0; %Preassigning CCC_sum and make sure CCC_sum doesn't change in each loop

%We want to find a station with largest moveout (i.e. the last station that
%detect the signal) and later on we pad zeros to the other stations 
%at the begining of the data to align them.
moveouts=max(moveouts)-moveouts; 

% if you want to scale CCC_sum to one assign chan_scale2one=1
if chan_scale2one==1
    weights=weights./n_c;
end

    for j=1:n_s %loop over stations
        for jj=1:n_c %loop over componets
            
            y=squeeze(templates(:,jj,j,:)); % get the template data for each channel
            sumy2 = sqrt(sum(y.^2)); %do some preprocessing for normalization, we need this later on 
            
            %buffer is a matlab function that divides the data into 
            %pieces with the length of k and with overlaps of m-1 
            %samples and make a matrix s from the data vector. 
            s=buffer(data(:,jj,j),k,m-1); 
            
            sumx2_t=sqrt(movsum(s.^2,[m-1,0])); % another preprocessing step for normalization

            %Calculating cross-correlation (CC) in frequency domain 
            y = y(end:-1:1,:);  %reverse the templates
            y(m+1:k,:) = 0;     %padding with zeros
            X=fft(s);           %transfering to the frequency domain for the data
            Y=fft(y);           %transfering to the frequency domain for templates

            for i=1:n_t %number of templates
                
                Z = X.*Y(:,i);   %do the dot product
                z = ifft(Z);     %going back to the time domain
                
                ccha=z(m:k,:)./(sumx2_t(m:k,:)*sumy2(i)); %devide by the normalization factor
                
                %calculate the CCC and sum over stations and components
                ccc_sum(:,i)=weights(j,i).*([zeros(1,moveouts(j,i)),ccha(m:l_data-moveouts(j,i))])'+ccc_sum(:,i);
                
            end

        end

    end

end

