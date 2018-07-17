%SEC-C_parfor_test

%Nader Shakibay Senobari, summer 2018

%A toy example fo making SEC-C parrelelized using parfor. For heavy cases,
%we recommend running SEC-C on separate MATLABs at the same time instead of
%using parfor

%test_data_producer and SEC-C functions should be added to your MATLAB
%path

n_threads=8; %number of CPU workers (i.e. threads) that you want to run SEC-C on them

k=2^12;  %k is a tunning parameter, we recommend assigning a power of two for k 
%(e.g. 2^12 for one day of 20 Hz data or 2^13 for one day of 50-100 Hz data) 
%for efficient performance. However, in all cases we advise running some 
%test cases to find values of k to optimize the run time. For more 
%information we refer you to the citation above.



%//////////////////////////////////////producing fake data/////////////
%in real-worl case you don't need this section as you already have your
%data and templates ready for matched filtering


%lets make a fake data directory, in real-world case you already have this
mkdir ('./Parkfield_data_dir')
%lets make a sub-directory, in real-world you would have many of this,
%possibly named by daily dates
data_dir='./Parkfield_data_dir/day_1/';
mkdir(data_dir);

%lets put some fake data on it, using test_FMF (Beauce et al., 2018)
%In this case, 5 stations, one channel, one template, 5 sec template
%duration, and 20 Hz

[data,templates,moveouts,weights]=test_data_producer(20,1,5,1,5);

%save the data
save([data_dir,'data.mat'],'data');
clear data
%save('./tmeplates.mat','templates');
%save('./moveouts.mat','moveouts');
%weights('./weights','weights');

% In this toy example we assume all 365 days include the same data
data_d = cell(1, 365);
data_d(:) = {[data_dir,'data.mat']};

%////////////////////////////////////////////////////////////////////////

%now create an output directory
CCC_dir='./Parkfield_CCC_out_dir/';
mkdir(CCC_dir);


%Start a paraller pool using n_threads
parpool('local',n_threads);


tic

parfor i=1:365

    data_ind=loading_data(data_d{i});
    
    tic
    [ccc_sum]=SEC_C(data_ind,templates,k,moveouts,weights,1);
    t{i}=toc;
    
    save_CCC_sum(ccc_sum,CCC_dir,i);
end
total_time=toc;

display(['matched filtering took ', num2str(total_time),'seconds']);
display(['SEC-C runtime per day took ' num2str(sum([t{:}])/(n_threads*365)),' seconds']);

delete(gcp('nocreate'))

