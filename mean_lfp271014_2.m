% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This program calculate the mean lfp of the experience for each 
% channel and save them in matriecs. (the calculte is after reduce
% the sampling rate (Hz)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load('Ligeti_271014all_recording_cttl_times.mat');
load('ligeti270314_order_of_playing.mat');

% computing the beginnig and the ending time of each version (+10 sec in
% the end and 5 sec at the beginning)
differences = diff(cttl_times); 
end_of_ver = find(differences > 30);
end_vers = [cttl_times(end_of_ver) cttl_times(end)] + 10;

start_vers = [cttl_times(1) cttl_times(end_of_ver + 1) - 5];
version_order = [find(runP == 1); find(runP == 2); find(runP == 3); find(runP == 4)];

clear end_of_ver differences
% mkdir Ligeti271014_2_mean

for chan=14:16
    load(['files_after_filter/Ligeti_271014_channel_' num2str(chan) 'file_85_after_filter.mat'], 't', ['lfp_channel_' num2str(chan) '_file_85']);
    if (chan == 1)
        timeline = t(1:10:end);    
    end
    all_files = [];
    eval(['all_files = lfp_channel_' num2str(chan) '_file_85(1:10:end);']);
    % input all files into one matrix and create timeline
    for ii=86:104
        ii;
        load(['files_after_filter/Ligeti_271014_channel_' num2str(chan) 'file_' num2str(ii) '_after_filter.mat'], ['lfp_channel_' num2str(chan) '_file_' num2str(ii)]);
        cur_file = eval(['lfp_channel_' num2str(chan) '_file_' num2str(ii) ';']);
        all_files = [all_files cur_file(1:10:end)];
        if (chan == 1)  % all the channel has the same timeline 
            load(['files_after_filter/Ligeti_271014_channel_' num2str(chan) 'file_' num2str(ii) '_after_filter.mat'], 't');
            timeline = [timeline t(1:10:end)];    
        end
    end
    
    if (chan == 14)
%         save('timeline.mat', 'timeline');
        max_len = 0;
        for ii=1:20
            cur_len = length(all_files(find(timeline >= start_vers(ii),1):find(timeline <= end_vers(ii),1 ,'last')));
            if (cur_len > max_len)
                max_len = cur_len;
            end
        end
        all_lfp_means = zeros(max_len, 4)';
    end

    % split to 20 versions
    all_vers = zeros(max_len, 20)';
    for ii=1:20
        ii;
        e = find(timeline <= end_vers(ii));
        one_ver = all_files(find(timeline >= start_vers(ii),1):e(end));
        all_vers(ii, 1:length(one_ver)) = one_ver;
    end

    % calculate the mean of each version
    for ver_ind=1:4
        eval(['ver_' num2str(ver_ind) ' = all_vers(version_order(' num2str(ver_ind) ', :), :);']);
        eval(['mean_' num2str(ver_ind) ' = [ver_' num2str(ver_ind) '; mean(ver_' num2str(ver_ind) ')];']);
        eval(['all_lfp_means(' num2str(ver_ind) ', :) = all_lfp_means(' num2str(ver_ind) ', :) + sum(ver_' num2str(ver_ind) ');']);
    end
    
    save(['Ligeti271014_2_mean/channel_' num2str(chan) '_mean_lfp.mat'], 'mean_1','mean_2','mean_3', 'mean_4', 'all_vers');
    for ii=35:54
        eval(['clear lfp_channel_' num2str(chan) '_file_' num2str(ii)]);
    end
    clear all_vers mean_1 mean_2 mean_3 mean_4 ver_1 ver_2 ver_3 ver_4;
end

all_lfp_means = zeros(605001, 4)';

for chan=1:16
    load(['Ligeti271014_2_mean/channel_' num2str(chan) '_mean_lfp.mat']);
    
    % calculate the mean of each version
    for ver_ind=1:4
        eval(['all_lfp_means(' num2str(ver_ind) ', :) = all_lfp_means(' num2str(ver_ind) ', :) + sum(mean_' num2str(ver_ind) '(1:5, :));']);
    end
    
    clear mean_1 mean_2 mean_3 mean_4;
end

for i=1:length(all_lfp_means(:, 1))
    all_lfp_means(i, :) = all_lfp_means(i, :) / 80;
end
save('Ligeti271014_2_mean/all_mean_lfp.mat', 'all_lfp_means');
