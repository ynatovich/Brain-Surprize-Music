%% file 1- creating the cttl timeline
%%

% this timeline is the same for all electrodes, therefor i made it from
% channel 6 and used the output for all other channels. 

cttl_times = [];
for ii = 85:104       % the file numbers; 105 doesn't have cttl data. 
    ii
    if ii<100
        load(['F271014-00100100' num2str(ii) '.mat'], 'CTTL_049_Up'); %the computer's UP ticks, in computer tick units.
        load(['F271014-00100100' num2str(ii) '.mat'], 'CTTL_049_TimeBegin'); % time of first tick in file, in seconds.
        load(['F271014-00100100' num2str(ii) '.mat'], 'CTTL_049_KHz'); % sampling freq. of the ticks.
        cttl_times = [cttl_times CTTL_049_TimeBegin+CTTL_049_Up/CTTL_049_KHz/1000]; % timeline is in seconds
    else
        load(['F271014-0010010' num2str(ii) '.mat'], 'CTTL_049_Up'); %the computer's UP ticks, in computer tick units.
        load(['F271014-0010010' num2str(ii) '.mat'], 'CTTL_049_TimeBegin'); % time of first tick in file, in seconds.
        load(['F271014-0010010' num2str(ii) '.mat'], 'CTTL_049_KHz'); % sampling freq. of the ticks.
        cttl_times = [cttl_times CTTL_049_TimeBegin+CTTL_049_Up/CTTL_049_KHz/1000]; % timeline is in seconds
    end
end

plot(cttl_times, '.') %there are indeed 20 rotations, this means that the following files 36-39
%are redundant, kept recording after the sound was over. this is why they don't have cttl data. 
save('Ligeti_271014all_recording_cttl_times.mat','cttl_times') %save the times







%% file 2- filtering the raw input to extract the LFP

%% run the filter on all files

b=fir1(270*3,150/27271*2);     

for ii = 1:16       % 16 channels overall
    ii
    for jj = 85:104     % the file numbers
        jj
        if jj < 100
            if ii<10
                load(['F271014-00100100' num2str(jj) '.mat']', ['CRAW_00' num2str(ii)])
                eval(['CRAW_00' num2str(ii) '_' num2str(jj) '=double(CRAW_00' num2str(ii) ');'])
                load(['F271014-00100100' num2str(jj) '.mat']', ['CRAW_00' num2str(ii) '_KHz']);
                load(['F271014-00100100' num2str(jj) '.mat'], ['CRAW_00' num2str(ii) '_TimeBegin']);
                l = length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                t = eval(['CRAW_00' num2str(ii) '_TimeBegin']) + ((1:length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)])))-1)/eval(['CRAW_00' num2str(ii) '_KHz'])/1000;
                temp = filtfilt(b,1,eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                eval(['lfp_channel_' num2str(ii) '_file_' num2str(jj) '= temp;'])
                save(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'],'t', ['CRAW_00' num2str(ii) '_' num2str(jj)],['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
                clear(['CRAW_00' num2str(ii)],['CRAW_00' num2str(ii) '_' num2str(jj)],['CRAW_00' num2str(ii) '_KHz'],...
                    ['CRAW_00' num2str(ii) '_TimeBegin'],'l','t','temp',['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
            else
                load(['F271014-00100100' num2str(jj) '.mat']', ['CRAW_0' num2str(ii)])
                eval(['CRAW_00' num2str(ii) '_' num2str(jj) '=double(CRAW_0' num2str(ii) ');'])
                load(['F271014-00100100' num2str(jj) '.mat']', ['CRAW_0' num2str(ii) '_KHz']);
                load(['F271014-00100100' num2str(jj) '.mat'], ['CRAW_0' num2str(ii) '_TimeBegin']);
                l = length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                t = eval(['CRAW_0' num2str(ii) '_TimeBegin']) + ((1:length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)])))-1)/eval(['CRAW_0' num2str(ii) '_KHz'])/1000;
                temp = filtfilt(b,1,eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                eval(['lfp_channel_' num2str(ii) '_file_' num2str(jj) '= temp;'])
                save(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'],'t', ['CRAW_00' num2str(ii) '_' num2str(jj)],['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
                clear(['CRAW_0' num2str(ii)],['CRAW_00' num2str(ii) '_' num2str(jj)],['CRAW_0' num2str(ii) '_KHz'],...
                    ['CRAW_0' num2str(ii) '_TimeBegin'],'l','t','temp',['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
            end
        else
            if ii<10
                load(['F271014-0010010' num2str(jj) '.mat']', ['CRAW_00' num2str(ii)])
                eval(['CRAW_00' num2str(ii) '_' num2str(jj) '=double(CRAW_00' num2str(ii) ');'])
                load(['F271014-0010010' num2str(jj) '.mat']', ['CRAW_00' num2str(ii) '_KHz']);
                load(['F271014-0010010' num2str(jj) '.mat'], ['CRAW_00' num2str(ii) '_TimeBegin']);
                l = length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                t = eval(['CRAW_00' num2str(ii) '_TimeBegin']) + ((1:length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)])))-1)/eval(['CRAW_00' num2str(ii) '_KHz'])/1000;
                temp = filtfilt(b,1,eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                eval(['lfp_channel_' num2str(ii) '_file_' num2str(jj) '= temp;'])
                save(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'],'t', ['CRAW_00' num2str(ii) '_' num2str(jj)],['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
                clear(['CRAW_00' num2str(ii)],['CRAW_00' num2str(ii) '_' num2str(jj)],['CRAW_00' num2str(ii) '_KHz'],...
                    ['CRAW_00' num2str(ii) '_TimeBegin'],'l','t','temp',['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
            else
                load(['F271014-0010010' num2str(jj) '.mat']', ['CRAW_0' num2str(ii)])
                eval(['CRAW_00' num2str(ii) '_' num2str(jj) '=double(CRAW_0' num2str(ii) ');'])
                load(['F271014-0010010' num2str(jj) '.mat']', ['CRAW_0' num2str(ii) '_KHz']);
                load(['F271014-0010010' num2str(jj) '.mat'], ['CRAW_0' num2str(ii) '_TimeBegin']);
                l = length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                t = eval(['CRAW_0' num2str(ii) '_TimeBegin']) + ((1:length(eval(['CRAW_00' num2str(ii) '_' num2str(jj)])))-1)/eval(['CRAW_0' num2str(ii) '_KHz'])/1000;
                temp = filtfilt(b,1,eval(['CRAW_00' num2str(ii) '_' num2str(jj)]));
                eval(['lfp_channel_' num2str(ii) '_file_' num2str(jj) '= temp;'])
                save(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'],'t', ['CRAW_00' num2str(ii) '_' num2str(jj)],['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
                clear(['CRAW_0' num2str(ii)],['CRAW_00' num2str(ii) '_' num2str(jj)],['CRAW_0' num2str(ii) '_KHz'],...
                    ['CRAW_0' num2str(ii) '_TimeBegin'],'l','t','temp',['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
            end
        end
    end
end


%% for each channel, concatenate results from all  14 files

% concatenate the timeline:
for ii = 1
    ii
    eval(['t_' num2str(ii) '=[];'])
    for jj = 85:104
        load(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'], 't');
        eval(['t_' num2str(ii) '=[t_' num2str(ii),',t];'])
        clear('t')
    end
    save(['Ligeti_271014_channel_' num2str(ii) '_timeline.mat'], ['t_' num2str(ii)]);  
    clear(['t_' num2str(ii)]);
end
% in fact, all the timelines are the same for all channels...
%%
%  concatenate the lfp
for ii = 1:16
    eval(['lfp_' num2str(ii) '=[];']);
    for jj = 85:104
        jj
        load(['Ligeti_271014_channel_' num2str(ii) 'file_' num2str(jj) '_after_filter.mat'], ['lfp_channel_' num2str(ii) '_file_' num2str(jj)]);
        eval(['lfp_' num2str(ii) '=[lfp_' num2str(ii),',[lfp_channel_' num2str(ii) '_file_' num2str(jj) ']];'])
        clear(['lfp_channel_' num2str(ii) '_file_' num2str(jj)])
    end
    save(['Ligeti_271014_channel_' num2str(ii) '_lfp.mat'], ['lfp_' num2str(ii)]);
    clear(['lfp_' num2str(ii)]);
end







%% file 3
%% create the mean_lfp of each electrode to each of the 4 versions. 

load('Ligeti_271014all_recording_cttl_times.mat')
figure
plot(cttl_times)
differences =diff(cttl_times); 
plot(differences); 
endings = []
for ii = 1:length(differences)
    if differences(ii) > 11
        endings(end+1) = ii;
    end
end

endings(20) = length(cttl_times)

ending_times = cttl_times(endings); 
ending_times_10 = ending_times + 10 %take 10 more seconds after cttl times, 
% since it's every 10 seconds and it could be that the recording kept playing there.

beg_time = zeros(1,20); ; 
beg_time(1) = cttl_times(1)

for ii = 2:20
    beg_time(ii) = cttl_times(endings(ii-1)+1); 
end

%take 5 secs before onset of sound, to see basal activity level:
beg_time_5 = zeros(1,20);
for ii = 1:20
    if ii == 1
        beg_time_5(ii) = beg_time(ii);  % since this is the beginning of the whole recording. 
    else
        beg_time_5(ii) = beg_time(ii) - 5;
    end
end

%divide the whole LFP into 20 playings, each one in the same length, which
%is:
lengths = ending_times_10 - beg_time_5; 

load('ligeti270314_order_of_playing') %runP contains the version order. 

%check length of trials
load('Ligeti_271014_channel_1_timeline.mat', 't_1'); 
load(['Ligeti_271014_channel_1_lfp.mat']);

t = t_1;
all_trial_lengths = zeros(1,20); 

for ii = 1:20
    beg_t = beg_time_5(ii); 
    end_t = ending_times_10(ii); 
    beg_dist = abs(t-beg_t); 
    [least_dist_beg, tbeg_ind] = min(beg_dist); % this is the place in the t
        % vector that is closest to the desired beginning time in beg_t
    end_dist = abs(t-end_t);
    [least_dist_beg, tend_ind] = min(end_dist); %place in the t vector that is 
        %closest to the desired end time in end_t
    lfp_trial = lfp_1(tbeg_ind : tend_ind); 
    trial_len = length(lfp_trial);
    all_trial_lengths(ii) = trial_len;
end
%%
%lenghts of the trials are 5940005 (first without 5 more seconds), 6050003, 
%6050004,6050005. 
%maximum of lengths is 6054805, so i create a matrix of zeros(6054805,5) 
% for each version, and put the trials in there.
max_val = 6050005;
min_val = 5940005;

for jj = 1:16 % go over 16 channels and divide the lfp
    jj
    load(['Ligeti_271014_channel_' num2str(jj) '_lfp.mat']);
    orig_1 = zeros(max_val, 5);
    orig_1_lengths = [];
    ident_2 = zeros(max_val, 5);
    ident_2_lengths = [];
    all_up_3 = zeros(max_val, 5);
    all_up_3_lengths = [];
    inverted_4 = zeros(max_val, 5);
    inverted_4_lengths = [];
    for ii = 1:20       % 20 playings of the Ligeti tune
        ii
        beg_t = beg_time_5(ii);
        end_t = ending_times_10(ii);
        beg_dist = abs(t-beg_t);
        [least_dist_beg, tbeg_ind] = min(beg_dist);
        end_dist = abs(t-end_t);
        [least_dist_beg, tend_ind] = min(end_dist);
        eval(['lfp_trial = lfp_' num2str(jj) '(tbeg_ind : tend_ind);']);
        trial_len = length(lfp_trial);
        if trial_len < max_val
            if trial_len == min_val
                lfp_trial = [zeros(max_val-min_val, 1); lfp_trial'];
            elseif trial_len == max_val - 1
                lfp_trial(end+1) = 0;
            elseif trial_len == max_val - 2
                lfp_trial(end+1) = 0;
                lfp_trial(end+1) = 0;
            end
        end
        play = runP(ii);
        if play == 1
            orig_1(:, length(orig_1_lengths)+1) = lfp_trial;
            orig_1_lengths = [orig_1_lengths trial_len];
        elseif play == 2
            ident_2(:,length(ident_2_lengths)+1) = lfp_trial;
            ident_2_lengths = [ident_2_lengths trial_len];
        elseif play == 3
            all_up_3(:, length(all_up_3_lengths)+1) = lfp_trial;
            all_up_3_lengths = [all_up_3_lengths trial_len];
        else
            inverted_4(:, length(inverted_4_lengths)+1) = lfp_trial;
            inverted_4_lengths = [inverted_4_lengths trial_len];
        end
    end
    clear(['lfp_' num2str(jj)]);
    %average over the 5 repeats
    mean_lfp = zeros(max_val,4);
    mean_lfp(:,1) = mean(orig_1,2);
    mean_lfp(:,2) = mean(ident_2,2);
    mean_lfp(:,3) = mean(all_up_3,2);
    mean_lfp(:,4) = mean(inverted_4,2);
    % save results
    save(['Ligeti_271014_channel_' num2str(jj) '4_lfp_matrixes_and_mean_lfp.mat'], 'orig_1',...
        'ident_2','all_up_3', 'inverted_4', 'mean_lfp');
    clear('orig_1','ident_2','all_up_3', 'inverted_4', 'mean_lfp');
end






%% file 4
%% find times and amplitude of the main 208 notes in the piece, for the first version:
% this is the same for all dates because the same pieces were played.
% therefore i can use the results i found for the original date. 
% they were saved here:

% save('ind&ampd_208.mat', 'ind208', 'amps208')

% save('three_versions_208_note_amps_and_indexes_27_4.mat', 'fsig1','fsig2', 'fsig3', 'fsig4',...
% 'ind2_208','ind3_208', 'ind4_208', 'amps2_208', 'amps3_208', 'amps4_208'); 

%% file 5
%% this file checks the times for the signal, so i don't need to do this again. 

%% file 6
%% this file gets the note order- i already have it, it's the same for all dates. 
% saved here:

% save('note_order.mat', 'notes')








%% file 7
%% extract the lfp reactions from 25 ms before to 110 ms
% after the time of the note. 
fs = 44100; 

load('ind&ampd_208.mat'); % for 1st version 
load('three_versions_208_note_amps_and_indexes_27_4.mat');
load('F271014-0010010085.mat', 'CRAW_006_KHz')

ind208_t  =ind208/fs;  % change indexes to seconds
ind208_t2  =ind208_t + 5;   % add the 5 seconds we took "extra" in the mean_lfp
%ind208_t3 = ind208_t2 * 9.856/10;   % the difference between the cttl times is 
        %not 10 seconds, it's 9.856, so here i normalize the timeline. 
ind1_208_t4 = ind208_t2 * 1000 * CRAW_006_KHz;

ind2_208_t  =ind2_208/fs; ind2_208_t2  =ind2_208_t + 5; 
%ind2_208_t3 = ind2_208_t2 * 9.856/10; 
ind2_208_t4 = ind2_208_t2 * 1000 * CRAW_006_KHz;

ind3_208_t  =ind3_208/fs; ind3_208_t2  =ind3_208_t + 5; 
%ind3_208_t3 = ind3_208_t2 * 9.856/10; 
ind3_208_t4 = ind3_208_t2* 1000 * CRAW_006_KHz;

ind4_208_t  =ind4_208/fs; ind4_208_t2  =ind4_208_t + 5; 
%ind4_208_t3 = ind4_208_t2 * 9.856/10; 
ind4_208_t4 = ind4_208_t2 * 1000 * CRAW_006_KHz;

% extract the reactions:
fs_ms = fs / 1000;
before = round(25 * fs_ms); % 25 ms space in index
after = round(110 * fs_ms);

for ii = 1:16
    ii
    reaction_peaks135 = zeros(round(before+after+1), 208,4);    %empty matrix to fill with reaction peaks
    load(['Ligeti_271014_channel_' num2str(ii) '4_lfp_matrixes_and_mean_lfp.mat'], 'mean_lfp');
    minus_mean_lfp = -mean_lfp;
    for jj = 1:4
        jj
        for kk = 1:length(eval(['ind' num2str(jj) '_208_t4']))
            index = round(eval(['ind' num2str(jj) '_208_t4(kk)']));
            start = index - before;
            finish = index + after ;
            section = minus_mean_lfp(start:finish,jj);
            reaction_peaks135(:,kk, jj) = section;
        end
    end
    save(['Ligeti_271014_channel_' num2str(ii) 'reaction_peaks.mat'], 'reaction_peaks135');
    clear( 'mean_lfp', 'reaction_peaks135');
end




%% file 8
%% plot the peaks of the reactions against the note amplitude,
% the surprise measure and the combination of the 2. 

% for plotting amps vs. reactions
load('note_order.mat')
% find frequency of each of the notes. 
freq_e = (sum(notes==1))/ 208; 
freq_f = (sum(notes ==2)) / 208; 
freq_hf = (sum(notes==3))/208; 
freq_g = (sum(notes==4))/208; 
load('three_versions_208_note_amps_and_indexes_27_4.mat', 'amps2_208', 'amps3_208', 'amps4_208')
load('ind&ampd_208.mat', 'amps208')
amps1_208 = amps208;
clear amps208
mid_notes = [notes(105:121),notes(144:176)];
amps = [amps1_208; amps2_208; amps3_208; amps4_208];
mid_amps = [amps(:,105:121), amps(:,144:176)];

% for plotting reactions against surprist measure 1 - -log(frequency of note)
surp_e1 = -log(freq_e); 
surp_f1 = -log(freq_f); 
surp_hf1 = -log(freq_hf); 
surp_g1 = -log(freq_g);
surps1 = [surp_e1, surp_f1, surp_hf1,surp_g1] 
surprise1_208 = zeros(1,208);
for ii = 1:208
    surprise1_208(ii) = surps1(notes(ii));
end
surprise_50 = [surprise1_208(105:121), surprise1_208(144:176)]; % this measure is the same for all versions

% for plotting reactions against combination of surprise and amplitude
comb = zeros(4,208);
for ii = 1:4
    for jj = 1:208
    comb(ii, jj) = surprise1_208(jj) * amps(ii, jj);  
    end
end

comb_50 = [comb(:,105:121), comb(:,144:176)];

        

for ii = 1:16
    load(['Ligeti_271014_channel_' num2str(ii) 'reaction_peaks.mat']);
    peaks = zeros(4,208);
    for jj = 1:4
        % plot peaks against amplitude
        peaks(jj,:) = max(reaction_peaks135(:,:,jj));
        mid_peaks = [peaks(:,105:121), peaks(:,144:176)];
        figure
        for kk = 1:length(mid_notes)
               if mid_notes(kk) == 4
                   col = [1 0 0];
               else
                   col = [0 1 0];
               end
               scatter(mid_amps(jj,kk), mid_peaks(jj,kk), [], col, 'filled')
               hold on
        end
        xlabel('amplitude')
        ylabel('reaction peak')
        r_amp_only = corrcoef (mid_amps(jj,:), mid_peaks(jj,:))  % 0.2896
        title(['correlation = ' num2str(r_amp_only(1,2)) '     channel ' ...
            num2str(ii) ' version ' num2str(jj)])
        savefig(['Ligeti_271014_peaks_vs_amps_channel_' num2str(ii) '_version_' num2str(jj)])
        close all
        %plot peaks against surprise measure 1: -log(p) of note
        figure
        for kk = 1:length(mid_notes)
            if mid_notes(kk) == 4
                col = [1 0 0];
            else
                
                col = [0 1 0];
            end
            plot(surprise_50(kk), mid_peaks(jj, :), '.', 'color', col)
            hold on
        end
        xlabel('surprise measure 1 (-log(probability of sound))')
        ylabel('reaction peak')
        r_prob_only = corrcoef(surprise_50, mid_peaks(jj,:))  %0.4716
        title(['correlation = ' num2str(r_prob_only(1,2))  '     channel ' ...
            num2str(ii) ' version ' num2str(jj)])
        savefig(['Ligeti_271014_peaks_vs_surp1_channel_' num2str(ii) '_version_' num2str(jj)])
        close all
        % plot against the combination of surprise and sound level
        for kk = 1:length(mid_notes)
            if mid_notes(kk) == 4
                col = [1 0 0];
            else
                
                col = [0 1 0];
            end
            plot(comb_50(jj,kk), mid_peaks(jj,kk), '.', 'color', col, 'markersize',20)
            hold on
        end
        r_comb = corrcoef(comb_50(jj,:), mid_peaks(jj,:))  %0.5033
        xlabel('surprise measure 1 (-log(probability of sound)) * amplitude')
        ylabel('reaction peak')
        title(['correlation = ' num2str(r_comb(1,2)) '    channel ' ...
            num2str(ii) ' version ' num2str(jj)])
        savefig(['Ligeti_271014_peaks_vs_amp_times_surp1_combination_channel_' num2str(ii) '_version_' num2str(jj)])
        close all

    end
end






%% file 9
%% make the following graph to validate the data:

% 3. plot each of the lfp mean responses on top of the ligeti signal, and
% on this put the peaks of the notes. right now it seems there is a
% constant 100-ms difference between the times and the reaction. 

% done for channel 6 only. 



%% plot each of the lfp mean responses on top of the ligeti signal, and
% on this put the peaks of the notes. right now it seems there is a
% constant 100-ms difference between the times and the reaction. 
load('F271014-0010010085.mat', 'CRAW_006_KHz')
load('Ligeti_271014_channel_64_lfp_matrixes_and_mean_lfp.mat', 'mean_lfp');
%load('indexes_new_t3');
load('note_order.mat')

minus_mean_lfp = -mean_lfp;
 % the lfp's sampling freq
lfp_timeline = ((1:length(minus_mean_lfp))-1)/CRAW_006_KHz/1000;
plot(lfp_timeline,minus_mean_lfp(:,1))
tremolo = 122:143;
hold on
for ii = 1:length(ind208_t2)
     if notes(ii) == 4
         if any(tremolo == ii)
            line([ind208_t2(ii), ind208_t2(ii)], [-1500, 3500], 'color', 'y')
         else
        line([ind208_t2(ii), ind208_t2(ii)], [-1500, 3500], 'color', 'r')
         end
     else
         line([ind208_t2(ii), ind208_t2(ii)], [-1500, 3500], 'color', 'g')    
     end
end
title('mean response and note times in version 1 - original')

hold on
load('three_versions_208_note_amps_and_indexes_27_4', 'fsig1')
fsig1_big = fsig1*5000;
fs = 44100;
fsig1_ts = ((1:length(fsig1))-1)/fs;

fsig1_ts5 = fsig1_ts + 5; %add 5 seconds of silence to the beginning of the ligeti signal
%fsig1_ts5 = fsig1_ts5 *9.856/10;

plot(fsig1_ts5, fsig1_big,'m')
savefig('Ligeti_271014 mean response and note times in version 1 - original')
close all
%%
%version 2
plot(lfp_timeline,minus_mean_lfp(:,2))
tremolo = 122:143;
hold on
for ii = 1:length(ind2_208_t2)
     if notes(ii) == 4
         if any(tremolo == ii)
            line([ind2_208_t2(ii), ind2_208_t2(ii)], [-1500, 3500], 'color', 'y')
         else
        line([ind2_208_t2(ii), ind2_208_t2(ii)], [-1500, 3500], 'color', 'r')
         end
     else
         line([ind2_208_t2(ii), ind2_208_t2(ii)], [-1500, 3500], 'color', 'g')    
     end
end
title('mean response and note times in version 2- ident')

hold on
load('three_versions_208_note_amps_and_indexes_27_4', 'fsig2')
fsig2_big = fsig2*5000;
fsig2_ts = ((1:length(fsig2))-1)/fs;
fsig2_ts5 = fsig2_ts + 5; %add 5 seconds of silence to the beginning of the ligeti signal
%fsig2_ts5 = fsig2_ts5 *9.856/10;
plot(fsig2_ts5, fsig2_big,'m')
savefig('Ligeti_271014 mean response and note times in version 2- ident')
close all
%%
%version 3
plot(lfp_timeline,minus_mean_lfp(:,3))
tremolo = 122:143;
hold on
for ii = 1:length(ind3_208_t2)
     if notes(ii) == 4
         if any(tremolo == ii)
            line([ind3_208_t2(ii), ind3_208_t2(ii)], [-1500, 3500], 'color', 'y')
         else
        line([ind3_208_t2(ii), ind3_208_t2(ii)], [-1500, 3500], 'color', 'r')
         end
     else
         line([ind3_208_t2(ii), ind3_208_t2(ii)], [-1500, 3500], 'color', 'g')    
     end
end
title('mean response and note times in version 3- all up')

hold on
load('three_versions_208_note_amps_and_indexes_27_4', 'fsig3')
fsig3_big = fsig3*5000;
fsig3_ts = ((1:length(fsig3))-1)/fs;
fsig3_ts5 = fsig3_ts + 5; %add 5 seconds of silence to the beginning of the ligeti signal
%fsig3_ts5 = fsig3_ts5 *9.856/10;
plot(fsig3_ts5, fsig3_big,'m')
savefig('Ligeti_271014 mean response and note times in version 3- all up')
close all
%%
%version 4
plot(lfp_timeline,minus_mean_lfp(:,4))
hold on
for ii = 1:length(ind4_208_t2)
     if notes(ii) == 4
         if any(tremolo == ii)
            line([ind4_208_t2(ii), ind4_208_t2(ii)], [-1500, 3500], 'color', 'y')
         else
        line([ind4_208_t2(ii), ind4_208_t2(ii)], [-1500, 3500], 'color', 'r')
         end
     else
         line([ind4_208_t2(ii), ind4_208_t2(ii)], [-1500, 3500], 'color', 'g')    
     end
end
title('mean response and note times in version 4- e,g inverted')

hold on
load('three_versions_208_note_amps_and_indexes_27_4', 'fsig4')
fsig4_big = fsig4*5000;
fsig4_ts = ((1:length(fsig4))-1)/fs;
fsig4_ts5 = fsig4_ts + 5; %add 5 seconds of silence to the beginning of the ligeti signal
%fsig4_ts5 = fsig4_ts5 *9.856/10;
plot(fsig4_ts5, fsig4_big,'m')

savefig('Ligeti_271014 mean response and note times in version 4- e,g inverted')


%% these graphs are good. 