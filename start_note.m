% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This file detect the start time of each note of ligeti. %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

hold all
% load('note_time_volume.mat')
% load('ind_ampd_208.mat', 'ind208')
% load('ts.mat')
% load('t.mat')
note_time = (ind208'/44.1+5000);
% plot(ts,data*2000+1000)
% hdata=hilbert(data);
% ahdata=abs(hdata);
b=fir1(1200,100/44100*2);
% fahdata=filtfilt(b,1,ahdata);
% plot([data(:) ahdata(:) fahdata(:)]) % plot(ts,hdata+1000)
% plot([data(:) ahdata(:)])
responses_ind = zeros(length(note_time), 1)';
for note=1:length(note_time)
    % find local minimum represent the responses for one note
    s = find(ts >= max(1, note_time(note) - 100), 1);  % find min x axis 
    e = find(ts <= note_time(note), 1, 'last');
    maxval = min(fahdata(s:e));
    ind = find(fahdata(s:e)==maxval,1,'last');
    % if the range is not enougth, continue to the local mini
    ind = ind + s - 1;
    responses_ind(note) = ind;
end
plot(ts,fahdata*2000+1000)
line((note_time)*[1 1],[500 1500],'col','g')
line((responses_ind'/44.1+5000)*[1 1],[500 1500],'col','r')

% plot(ts, fahdata, ts(responses_ind), fahdata(responses_ind), 'r*')
% plot(t, fahdata, t(responses_ind), fahdata(responses_ind), 'r*')
%%
responses_ind_t = zeros(length(note_time), 1)';
for itime=1:length(responses_ind)
    % find local minimum represent the responses for one note
    responses_ind_t(itime) = find(t >= ts(responses_ind(itime)), 1);
end





%%
% hold on
% plot(ts,fahdata*2000+1000)
% line((note_time)*[1 1],[500 1500],'col','g')
% line((responses_ind'/44.1+5000)*[1 1],[500 1500],'col','r')
% hold off