% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This program find where is the response of each ligeti`s notes and
% plot it.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% t = (1:length(mean_1))/2.2
% ts=(1:length(data))/44.1+5000

% load('ind_ampd_208.mat', 'ind208')
note_time = (ind208'/44.1+5000);
responses_ind = zeros(length(note_time), 1)';
% load('Ligeti271014_2_mean/channel_6_mean_lfp.mat', 'mean_1');
figure
for i=1:6
    y = mean_1(i,:);
    t = (1:length(y))/2.2;
    for note=1:length(note_time)
        % find local minimum represent the responses for one note
        s = find(t >= max(1, note_time(note) - 50), 1);
        e = find(t <= min(note_time(note) + 50, note_time(end)), 1, 'last');
        [x, ind] = min(y(s:e));
        % if the range is not enougth, continue to the local mini
        ind = ind + s - 1;
        if and(ind == s, ind > 1)
            while y(ind-1) < y(ind)
                ind = ind-1;
            end
        end 
        if and(ind == e, ind < length(y))
            while y(ind+1) < y(ind)
                ind = ind+1;
            end
        end
        responses_ind(note) = ind;
    end
    subplot(2, 3, i)
    responses_values = y(responses_ind);
    hold all
    line((note_time)*[1 1],[-4000 1500],'col','r')
    plot(t, y, t(responses_ind),y(responses_ind),'y*')

    hold off
end


for i=1:6
    y = mean_1(i,:);
    t = (1:length(y))/2.2;
    for note=1:length(note_time)
        % find local minimum represent the responses for one note
        s = find(t >= max(1, note_time(note) - 50), 1);
        e = find(t <= min(note_time(note) + 50, note_time(end)), 1, 'last');
        [x, ind] = min(y(s:e));
        % if the range is not enougth, continue to the local mini
        ind = ind + s - 1;
        if and(ind == s, ind > 1)
            while y(ind-1) < y(ind)
                ind = ind-1;
            end
        end 
        if and(ind == e, ind < length(y))
            while y(ind+1) < y(ind)
                ind = ind+1;
            end
        end
        responses_ind(note) = ind;
    end
    responses_values = y(responses_ind);
    figure(i)
    hold on
    line((note_time)*[1 1],[-4000 1500],'col','r')
    plot(t, y, t(responses_ind),y(responses_ind),'y*')
    hold off
end
t=(1:length(mean_1))/2.2;
data=data(85020:end);
ts=(1:length(data))/44.1+5000;
hold off
plot(t, mean_1(6,:))
hold all
line((ind208'/44.1+5000)*[1 1],[-4000 1500],'col','r')
plot(ts,fahdata*2000+1000)
clf
plot(t, mean_1(6,:))
line((ind208'/44.1+5000)*[1 1],[-4000 1500],'col','r')

%%
load('ind_ampd_208.mat', 'ind208')
note_time = (ind208'/44.1+5000);
responses_ind = zeros(length(note_time), 1)';
% responses_ind_far = [];
load('Ligeti271014_2_mean/channel_6_mean_lfp.mat', 'mean_1');
figure;clf;
y = mean_1(1,:);
t = (1:length(y))/2.2;
for note=1:length(note_time)
    % find local minimum represent the responses for one note
    s = find(t >= max(1, note_time(note) - 50), 1);
    e = find(t <= min(note_time(note) + 50, note_time(end)), 1, 'last');
    [x, ind] = min(y(s:e));
    % if the range is not enougth, continue to the local mini
    ind = ind + s - 1;
    if and(ind == s, ind > 1)
        while y(ind-1) < y(ind)
            ind = ind-1;
        end
    end 
    if and(ind == e, ind < length(y))
        while y(ind+1) < y(ind)
            ind = ind+1;
        end
    end
%     if or(ind-s < 20, e-ind < 20)
%         responses_ind_far = [responses_ind_far ind];
%     end 
    responses_ind(note) = ind;
end
responses_values = y(responses_ind);
hold on
line((note_time)*[1 1],[-4000 1500],'col','r')
% plot(t, y, t(responses_ind),y(responses_ind),'g*')
plot(t, y, t(responses_ind_far),y(responses_ind_far),'g*')
hold off