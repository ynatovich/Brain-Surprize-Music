load('ind_ampd_208.mat', 'ind208')
hold all
note_time = (ind208'/44.1+5000);
for i=1:4
    mean = sprintf('mean_%d', i);    
    responses_ind = zeros(length(note_time), 1)';
    m = load('Ligeti271014_2_mean/channel_6_mean_lfp.mat', mean);
    y = m.(mean)(6,:);
    t = (1:length(y))/2.2;
    for note=1:length(note_time)
        % find local minimum represent the responses for one note
        s = find(t >= max(1, note_time(note) - 50), 1);
        e = find(t <= min(note_time(note) + 50, note_time(end)), 1, 'last');
        [x, ind] = min(y(s:e));
        % if the range is not enougth, continue to the local mini
        if and(ind > 1, ind < note_time(end))
            if and(ind == s, y(ind-1) < y(ind)) 
            while y(ind-1) < y(ind)
                ind = ind-1;
            end
            end 
            if and(ind == e, y(ind+1) < y(ind))
                while y(ind+1) < y(ind)
                    ind = ind+1;
                end
            end 
        end
        responses_ind(note) = ind + s - 1;
    end
    subplot(2, 2, i)
    responses_values = y(responses_ind);
    hold all
    line((note_time)*[1 1],[-4000 1500],'col','r')
    plot(t, y, t(responses_ind),y(responses_ind),'y*')
    tit = sprintf('mean of version %d', i);
    title(tit)
    hold off

%     file = sprintf('response_mean_%d.mat', i);
%     save(file, 'responses_ind', 'responses_values');

end
