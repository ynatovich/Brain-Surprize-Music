% make histogram of distance between the note and the react

load('ind_ampd_208.mat', 'ind208')
load('t.mat')
delta = [];
note_time = (ind208'/44.1+5000)';
for i=1:4
    m = sprintf('figures/response_mean_%d.mat', i);
    x = load(m, 'responses_values');
    delta = [delta note_time - t(responses_ind)];
end 
h = histogram(delta, 100)
