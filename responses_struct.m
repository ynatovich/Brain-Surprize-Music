% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This program build responses struct represent a response for each note. 
% Each note response represented with:
% x: The range of time -50 ms till +150 ms from where the note started
% y: The range of response electrodes -50 ms till +150 ms from where the
%    note started
% note number: the number of the note (1-4)
% isSurprize: 1 if the note is 4, 0 else. Represent whether the note is the
%             surprized note.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/Ligeti271014_2_mean/channel_6_mean_lfp.mat")
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/note_order.mat")
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/ind_ampd_208.mat", 'ind208')
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/t.mat")
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/ts.mat")
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/responses_ind.mat")
load("/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/matrix_data/fahdata.mat")
N = 208;
zerosfill = zeros(100, 1)';
dates = ['040514', '051114', '110514', '150714_2', '271014', '271014_2'];
% for date=dates
    date = '271014_2';
    mkdir(sprintf('Ligeti%s/responses', date))
    for channeli=1:16
        load(['/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti' date '/Ligeti271014_2_mean/channel_' num2str(channeli) '_mean_lfp.mat']);
        mkdir(sprintf('Ligeti%s/responses/channel%d', date, channeli))
        oldfile = cd(sprintf('Ligeti%s/responses/channel%d', date, channeli));
        resp = mean_1(1, :);
        consts = repmat(struct('x', [],  'xfah', [], 'fahdata', []), N, 1 );
        for ii=1:N
            consts(ii).note = notes(ii);
            consts(ii).x = t(responses_ind_t(ii)-50:responses_ind_t(ii)+150);
            consts(ii).notetime = ts(uint64(ind208(ii)));
            consts(ii).maxvolume = fahdata(uint64(ind208(ii)))*2000+1000;
            if ii>1
                consts(ii).xfah = ts(responses_ind_ts(ii)-2000:20:responses_ind_ts(ii)+5000);        
                consts(ii).fahdata = fahdata(responses_ind_ts(ii)-2000:20:responses_ind_ts(ii)+5000)';
                consts(ii).isSurprize = (((notes(ii) == 4)&notes(ii-1)~=4)*true + (notes(ii) ~= 4)*false);
            else
                consts(ii).xfah = [zerosfill ts(1:20:responses_ind_ts(ii)+5000)];
                consts(ii).fahdata = [zerosfill fahdata(1:20:responses_ind_ts(ii)+5000)'];
                consts(ii).isSurprize = ((notes(ii) == 4)*true + (notes(ii) ~= 4)*false);
            end
        end
        save('consts.mat','consts')
        %%
        % this section specific calculate the fahdata for eac channel and version
        versions = repmat(struct('version', []), 4, 1 );
        for veri=1:4
            versions(veri).version = repmat(struct('repeat', []), 6, 1 );
            for meani=1:6
                resp = eval(['mean_' num2str(veri) '(meani, :)']);
                versions(veri).version(meani).repeat = repmat(struct('responses_one_note', []), N, 1 );
                for ii=1:N
                    versions(veri).version(meani).repeat(ii).responses_one_note = resp(responses_ind_t(ii)-100:responses_ind_t(ii)+250);
                end
            end
        end
        save('responses.mat','versions')
        cd(oldfile);
    end 
% end
