% % plot(t, res)
% plot(ts, fahdata*2000+1000)
% hold on
% for ii=1:N
%     if responses(ii).isSurprize
%         plot(responses(ii).x, responses(ii).res, 'b')
%         plot(responses(ii).xfah, responses(ii).fahdata*2000+1000, 'r')
%         continue
%     end 
%     plot(responses(ii).x, responses(ii).res, 'y')
%     plot(responses(ii).xfah, responses(ii).fahdata*2000+1000, 'g')
% end
% 
% hold all
% plot(ts,fahdata*2000+1000)
% scatter(ts(uint64(ind208)), fahdata(uint64(ind208))*2000+1000, 'filled')
allchannels = repmat(struct('channel', []), 16, 1);
for ii=1:16
    load(['/GoodmanHome/global/ligeti/all_trials_analyzed/Ligeti271014_2/responses/channel' num2str(ii) '/responses.mat']);
    allchannels(ii).channel = versions;
end
save('allresponses271014_2.mat', 'allchannels') 