function MI = mutualInformation(signal1, signal2)
% calcutlate the mutual information between signals
% signal = [timeBins, channels]
% MI = [signal1Channels, signal2Channels]
spikeLength = length(signal1);
p_s1_0   = sum(signal1   == 0) / spikeLength;
p_s1_1   = sum(signal1   == 1) / spikeLength;
p_s2_0 = sum(signal2 == 0) / spikeLength;
p_s2_1 = sum(signal2 == 1) / spikeLength;

p_00 = (signal1==0)'*(signal2==0) / spikeLength;
p_01 = (signal1==0)'*(signal2==1) / spikeLength;
p_10 = (signal1==1)'*(signal2==0) / spikeLength;
p_11 = (signal1==1)'*(signal2==1) / spikeLength;

MI = p_00.*log2(p_00./(p_s1_0'*p_s2_0)) + ...
     p_01.*log2(p_01./(p_s1_0'*p_s2_1)) + ...
     p_10.*log2(p_10./(p_s1_1'*p_s2_0)) + ...
     p_11.*log2(p_11./(p_s1_1'*p_s2_1));
end
