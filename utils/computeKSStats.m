function [ Z, U, xAxis, KSSorted, ks_stat ] = computeKSStats( spikeTrain, lambda, opt )
%COMPUTEKSSTATS Compute related KS statistics.
%   Input: 
%       spikeTrain    - 1 * #timeslot
%       lambda        - 1 * #timeslot
%       opt           - need opt.sampleRate, opt.DTCorrelation
%   Output:
%       Z: rescaled spike times
%       U: Zjs tranformed to be uniform(0,1) 
%       xAxis: x-axis of the KS plot
%       KSSorted: y-axis of KS plot
%       ks_stat: the KS statistic. Maximum deviations from the 45
%               degree line for each conditional intensity function.

if (opt.DTCorrelation == 1)
    pk = lambda .* (1/opt.sampleRate);
    pk = max(pk, 1e-10);
    minDim = min(length(spikeTrain), length(lambda));
    pk = pk(1:minDim);
    spikeTrain = spikeTrain(1:minDim);
    
    pk = nanmin(nanmax(pk, 0), 1);
    intValues = ksdiscrete(pk, spikeTrain, 'spiketrain');
    intValues = intValues';
else 
    pk = lambda .* (1/opt.sampleRate);
    pk = max(pk, 1e-10);
    minDim = min(length(spikeTrain), length(lambda));
    pk = pk(1:minDim);
    spikeTrain = spikeTrain(1:minDim);
    
    spikes = find(spikeTrain ~= 0);
    spikes = [1 spikes'];
    intValues = zeros(1, length(spikes)-1);
    
    for spikeIdx = 1:length(spikes) -1
        accum = sum(pk(spikes(spikeIdx):spikes(spikeIdx+1)));
        intValues(1, spikeIdx) = accum;
    end
end

Z = intValues;
U = 1 - exp(-Z);
KSSorted = sort(U, 'ascend');
N = size(KSSorted, 2);
if (N ~= 0)
    xAxis = (([1:N]-.5)/N);
    ks_stat = max(abs(KSSorted - (([1:N]-.5)/N)));
end

end


function [rst,varargout] = ksdiscrete(pk,st,spikeflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ksdiscrete.m 
% written by Rob Haslinger, December 2009
%
% This function performs time rescaling of ISIs based upon the discrete
% time version of the time rescaling theorem as described in Haslinger,
% Pipa and Brown (2010?).  This method corrects for biases in the KS plot 
% caused by the temporal discretization.
%
% This function can be called in two ways
%
% 1) input the discrete time conditional probabilities "pk"  where 0<=pk<= 1
% and the spike train "spiketrain" which has elements either equal to 0 (no
% spike) or 1 (spike). There is also a flag 'spiketrain' to indicate that
% it is the full spike train.
%
% [rst,rstsort,xks,cb,rstoldsort] = ksdiscrete(pk,spiketrain,'spiketrain')
%
% 2) input the discrete time conditional probabilities "pk" and a list of
% the indicies "spikeind" of the bin indicies that the spikes are locaed in. 
% There is also a flag 'spikeind' to indicate that the indicies are
% being given, not the full spike train
%
% [rst,rstsort,xks,cb,rstoldsort] = ksdiscrete(pk,spikeind,'spikeind');
%
% required output:
%
% rst : a vector of unsorted uniformly distributed rescaled times. This is
% the only output that is required.
% 
% optional output, given in the order they appear in the list function
% outputs :
%
% rstsort : a vector of rescaled times sorted into ascending order
% xks : a vector of x axis values to plot the sorted rescaled times against
% cb : the value of the 95% confidence bounds
% rstoldsort : a vector of sorted rescaled times done without the discrete
% time correction
%
% To make a KS plot one would do
% plot(xks,rstsort,'k-');
% hold on;
% plot(xks,xks+cb,'k--',xks,xks-cb,'k--');
%
% To make a Differential KS plot one would do
% plot(xks,rstsort-xks,'k-');
% hold on;
% plot(xks,zeros(length(xks))+cb,'k--',xks,zeros(length(xks))-cb);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with determining the inputs and some basic input error checking

    if nargin < 3 || nargin > 3;  
        error('Number of input arguments must be equal to 3'); 
    end;

    % make pk into a column vector;

    [m1,m2]=size(pk);
        if (m1 ~=1 && m2 ~=1); error('pk must be a vector'); end;
        if (m2>m1); pk=pk'; end;
        [m1,m2]=size(pk);

    % make sure pk's are within [0,1]
    index=find(pk<0);
    if isempty(index) ~=1; 
        error('all values for pk must be within [0,1]'); 
    end;
    index=find(pk>1);
    if isempty(index) ~=1; 
        error('all values for pk must be within [0,1]'); 
    end;
    clear index;    

    % make column vector of spike indicies

    if strcmp(spikeflag,'spiketrain'); % spike train input

        [n1,n2]=size(st);
          if (n1 ~=1 && n2 ~=1); error('spike train must be a vector'); end;
        if (n2>n1); st=st'; end;

        if m1 ~= n1; error('pk and spike train must be same length'); end;

        spikeindicies=find(st==1);

        Nspikes=length(spikeindicies);

    elseif strcmp(spikeflag,'spikeind'); % spike index input

        [n1,n2]=size(st);
          if (n1 ~=1 && n2 ~=1); error('spike indicies must be a vector'); end;
        if (n2>n1); st=st'; end;

        spikeindicies=unique(st);
        Nspikes=length(spikeindicies);

    end;

    % check that those indicies are in [1:length(pk)];

    if(isempty(spikeindicies))
    	rst = pk;
        return;
    end
    if spikeindicies(1)<1; 
         error('There is at least one spike with index less than 0'); 
    end;
    if spikeindicies(Nspikes)>length(pk); 
         error('There is at least one spike with a index greater than the length of pk'); 
    end;    

    % error checking done

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now do the actual discrete time KS test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % initialize random number generator
    s = RandStream('mt19937ar','Seed', sum(100*clock));
    RandStream.setGlobalStream(s);
    %rand('twister',sum(100*clock));

    % make the qk's

    qk=-log(1-pk);

    % make the rescaled times

    rst=zeros(Nspikes-1,1);
    rstold=zeros(Nspikes-1,1);

    for r=1:Nspikes-1;

        total = 0;

        ind1=spikeindicies(r);
        ind2=spikeindicies(r+1);

        total=total+sum(qk(ind1+1:ind2-1));

        delta=-(1/qk(ind2))*log(1-rand()*(1-exp(-qk(ind2))));

        total=total+qk(ind2)*delta;

        rst(r)=total;

        rstold(r)=sum(qk(ind1+1:ind2));

    end;

%     rst=1-exp(-rst);
%     rstold=1-exp(-rstold);

    % optional outputs

    rstsort=sort(rst);
    varargout{1}=rstsort;

    inrst=1/(Nspikes-1);
    xrst=(0.5*inrst:inrst:1-0.5*inrst)';
    varargout{2}=xrst;

    cb=1.36*sqrt(inrst);
    varargout{3}=cb;    

    varargout{4}=sort(rstold);
end