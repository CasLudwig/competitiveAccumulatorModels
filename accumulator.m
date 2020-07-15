function [choice,rt]=accumulator(inp,s,sv,parms)

% Integrates the internal responses over time. It receives two vectors that
% determine the mean and variability of the start points: 's' is a 1 x M
% vector with the mean starting point of the target and non-targets; 'sv'
% is a 1 x M vector with half the range around those means. The parameters
% of the accumulation process are: beta (lateral inhibition, set to 0 for
% an independent race); theta (threshold); mu_nd (mean of Gaussian
% non-decisional delay) and sigma_nd (standard deviation of Gaussian
% non-decisional delay).
% 
% Bear in mind that some parameters will scale together, so it is customary
% to fix one of them. For instance, an increase in the gain of the internal
% response may be compensated for by scaling the variability in the gain 
% and the response threshold. In this particular application I have
% fixed the noise within the accumulation process. For ballistic
% accumulation simply set intvar to 0, but in that case I recommend fixing
% the threshold (at some arbitrary value, e.g. 100).

global h;

beta=parms(1); %strength of lateral inhibitory interactions
theta=parms(2); %response threshold
t0=parms(3); %non-decisional delay
t0_unif=parms(4); %variability in the non-decisional delay (0=no, 1=yes)

N=size(inp,1);
M=size(inp,3);
t=0:h:h*(size(inp,2)-1);
intvar=0; %variance of the noise added to the accumulator (i.e. diffusion noise); I have fixed this (but see comments at the top)
intsd=sqrt(intvar*h); %sd of Gaussian noise scaled with the time step

choice=zeros(N,2);
rt=repmat(NaN,N,1);

s=repmat(s,N,1);s=reshape(s,N,1,M);
sv=repmat(sv,N,1);sv=reshape(sv,N,1,M);
allstartpoints=(s-sv)+rand(N,1,M).*(2*sv); %N x 1 x M matrix of uniformly distributed start points
noisy_inp=inp+normrnd(0,intsd,N,length(t),M); %noise added to accumulator
noisy_int=zeros(size(noisy_inp));
noisy_int(:,1,:)=allstartpoints+noisy_inp(:,1,:);
noisy_int(noisy_int(:,:,:)<0)=0; %rectification
competitors=zeros(N,1); %use this vector to keep track of the sum of the activations of competing accumulators
%now we start integrating
for i=2:length(t) %i indexes the time step
    for j=1:M %run through each accumulator in turn: j indexes the current accumulator of interest
        for k=1:M %we first determine the sum of the competing accumulators (i.e. k indexes the competing accumulators)
            if(k~=j) %this defines what is a competitor
                competitors=competitors+noisy_int(:,i-1,k);
            end
        end
        noisy_int(:,i,j)=noisy_int(:,i-1,j)+noisy_inp(:,i,j)-beta*competitors; %e.g. see equation 3 in Usher & McClelland (2001)
        noisy_int(noisy_int(:,:,:)<0)=0; %rectification
        competitors(1:N,1)=0; %re-set to 0 for the next response alternative or, indeed, time step
    end
end

%now obtain threshold crossing data
logmatrix=noisy_int>=theta; %logical matrix with binary numbers to indicate whether a value is greater than threshold
logmatrix=double(logmatrix); %turn into double so we can take cumulative sum
logmatrix=cumsum(logmatrix,2);
logmatrix(logmatrix==0)=NaN;
[vals,ind]=min(logmatrix,[],2); %both vals and ind are Nx1xM matrices: ind gives specifies the time of the first threshold crossing for each accumulator
ind(isnan(vals))=t(end)+1; %processes that never reached threshold are given a time greater than the maximum allowed
allrts=zeros(N,M);
for i=1:M %now put all threshold crossing indices in one 2D matrix
    allrts(:,i)=ind(:,1,i)*h;
end
[rtvals,rtinds]=min(allrts,[],2); %rtvals gives the winning decision times; rtinds gives the winning accumulators (1...M)
rtvals(rtvals>=t(end))=NaN; %if the winning RT is greater than the maximum allowed, it was a process that never reached the threshold
rtinds(isnan(rtvals))=NaN;
choice(:,1)=rtinds;
choice(:,2)=rtinds;
%first column of the choice matrix codes the accuracy of the response: if the first accumulator did NOT win, the observed response was an error
choice(choice(:,1)>1,1)=0;

%deal with possible ties - only if the tie is between the winning rt and one or
%more other accumulators (otherwise we don't care)
tie=zeros(N,M+1);
for i=1:M
    tie(:,i)=rtvals==allrts(:,i);
end
tie(:,M+1)=sum(tie,2); %final column now simply indicates whether there is a tie (values > 1) or not (1: the match of the winning rt with itself)
tie(tie(:,M+1)==1,M+1)=0; %turn final column into a "logical" identifier
tie(tie(:,M+1)>1,M+1)=1;
tieinds=find(tie(:,M+1)); %indices to the trials with a tie
if ~isempty(tieinds) %if there are ties, deal with each one in turn
    for i=1:length(tieinds)
        tindex=allrts(tieinds(i),choice(tieinds(i),2))/h; %column index to the noisy_int matrix
        %we resolve the tie by choosing the accumulator with the largest value at the time of the winning RT!
        [maxval,maxind]=max(noisy_int(tieinds(i),tindex,:));
        if(maxind==1)
            choice(tieinds(i),1:2)=1; %target accumulator has the largest value
        else choice(tieinds(i),1:2)=[0 maxind]; %competitor has the largest value
        end
    end
end

%the rt is a combination of the decision time to threshold and a
%non-decisional delay that is either constant (t0) or uniformly
%distributed (in which case t0 is the upper limit)
if t0_unif==1
    rt=rtvals+t0*rand(N,1); %add uniformly distributed non-decisional component to the winning rt
else rt=rtvals+t0; %add constant non-decisional delay
end

%plot the average accumulator paths - again, you will probably want to
%comment this out at some point
figure;hold on;
lines={'k-','b-','y-','r-','g-'}; %assuming we will never use more than 5 alternatives
for i=1:M
    plot(t,mean(noisy_int(:,:,i)),lines{i},'linewidth',2);
end
set(gca,'FontSize',12,'PlotBoxAspectRatioMode','manual','XLim',[0 500]); %just show the first 500 ms
xlabel('time (ms)');
ylabel('activity (arbitrary units)');