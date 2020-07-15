function [y]=compute_fit(observed,model,varargin)
% Takes a vector of observed latencies and a vector of simulated latencies
% and computes some measure of fit: How well does the simulated
% distribution approximate the observed data? To call, type for example
% 'y=compute_fit([observed distribution], [simulated
% distribution],'deviance',[N1 N2],[quantiles]);
% 
% The last 3 input arguments are, in turn:
% 
% Argument 3: The relevant metric to compute.
% 1. Deviance: -2 times the log-likelihood; this measure can be
% computed for binned (quantiles) and unbinned data. The quantile version
% is chosen when the user supplies a set of quantiles as a fifth argument.
% 2. Chi-square: can only be computed on binned data, so requires the user
% to specify a set of quantiles at which to evaluate the cumulative
% probability distribution.
% 3. G-square: another commonly used metric for binned data. As with
% Chi-square, the user will need to specify the set of quantiles.
% 
% Argument 4: the total NUMBER of observed and simulated trials (1 x 2
% vector). These may not necessarily correspond to the length of the
% 'observed' and 'model' vectors, for instance, when you specifically want
% to fit the proportion of correct responses. In that case, the model fits
% a defective latency distribution: one that integrates to the probability
% of the response (i.e. 0.8 if the error proportion was 0.2). Implicitly
% then, we fit the latency distribution of the dominant (e.g. correct)
% response type and the overall proportion of the alternative (e.g. error)
% response type. If the frequencies supplied in this vector do match the
% length of the latency vectors, the distributions will integrate to 1 and
% you would only be fitting the latency distribution of the dominant
% response.
% 
% Optional argument 5: quantile vector. If and when the distributions need
% to be binned (essentially a comparison of the observed and predicted
% CDFs, or distribution functions), you need to specify a vector of
% quantiles as a fifth argument. Reasonable quantiles are, for example,
% [0.1; 0.3; 0.5; 0.7; 0.9]. Actually, for simulated data, working with
% quantiles is preferable because the alternative involves the additional
% step of estimating a full probability density function from the simulated
% data (through a kernel density estimator).
% 
% NOTE 1: if you decide to work with the full probability density (i.e.
% unbinned data) the goodness-of-fit metric will always be the deviance,
% regardless of the third argument.
% 
% NOTE 2: If this function returns some strange value (e.g. -Inf or NaN), it
% is likely that there is something wrong with the model parameters. That
% is, if the predicted distribution is way off the empirical one, this will
% cause problems for all measures of fit used here. For instance, for
% Chi-square, the usual caveat applies in that you need to make sure that
% the predicted frequency for each bin is at least 5. For the deviance and
% G-square measures you need to make sure that you are not trying to take
% the logarithm of 0 or negative numbers.
% 
% NOTE 3: The input distributions will be monitored for the presence of
% NaNs. For example, simulated trials that return NaN correspond to those
% trials in which the decision threshold was never reached. The presence of
% NaNs and their subsequent removal in the code, means that the length of
% the latency vector may not match the specified number of trials in the N
% vector. As a result, the code below would fit a defective distribution.
% If you do not want this, you need to make sure that either you do not
% include the NaN trials in your specified trial numbers (the N vector).
% You do not have to filter the NaN trials out from the latency vectors,
% the code below does this for you.

CDF=0; %default: assume we are working with the whole distribution
%Make sure latency vectors are columns
if(size(observed,2)>size(observed,1))
    observed=observed';
end    
if(size(model,2)>size(model,1))
    model=model';
end
observed=observed(~isnan(observed)); %filter out possible NaNs - see Note 3
model=model(~isnan(model));

%Now deal with the various input arguments
if(size(varargin,2)>=2)
    metric=varargin{1}; %which metric to compute?
    N=varargin{2}; %total frequencies of observed and simulated responses
    defect_p=length(observed)/N(1); %probability of the observed responses
else
    error('Error: insufficient input arguments! Type ''help compute_fit'' for more information');
end
if(size(varargin,2)==3)
    CDF=1;
    q=varargin{3};
    %make sure q is a column vector
    if(size(q,2)>size(q,1))
        q=q';
    end
end

if(CDF==1) %Compare distribution functions
    %Determine the quantile latencies for the observed data. Make sure the
    %quantile latencies are in a column vector.
    qobs=quantile(observed,q);
    qobs(:,2)=q;
    qobs(end+1,1:2)=[inf 1]; %Final quantile (corresponding to an infinitely long latency)
    qobs(:,2)=qobs(:,2)*defect_p; %normalise to the defective probability
    %if the distribution does not add to 1, create an additional bin for
    %the remaining alternative responses (e.g. errors)
    if(defect_p<1) 
        qobs(end+1,1:2)=[NaN 1];
    end
    %Third column for the probability mass in each bin
    qobs(:,3)=qobs(:,2);
    qobs(2:end,3)=diff(qobs(:,2));
    %Fourth column for the observed frequencies in each bin
    qobs(:,4)=qobs(:,3)*N(1); 
    
    if(defect_p<1)
        nmodel=histc(model,[0;qobs(1:end-1,1)]); %count the number of simulated latencies in between the quantile boundaries
        pmod=nmodel(1:end-1)/N(2); %convert count into probabilities, ignoring the last bin as it will have a count of 0 (corresponding to inf)
        pmod(end+1)=1-(length(model)/N(2)); %probability of the alternative response type (e.g. error responses)
    else
        nmodel=histc(model,[0;qobs(:,1)]); %count the number of simulated latencies in between the quantile boundaries
        pmod=nmodel(1:end-1)/N(2); %convert counts into probabilities, ignoring the last bin as it will have a count of 0 (corresponding to inf)
    end
    pmod(pmod==0)=eps; %avoid numerical problems with predicted probabilities of 0
    
    %We now have the observed and predicted probabilities for each bin. We
    %are ready to compute the relevant goodness-of-fit measures.
    if strcmpi(metric,'Deviance') %twice the negative log-likelihood; my preferred measure as it can easily be used to compute AIC or BIC
        y=-2*sum(qobs(:,4).*log(pmod));
    elseif strcmpi(metric,'Chisquare')
        y=sum((N(1)*(qobs(:,3)-pmod).^2)./pmod);
    else %alternatively, compute G-square statistic
        y=2*sum(qobs(:,4).*log(qobs(:,3)./pmod));
    end        
    
else %use the PDF to obtain the likelihood of the observed data
    pmod=kernel_estimator(model,N(2),max(observed)+1); %returns probability density for t=1...[max(observed)+1]
    flobs=floor(observed); %lower integer bounds
    frobs=observed-flobs; %%fractions remaining after subtracting the lower bounds
    clobs=flobs+1; %upper integer bounds
    flp=pmod(flobs); %probability densities associated with the lower bounds
    clp=pmod(clobs); %probability densities associated with the upper bounds
    pobs=flp+(frobs.*(clp-flp)); %linear interpolation over the 1 ms interval between the lower and upper bounds
    pobs(pobs==0)=eps; %avoid numerical problems with densities of 0
    y=-2*sum(log(pobs)); %-2 times the sum of the log likelihoods
end

end

function [f]=kernel_estimator(simdat,N,maxT)
%Gaussian kernel estimator, as set out in Van Zandt, 2000, PB&R

t=1:maxT;
z=sort(simdat);
q=iqr(z);
s=std(z);
h=(.51/N^.2)*min([q/1.349 s]); %standard deviation of smoothing Gaussian kernel
f=zeros(length(t),1);
for i=1:length(t)
    w=(1/sqrt(2*pi))*exp(-0.5*((z-t(i))/h).^2);
    f(i)=sum(w)/(N*h);
end

end