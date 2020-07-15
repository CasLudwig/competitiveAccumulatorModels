function inp=gen_inputs(N,resp_type,gain,resp_parms,delay)

% Generates a 3-D matrix for N trials (rows), specifying the internal
% response over time (columns), for each of M response alternatives
% (depth). The number of alternatives is specified by the length of the
% gain vector. The gains determine the amplitudes of the internal responses,
% and - in a sense - the discriminability between response alternatives.
% Note that the first element of the gain vector corresponds to the
% "target" and so will typically have a larger value (to ensure it tends to
% win the race to threshold).
% 
% The internal response may simply be a constant (as it is in any time
% homogenous accumulator model), or it may follow some more
% realistic physiological profile (e.g. transient or sustained response).
% In the latter case, the drift rate in the subsequent accumulation stage
% would vary over time (i.e. a time inhomogenous process). We use a gamma
% function or its integral to model the response of transient and sustained
% channels respectively. In this case, the program expects resp_parms to
% contain 2 additional parameters.
% 
% Finally, we allow for a temporal asynchrony between the onset of the
% target and non-target responses (with more than 1 non-target, the
% responses to all non-targets start together). The delay parameter is
% defined with respect to the onset of the target response. Negative values
% then indicate that the non-target(s) were responded to earlier than the
% target. Positive values indicate that the target response precedes the
% non-target responses. Best to make sure that delay (in milliseconds) is
% divisible by the time step: rem(delay,h)=0. Any delay is fixed - there is
% no trial-to-trial variability in the delay. The earliest response always
% starts at t=0.

global h;
global MAXT;

gain_sd=resp_parms(1); %variability in the gain of the visual mechanism: translates into variability in the accumulation rate (across alternatives and trials)
if(resp_type>0) %we only need to know the following parameters if the internal response is not constant
    n_stages=resp_parms(2); %number of stages in leaky, cascaded integration equation (Watson, 1986; Smith, 1995; Ludwig, 2009)
    scale=resp_parms(3); %1/scale is the rate constant of each stage of the transient or sustained channel
end

M=length(gain); %number of response alternatives
t=0:h:MAXT; %time vector for the whole trial duration
t_short=0:h:MAXT-abs(delay); %time vector for the restricted internal response(s) in the presence of a delay
inp=zeros(N,length(t),M); %THE matrix to be filled
gains=normrnd(repmat(gain,N,1),repmat(gain_sd,N,M)); %actual noisy gain values
delay_index=ceil(abs(delay)/h)+1; %index at which the later response will start

switch resp_type
    case 1 %transient internal response: use the functional form of a Gamma pdf
        y=gampdf(t,n_stages,scale);
        y_short=gampdf(t_short,n_stages,scale);
    case 2 %sustained internal responses: use the functional form of a Gamma cdf
        y=gamcdf(t,n_stages,scale);
        y_short=gamcdf(t_short,n_stages,scale);
    otherwise %constant internal responses
        y=ones(1,length(t));
        y_short=ones(1,length(t_short));
end

%now everything is in place to fill the input matrix
for i=1:M
    if((delay<0)&&(i==1)) %non-target(s) first, then target
        inp(:,delay_index:end,i)=repmat(y_short,N,1).*repmat(gains(:,i),1,length(t_short));
    elseif((delay>0)&&(i>1)) %target first, then non-target(s)
        inp(:,delay_index:end,i)=repmat(y_short,N,1).*repmat(gains(:,i),1,length(t_short));
    else %either no delay or dealing with the earliest response (i.e. the one that starts at t=0)
        inp(:,:,i)=repmat(y,N,1).*repmat(gains(:,i),1,length(t));
    end
end
inp=inp*h; %scale with the magnitude of the time step (in anticipation of temporal integration in the subsequent stage)

%illustrate the profile of the internal response - you will probably want
%to comment this out at some point
plot(t,y*gain(1),'k-','linewidth',2); %show the response profile
set(gca,'FontSize',12,'PlotBoxAspectRatioMode','manual','XLim',[0 500]);
xlabel('time (ms)');
ylabel('internal response (arbitrary units)');