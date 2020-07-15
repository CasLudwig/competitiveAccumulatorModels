function [choice,rt]=run_competitive_accumulator

%Simulates one particular instantiation of the competitive accumulator
%model. Outputs: Nx2 choice matrix, with the first column indicating the
%accuracy for each trial (1 = correct, 0 = error) in column 1 and the
%winning accumulator (1...M) in column 2; Nx1 rt vector with the
%corresponding reaction time (including non-decisional delay). Trials in
%which no accumulator reached threshold within the maximum allowed decision
%time are denoted NaN. If there are many of those you probably should
%adjust the parameter values.

global h; %time step in ms; global because it is used in all subroutines
global MAXT; %maximum latency

%seed the random number generators; note that when examining different
%parameters you will probably want to keep the seed constant - otherwise
%variations from one simulation to the next are not only a result of the
%differences in parameter values, but also of the different initialisation
%of the RNG (which determine the various forms of internal noise).
seed=sum(clock*100);
rand('twister',seed);
randn('state',seed);

h=5;
MAXT=1000; %maximum decision latency
N=1000; %number of trials to simulate

%model parameters of interest - not all need to be free parameters of
%course!
%response: 0 - "standard" accumulator model with linear rate of rise (i.e. the internal
%response is a constant); 1 - time inhomogenous accumulator model with a transient internal response;
%2 - time inhomogenous accumulator model with a sustained internal response
response=0;
%gain: used to scale the internal response, and so determines the mean
%accumulation rate; the length of this vector corresponds to the number of
%alternatives, the first element corresponds to the "target" gain
gain=[1 .25]; %e.g. for 3-AFC: [2 1.5 1.5]
gain_sd=0.2; %variability around the gain, resulting in variation in drift across accumulators and trials
n=9; %first parameter of the internal response function (only for transient and sustained responses)
s=10; %second parameter of the internal response function (only for transient and sustained responses)
beta=.01; %weight of lateral inhibition
theta=200; %threshold - this value will strongly depend on the nature of the internal response and the value of the gain(s)
start=[25 0]; %e.g. for 3-AFC: [25 0 0]; as before, the first element corresponds to the target starting point
startvar=[0 0]; %uniform variability around the mean start points: mean +/- startvar (make sure [mean-startvar]>=0)
t0=100; %(upper limit of) non-decisional delay
t0_unif=0; %binary variable to indicate whether the non-decisional delay is variable (if ndvar==1, t0 is the upper limit of a uniform distribution)
delay=0; %difference in onset (ms) in target/non-target responses: negative - non-target(s) before target; positive - non-target(s) after target (see 'gen_inputs')

%first step is to simulate some internal responses
int_responses=gen_inputs(N,response,gain,[gain_sd n s],delay);
%'gen_inputs' will generate a graph of the "average" internal response over
%time. At some point you will probably want to comment this out.

%next, we pass these responses on to the accumulator stage, where they are
%further noise-perturbed and integrated over time
[choice,rt]=accumulator(int_responses,start,startvar,[beta theta t0 t0_unif]);
%'accumulator' will generate a graph of the "average" accumulator paths
%over time (not limited by the response threshold); this is useful to come
%up with a reasonable value for the threshold.