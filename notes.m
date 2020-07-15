% C.J.H. Ludwig, September 2009 - July 2013
% Last update: 17 July 2013 (change of assumptions concerning the
% variability in the non-decisional delay).

% The code for a particular instantiation of an (competitive) accumulator
% model consists of three functions:
% 'run_competitive_accumulator': sets a number of parameters that defines a
% particular instantiation and runs the model.
% 'gen_inputs': generates some internal responses as a function of time
% 'accumulator': integrates the internal responses over time
% 
% An additional function, 'compute_fit', is provided to enable you to
% compare an empirical latency distribution with a simulated one. Check the
% help documentation for that function for an overview of how to use it.
% 
% To see what a function does and how to use it, either just read the
% introductory information provided in the relevant m-file, or type 'help
% [function name]' in the command window. In addition, the code itself is
% extensively commented. 
% 
% The main function (run_competitive_accumulator) only simulates a set of
% rts and choices. It does not compare this to any observed data or,
% indeed, fit it to observed data. This is because certain sub-sets of
% simpler models will have analytic solutions (e.g. Brown and Heathcote's
% LBA model is contained within this code by choosing a constant internal
% response and setting the "within-trial" noise to 0). The more sources of
% noise you include the more likely it is that fitting a model will involve
% simulation of a large number of trials (as in Usher & McClelland, 2001;
% Ratcliff & Smith, 2004; Ludwig, 2009), in which case this code could
% potentially be used for that purpose. For this reason I have
% included the 'compute-fit' function which works for full probability
% densities as well as for quantile distributions.
% 
% Common features of these models that have not (yet) been implemented are:
% - leakage (Usher & McClelland, 2001; Brown & Heathcote, 2005);
% - recurrent self-excitation (Usher & McClelland, 2001);
% - temporal variation in the threshold (Smith, 2000; Ditterich, 2006;
%   Ludwig, 2009);
% 
% With regard to timing, this will strongly depend on a) the speed of your
% computer; b) the number of trials you are simulating; c) the size of the
% time step used for the simulation; d) the number of choice alternatives.
% For instance, on my laptop I get the following approximate timings (all
% for 1000 trials):
% time step = 1 ms; M-AFC = 2; running time ~70 seconds
% time step = 1 ms; M-AFC = 3; running time ~160 seconds
% time step = 5 ms; M-AFC = 2; running time ~ 3 seconds
% time step = 5 ms; M-AFC = 3; running time ~ 7 seconds
% As you can see, both the time step and the number of choice alternatives
% make a huge difference, and things scale non-linearly. In the current
% code I have set the time scale to 5 ms.
% 
% So, to get started just type: '[choice,rt]=run_competitive_accumulator'
% in the command window. 
% To then determine choice accuracy compute sum(choice(:,1))/N. Bear in
% mind that choice may contain some NaNs if the parameters are not very
% well-chosen (as a result of trials in which there is no winning
% accumulator within the maximum time allowed). You probably want to filter
% these out. To plot the distribution of correct RTs you just need to
% select the rts of interest (e.g. correct, error: choice(:,1) = 1 or 0
% respectively), again after filtering out any NaNs (if present).