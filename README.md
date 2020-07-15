**Disclaimer: the last update of this toolbox was made in 2013. Since then plenty of other code has been developed (in open source languages such as R) and been made openly available for simulating and fitting sequential sampling and accumulator models. _I_ don't even use this code anymore, but I have given it to students and post-docs occasionally, typically just to give them a starting point for developing their own code. So I'm documenting it here, just so I don't lose it.**

Sequential sampling and accumulator models form a powerful framework to model choice and latency, and can be successfully applied to account for the 'when' and 'where' of saccadic eye movement behaviour (e.g. Carpenter's LATER model and some of our own work). The basic idea is that "evidence" in favour of different potential saccade targets is accumulated over time towards a response threshold. These models come in many flavours, but they all share the central concept of temporal integration of some noisy decision variable to criterion. Models differ in what types of noise they assume (noise in drift rates, starting point variability, etc.) and in whether they include components like lateral inhibition between accumulators, leakage, and self-excitation. Quite clearly then, the space of all possible models in this domain is really quite large.

I have written some basic (Matlab) code here that will allow you to simulate quite a wide range of models within this general class. The code, in its current form, simulates a competitive accumulator that includes:

- Lateral inhibition;
- Noise in the mean accumulation rate (or drift rate), from trial-to-trial, and from one accumulator to the next;
- Noise within the accumulation process (often referred to as 'within-trial noise');
- Starting point variation between accumulators;
- Uniform variability around the mean starting point(s);
- Three different 'input profiles':
  - Constant input: this is the standard model in which the mean accumulation rate is constant within an accumulation process;
  - Transient input: a physiologically more plausible input that peaks rapidly and then decays to baseline;
  - Sustained input: a physiologically more plausible input that rises to some elevated level and then maintains this activation over time;
  - Temporal asynchronies between the onset of accumulators;
- As many choice alternatives as your machine can cope with.

In its current form the model(s) does not include:

- leakage (Usher & McClelland, 2001; Brown & Heathcote, 2005; Ludwig et al., 2007);
- recurrent self-excitation (Usher & McClelland, 2001);
- temporal variation in the threshold (Smith, 2000; Ditterich, 2006; Ludwig, 2009).

In many model fits, leakage is cancelled out by self-excitation, so it seemed to make sense to omit these components. Although some evidence for temporal variation in the threshold exists (see Ditterich, 2006; Ludwig, 2009), behavioural data are often well accounted for without such additional complexity. 

By changing some of the parameters, some well-established models are nested within this architecture. For instance, using a constant input and setting the within-trial noise and lateral inhibition both to 0, the model corresponds to the Linear Ballistic Accumulator (Brown & Heathcote, 2008). Likewise, setting the noise in the mean accumulation rate to 0 gives you a version of the Leaky Competitive Accumulator (Usher & McClelland, 2001), with leakage cancelled out by recurrent self-excitation (Bogacz et al., 2006). This model, in turn, is capable of mimicking the classic diffusion model (e.g. Ratcliff & Smith, 2004).

This code is really just intended for those initial simulations that often precede a formal model fitting exercise, further experimentation or the formulation of more qualitative theories (e.g. just as a check to see whether this kind of model could, in principle, account for your findings). Although the code could be used to actually fit observed data, it may be that one of the simpler, nested models is sufficient for this purpose. In that case analytic or at least more efficient solutions may be available that will work much faster.

However, with many different sources of noise such models can often become unwieldy and simulation is the only feasible way to go (Brown, Ratcliff, & Smith, 2006). For this reason, the "toolbox" also contains a function that allows you to compare a simulated with an empirically observed (defective) latency distribution. The goodness-of-fit measures include the deviance (-2 times the log-likelihood), Chi-square, and G-square. The deviance may be computed for both binned and unbinned distributions. That said, with model predictions derived from simulation I recommend using bins or quantiles. My own preferred metric for goodness-of-fit is the deviance, as it is easily transformed into AIC or BIC, which enables straightforward comparison between models with different numbers of free parameters. Note that any set of parameters that minimises the deviance will also minimise G-square, but not Chi-square. For all measures lower values indicate a better fit.

To get started, download the matlab files to a folder on the Matlab path and begin by reading 'notes.m' (or type 'help notes' in the command window). To check how to compare an observed with a simulated distribution type 'help compute_fit'.

Finally:

- There may well be lots of errors in this code, so do let me know if you come across any.
- You will need the statistics toolbox for the code to run, though the only function that relies on this toolbox (I think) is 'normrnd'. If you haven't got the stats toolbox, you can use the standard Matlab 'randn' function to draw random numbers from a Gaussian distribution.
- I have not done any compatibility checking. As far as I am aware this code should work in Matlab 2007b onwards (but quite possibly before that as well). The 2009 releases appear to prefer a different way of initialising the random number generators, but are still backward compatible with the 'twister' and 'state' calls used here.
- This work is ongoing and so stuff will be added and changed as we go along depending on our needs and whims!
