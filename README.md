
Aimming for another review in early April among group members and Davis collaborators for finalizing the methodology and results sections.

Effort to debug Will 's MCMC code 
====
[ ] fit the above catalogs with two different methods (to be
done by week of 04/01)      		
	[ ] least square fitting using bfgs	 
	[ ] MCMC			 
[ ] fix shear to be calculated more exactly (See Issue
wadawson/ResearchCode#10)			
[ ] use James' mock catalog with two halos, fix the mass of one, just
infer mass of the other halo		 

Computation Tasks sorted by priority  
=====
[x] Test sensitivity of polarization prior (over the weekend)   	
[x] rerun MCMC with suitable acceptance rate / variance (early next week)  	 
[ ] rerun MCMAC with updated data (late next week)		 
[ ] Debug MCMC code		 
[ ] modify MCMAC to output the location of the centroids - (which coordinate
 do we want them in???)			
[ ] Grab Will 's code / modify it for plotting all the simulation outputs
(late next week)	 
[ ]  make picture with locations of subclusters, center of mass and
 pericenters    	
[ ] make use of SVM for classification of subcluster members (optional)   	
[ ] finalize all the figures and numerical results before working on the
the writing tasks in detail        	
[ ] zoom in on the merging scenario figure and see which scenario is
  actually more favored for beta = 0.9 with the polarization prior applied
  
  	
	

Writing Tasks sorted by priority  
=====
[ ] incorporate comments from Tues group meetings to paper  
[x] Modify polarization section (over the weekend)   
[ ]  Write result sections  
[ ]  add stuff to appendix this should be quick   
[ ] add Wright and Brainerd to reference   

Table of Content
====
1. Introduction 
2. Method - Monte Carlo Simulation 
	* 2.1 Inputs of the Monte Carlo simulation 
		* 2.1.1 Membership selection and redshift estimation of subclusters 
		* 2.1.2 Weak lensing mass estimation 
		* 2.1.3 Estimation of projected separation 
	* 2.2 Outputs of the Monte Carlo Simulation 		
	* 2.3 Design and applications of priors 
		* 2.3.1 Monte Carlo filters based on the integrated polarization of the radio relic 
	* 2.4 Extension to the Monte Carlo simulation - Determining merger
	scenario with radio relic position 
3. Results 
	* 3.1 Relative merger speed  
	* 3.2 Time-since-collision
	* 3.3 Effects of applied prior on output variables 
	* 3.4 Three-dimensional configuration of El Gordo 
4. Discussion 
	* 4.1 Our findings in the context of other studies of El Gordo
	* 4.2 Comparison to other merging clusters of galaxies 
	* 4.3 Possible improvements for radio relic priors  
	* 4.4 Other limitations and future work 
5. Appendices 
	* App A. Details and output diagnostics for the MCMC mass inference 
	* App B. Plots of the outputs of the Monte Carlo simulation 
