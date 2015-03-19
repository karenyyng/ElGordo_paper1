
Title: The return of the merging galaxy subclusters of El Gordo?
Manuscript ID: MN-14-3494-MJ   
Authors: Ng, Karen (contact); Dawson, William; Wittman, David; Jee, Myungkook; Hughes, John;
Menanateau, Felipe; Sifon, Cristobal   

This	paper is	investigating	the	merger	properties	of	the	“El	Gordo” merging	cluster	using	
a	wide	range	of	data	set.	Although	the	analysis	sounds	convincing,	it	would	benefit	of	a	
better	presentation so	that	the	reader	have	a	better	understanding	on	the	importance	of	
various	constraints.	Once	this	is	corrected,	the	paper	can	certainly	be	accepted.

# Introduction:

> more	references on cluster	merger	and	in	particular	regarding	the typical	speed of	a	merger	should	be	 given.

Accepted. Added in section 

> Figure	1	is	key	for	the	understanding	of	the	paper.	I	would	suggest	to	make	it	
larger	(by	rotating	it	by	~45-50	deg	in	order	to	have	the	elongation axis	of	the	
cluster	horizontal)	

Accepted. After tiling the elongation axis, the RA and DEC coordinates are going to be diagonal.
To remove any potential confusion about the coordinate grid, we have removed the axis ticks.

> showing	possibly	the	galaxy	luminosity	contours.	

Maybe. There may be an overplotting issue.

> Can	the	relic	radio	data	be plotted	on	top	instead	of	some	schematic	of	it? (it	looks	like	the	size	of	the	relic	does	not	match	the	size given	in	Lindner	et	al	2014).

It depends on which wavelength of the radio relic image you are looking at. We
referred to Figure 5 and Figure 8. 



> The	concept	of	time-since-pericenter	is	interesting,	but	it	need	to	be 	clarified (the	
time	at	pericenter	is	only	defined	in	section	3.3).	I	would	suggest	that	you	draw	a	
diagram	of	the	merger	(may	be	at	different	time	step)	so	that	the	reader	can	have	
a	clearer	idea	of	the	geometry	and	evolution	of	the	system (on	the	plot	all	the	
quantities	used	such	as	distances	and	velocities	must	be	indicated).	

Accepted 

> There	you	can	also	define	what	 are	the	different	merger	scenarios that	will	thereafter	 discussed	(the	outgoing	and	the	returning	scenario).	At	the	moment,	the	reader	needs	to	go	through	the	literature	to	understand	what	 are	the	effective	geometric assumptions	used.

Maybe 

# Section	2
> The	presentation	of	the	data	is	not	fully	quantitative. A	table	is	given	for	the	WL	
data,	but	the	paper	is	also	using	the	radio	data	relic,	and	a	summary	table	would
be	great.	This	table	should	underline	the	key	number	from	the	radio	data	used	in	
the	analysis (position,	Mach	Number/velocity?,	polarisation?	...)

The radio data are not inputs like the other data in Table 1. The radio data helps determine the cuts to the inputs.

> What	about	the	velocity	data?	It	is	described	in	section	3.1.1	but	should	it	not	be	
moved	to	section	2?

We have moved some of the descriptions to section 2.

# Section	3

> I	 believe	the	vector	D	is	representing	the	data,	but	this	could	be	clearer,	and	it	
would	be	good	to	clarify	which	data	is	effectively	used.	Table	1	seems	to	only	give	
part	of	the	data	used	in	the	analysis.

We have added descriptions in Table 1 for clarification. 

> It	would	be	helpful	to	describe	a	little	bit	more	the	MC	simulation	code.	Does	the	
simulation	use	a	large	number	of	particules?	Or	is	it	just	using	2	“particles”	with	a	
NFW	mass profile.	It	would	be	good	to	remind	the	reader	of	some	of	the	key	
element	in	D13.

Agreed. 

> Are	the	galaxies	introduces	as	test	particles	in	the	simulations?	It	seems	not,	but	
would	this	be	a	way	to	better	model	and	possibly	constraints the	merger?
I	computed	a	Delta_v_rad	of	463	km/s	(based	on	the	2	redshift:	0.8684	and	
0.8713)	instead	of	the	number	of	476km/s	given	in	the	text,	can	you	explain	
why?
We thank the referee for double-checking our calculation.
The calculations that led to the results were done with the entire PDF distribution, then we take the biweight location as the estimate. In the referee 's calculation, only the best estimate was used, therefore, the discrepancies arise from the different inputs. 
Also, the uncertainty that we give is much bigger than this discrepancy.


> Similarly	the	d_proj	I	found	is	0.744	instead	of	.74	(based	on	the	RA,	DEC	given)

We also used the full PDF for computing $d_{proj}$ which is a more robust calculation.

> It	would	be	much	better	to	define	the	output	parameters	using	a	diagram	of	the	
different	merger	scenarios.

> Would	there	be	a	correlation	between	beta	and	TSP?

No. The parameter $\beta$ is not a property, nor an output of each Monte Carlo simulation, it is a choice of parametrization of our uncertainty, i.e. you cannot find one $\beta$ value for each simulation so you cannot compute a correlation between $\beta$ and $TSP$, since $TSP$ is an output and has one value for each Monte Carlo simulation.  

# Section	4
> It	would	be	good	to	move	some	of	the	key	likelihood	plot	from	the	annex	to	the	
main	text.	Indeed,	the	start	of	section	4,	starts	with	results	that	are	basically	
described	in	the	appendix,	which	is	not	ideal.

Which specific likelihood plot would you like to see in the main text? We feel like there are already many plots in the main text.

# Section	5

> For	the	comparison	with	other	analysis,	I	would	suggest	that	the	author	
summarizes	the	comparison	in	two tables. One	table	comparing	results	with	
other	El	Gordo	modelling.	
A lot of the comparison has been done in the main text.

> And	one	table	comparing	the	properties	of	the	
different	merger	described	with	the	same	modelling	principle,	in	particular	
comparing	El	Gordo	to	the	bullet	cluster	and	the	musket	cluster.

Accepted.

> Written with [StackEdit](https://stackedit.io/).
