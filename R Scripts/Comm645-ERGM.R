
#=============================================================#
#                                                             #
#  Basics of Exponential Random Graph Models in R             #
#  (Adapted from Statnet.org papers)                          #
#  COMM 645 - Communication Networks                          #
#  Katya Ognyanova, 10/10/2012                                #
#                                                             #
#=============================================================#


# This lab will use the R packages ergm & sna:


install.packages("ergm")
install.packages("sna") 

library(ergm)
library(sna) 


# ERGM on a symmetric network: the Florentine families
#=============================================================#

# As is traditional in ergm tutorials, we'll use the Padgett florentine marriage & business ties
# dataset included with the ergm package:
 
data(florentine) 

# The data contains two network objects - one with marital and another one
# with business relationships between Florentine families.

flobusiness; flomarriage; 

# Exponential random graph models - what terms can we use in a model?

help('ergm-terms')

# Let's estimate a simple  model which only examines density (edge term)
# The format of the ergm command is ergm(YourNetwork ~ Signature1 + Signature2 + ...) 
# where YourNetwork can be a matrix or a network object.

flo.mar.1 <- ergm(flomarriage ~ edges)  		 
flo.mar.1
summary(flo.mar.1)

# We get a negative edge parameter since the network is rather sparse.
# The edge parameter here is the log of the edge odds, i.e. log(#dyads-w-edge/#dyads-no-edge)
# The network has 20 ties of 120 possible ties. Let's calculate the log odds ourselves:
# [ remember that an event with probability p has odds of p/(1-p) and log odds of log(p/(1-p)) ]

log(20/(120-20)) # We get -1.609, the same as the edge parameter in the erg model.

# The corresponding probability is .167:
exp(-1.609)/(1+exp(-1.609)) # you can also get that using inv.logit() from package "boot"
 
# Next we look at a fancier model that includes triangles in addition to edges:

flo.mar.2 <- ergm(flomarriage ~ edges + triangles, seed=1)    	 
flo.mar.2
summary(flo.mar.2)

# The triangle coefficient is not significant - so this is not a signature 
# driving the network formation. What do the coefficients tell us?
# Conditional log-odds of a tie between two actors here =
# = -1.675*(change in the number of ties) + 0.158 * (change in the number of triangles)
# = -1.675*1 + 0.158*(change in the number of triangles)
# 
# if the tie will not add any triangles to the network, its log-odds = -1.675.
# if it will add one triangle to the network, its log-odds = -1.675 + 0.158 = -1.517
# if it will add two triangles to the network, its log-odds = -1.675 + 0.158*2 = -1.359
# The corresponding probabilities are 0.158, 0.180, and 0.204.
#
# (note: we're using a stochastic algorithm - so you will get slightly different estimates
# if you rerun the model. Here we use seed=1 to make sure we'll get the same results every time)

# There are a large number of other structural signatures you could add paramteters for.
# For instance 2-stars: kstar(2), 3-stars: kstar(3) isolates: isolates, etc. 

# Let's run a model checking whether edges in the Florentine business network are predicted by
# edges in the marriage network. To do that, we can use an edge covariate parameter edgecov()
# As in: ergm(MyNetwork ~ Signature1 + Signature2 + ... + edgecov(AnotherNetwork))

flo.mar.3 <- ergm(flobusiness ~ edges + edgecov(flomarriage))       
flo.mar.3
summary(flo.mar.3)

# We can also use node attributes in an erg model.
# For the Florentine families, we have an attribute called "wealth" in the network object.

w.vec <- flomarriage %v% 'wealth'  # Store the node wealth in a numeric vector.
w.vec
gplot(flomarriage, vertex.cex=w.vec/20)	# plot the network with vertex size proportional to wealth

# Let's test whether the edge probabilities are a function of wealth:
# Are wealthy families more likely to form ties?  

flo.mar.4 <- ergm(flomarriage ~ edges + nodecov("wealth"))       
flo.mar.4
summary(flo.mar.4)

# Yes, there is a significant positive main effect for wealth:
# - The p-value for the wealth parameter makes it significant at the .05 level.
# - It's positive, which means we see more of that configuratoin than we'd expect by chance.


# ERGM on a directed network: Sampson Monastery
#=============================================================#

# ERG model of a directed network - the liking relations between monks
# in Sampson's dataset.

data(samplk)
samplk1; samplk2; samplk3
plot(samplk3)

# Is there a statistically significant tendency for ties to be reciprocated?

samp.mod.1 <- ergm(samplk3 ~ edges + mutual)	
summary(samp.mod.1)				

# Conditional log-odds of two actors forming a tie =
# = -2.15 * change in the number of ties + 2.3 * change in number of mutual dyads
# If adding the tie will not make a dyad reciprocal, its log-odds = -2.15
# if it will add a mutual dyad to the network, its log-odds = -2.15 + 2.3 = 0.15  

 


# ERGM with node attributes: Faux Mesa High
#=============================================================#

# Faux mesa high is simulated data representing a high-school friendship network.
# Attributes for each node (student) include gender, race, and grade.


data(faux.mesa.high)  			
fmh.net <- faux.mesa.high
plot(fmh.net)						
fmh.net


# Taking a look at gender 
plot(fmh.net, vertex.col='Sex')

# Taking a look at the grade of the students
plot(fmh.net, vertex.col='Grade') 

# Taking a look at the race of the students
plot(fmh.net, vertex.col='Race') 


# A simple model that includes just the edge (density) parameter:
fmh.mod.1 <- ergm(fmh.net ~ edges)
summary(fmh.mod.1)


# NODEMATCH 
# Are nodes with the same attribute levels more likely to be connected?
# Do high-school students tend to have friends of the same grade?

fmh.mod.2 <- ergm(fmh.net ~ edges + nodematch("Grade"))
summary(fmh.mod.2)

# We can add an attribute diff=T to nodematch to get a separate parameter for
# each level of the categorical variable. 
# Here, a separate parameter for each grade:

fmh.mod.3 <- ergm(fmh.net ~ edges + nodematch("Grade", diff=T))
summary(fmh.mod.3)


# How about gender and race?  

fmh.mod.4 <- ergm(fmh.net ~ edges + nodematch("Grade") + nodematch("Race") + nodematch("Sex"))
summary(fmh.mod.4)


# NODEMIX
# Nodemix will add a parameter for each combination of levels for the categorical variable.
# Let's look at the parameters for edges between students from different race groups:

fmh.mod.5 <- ergm(fmh.net ~ edges + nodemix("Race"))
summary(fmh.mod.5)

table(fmh.net %v% "Race")  			# Check out race frequencies
mixingmatrix(fmh.net, "Race")   # Check out # of links between/within groups

# Note that we got -Inf parameters in the model for configurations 
# that don't exist in the observed network at all.


# NODEFACTOR
# Main effect of a categorical attribute.
# Are some types of nodes more likely to form ties than others?
# For example, are boys forming friendship ties more actively than girls?

fmh.mod.6 <- ergm(fmh.net ~ edges + nodematch("Grade", diff = T) + nodefactor("Sex"))
summary(fmh.mod.6)

# Negative parameter for males means females are more actively forming friendships.

# NODECOV
# Main effect of a continuous attribute (we'll treat grade as continuous here).
# Are nodes with high levels on a continuous attribute more likely to form ties?
# Let's check if students with higher values on attribute "Grade" tend to form more friendships.

fmh.mod.7 <- ergm(fmh.net ~ edges + nodecov("Grade") + nodematch("Sex"))
summary(fmh.mod.7)

# Note that this is the parameter version for undirected networks.
# For directed networks, we have nodeicov (for incoming links)
# and nodeocov (for outgoing links).
# Similarly nodefactor has directev versions nodeifactor & nodeofactor.

 
# ABSDIFF
# For continuous attributes: are people more likely to be connected to others
# who have similar values on an attribute? Absdiff = abs(ValueNode1-ValueNode2)
# Here, are students more likely to have friends close to their own grade?
# (that is, links i->j are more likely for smaller values of abs(grade of i - grade of j))

fmh.mod.8 <- ergm(fmh.net ~ edges + absdiff("Grade") + nodematch("Sex"))
summary(fmh.mod.8)




# Simulating networks based on a model
#=============================================================#
 
# After we have estimated model coefficients, we can draw graphs from
# the probability distribution defined by those parameter values.
# If our model was good, the graphs we draw from this distribution
# should be similar to our observed data.

# Simulate 15 networks based on the fmh.mod.6 model:

fmh.mod.8.sim <- simulate(fmh.mod.8, nsim=15)
summary(fmh.mod.8.sim)

# All the simulated network are stored in the returned object:
class(fmh.mod.8.sim)

# We can access any of them and take a look at it:
fmh.mod.8.sim[[1]]




# Goodnes of Fit and MCMC diagnostics
#=============================================================#

# After estimating parameters for your mode, you want to know how well
# it fits the observed data.

# Let's check the goodness of fit for one of our initial models of the Padgett network.
# Check how well the degree distribution of the networs generated from our model
# match the degree distribution of the observed network:

summary(flo.mar.4) # Take a look at the model

flo.mar.4.gof <- gof(flo.mar.4 ~ degree) # goodness of fit for degree distribution

# If this was a directed network, we could check gof for in- or out-degree instead
# using gof(flo.mar.4 ~ idegree) or gof(flo.mar.4 ~ odegree)

flo.mar.4.gof # Take a look at the observed & simulated values
plot(flo.mar.4.gof) # plot the observed & simulated values

# The resutls contain 1 row for each possible node degree (e.g. row 0 - number of isolates,
# row 1 - number of nodes with only 1 link, row 2 - # of nodes with 2 links, etc.)
# The first column contains the counts from the observed network, the other give staticstics
# from the simulated networks.  

# The fit is not bad - observed values are within the confidence interval (see the plot).
# P values are high (the observed & simulated values do not differ significantly).
# This is one of those rare cases where a high p value is a good thing :)

# We can check the goodness of fit with regard to other network statistics.
# For instance geodesic distance.
# Compare our network with 20 simulated networks based on the flo.mar.4 model:

flo.mar.4.gof2 <- gof(flo.mar.4 ~ distance, nsim=20) # gof based on 20 simulated nets
summary(flo.mar.4.gof2)
plot(flo.mar.4.gof2)

# Here each row in the summary is the number of geodesics with a particular length.
# For instance, we have 20 node pairs in the observed network with shortest paths of 1 
# (those correspond to the 20 edges in our observed network).


# Model diagnostics (for MCMC)
# Information about the model that can help diagnose problems.
# Note we can't get (and don't need) these diagnostics for flo.mar.4 since
# it was not estimated using MCMC. This is because it was simple enough 
# (i.e. a dyadic independence model) that we did not need MCMC estimation.

mcmc.diagnostics(flo.mar.2)

# It's easier to go through the charts if we save in a PDF file:

pdf("flo_mar_model2.pdf")
mcmc.diagnostics(flo.mar.2)
dev.off()

# We can examine the diagnostics to see whether our model looks ok,
# check for model degeneracy, see if MCMC sample size and burn-in are large enough, etc. 
# Read more about interpreting model diagnostics in Section 5 of this document:
# http://statnet.csde.washington.edu/trac/raw-attachment/wiki/Resources/ERGMtutorial.pdf



#=============================================================#
 
