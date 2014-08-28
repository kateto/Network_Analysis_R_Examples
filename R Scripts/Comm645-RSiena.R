
#=============================================================#
#                                                             #
#  Basics of Actor-Based Modeling with RSiena                 #
#  COMM 645 - Communication Networks                          #
#  Katya Ognyanova, 10/24/2012                                #
#                                                             #
#=============================================================#


# Check out http://www.stats.ox.ac.uk/~snijders/siena/
# for more detailed RSiena scripts, tutorials, and data sets.
 
# Install RSiena:
install.packages("RSiena")

# Load the package
library(RSiena) 

# For the following examples, we'll use data from the  Teenage Friends 
# and Lifestyle Study. We have 3 networks representing the friendships
# between 50 girls at 3 points in time:

net.t1 <- as.matrix(read.table("s50-network1.dat"))
net.t2 <- as.matrix(read.table("s50-network2.dat"))
net.t3 <- as.matrix(read.table("s50-network3.dat"))


# We also have some behavioral data measured at each of the 3 points in time:
# Smoking: 1 (non), 2 (occasional) and 3 (regular, i.e. more than once per week).
# Alcohol use: 1 (non), 2 (once/twice a year), 3 (once a month), 4 (once a week) 
# and 5 (more than once a week). 

drink <- as.matrix(read.table("s50-alcohol.dat"))
smoke <- as.matrix(read.table("s50-smoke.dat")) 


# VARIABLES used in Siena models: network, behavior, covariates  
#=================================================================#

# Before you can estimate your model, you have to transform 
# your original variables into objects that can be used by Siena.


# Network variables (DV)
#====================================#

# Put network data in an array with dimensions N x N x T where:
# N is the number of nodes in the network (here 50) 
# and T is the number of time points (here 3)

friend.t123 <- array(c(net.t1, net.t2, net.t3), dim=c(50, 50, 3))

# Create a Siena network object with sienaNet()
# The type could be "oneMode" for one mode networks, 
# "bipartite" for 2-mode networks, or "behavior", 
# for behavioral dependent variables (see next section)

friend.net <- sienaNet(friend.t123, type="oneMode") 

# Take a look at the SienaNet object we created:

class(friend.net)
dim(friend.net)
attributes(friend.net)


# Behavior variables (DV)
#====================================#

# We can also create behavioral variables for a model using sienaNet.
# The difference between covariates (see next section) and behavior is that
# behavior is going to be used as dependent variable, to be explained by the model.
# Other covariates are used as independent variables, to predict network structure.
# So "behavior" here should not be taken literally - those could be attitudes, etc.
  
drink.beh <- sienaNet(drink, type = "behavior")
smoke.beh <- sienaNet(smoke, type = "behavior")


# Constant & changing covariates (IV)
#====================================#

# Covariates are node or edge attributes used as IV (explanatory variables) in the model.

# Constant actor covariates: attributes of the nodes that do not change over time.
# Here all of our attributes change over time, but for the sake of the example, we'll
# create constant covariates from the variables at time 1. We do this with coCovar:

smoke.cc <- coCovar(smoke[,1])
drink.cc <- coCovar(drink[,1])

# Changing (varying) actor covariates: as the name suggests, those are 
# actor attributes that do change over time. We define those with varCovar:

smoke.vc <- varCovar(smoke)  
drink.vc <- varCovar(drink)
 
class(drink.cc)
class(drink.vc)  

# If you want to use another network (or edge attribute) as a covariate,
# you can do that using coDyadCovar for constant dyadic covariates, or
# varDyadCovar for changing dyadic covariates.



# DATA SPECIFICATION: combine dependent & independent variables
#=================================================================#

# Create a siena data object that contains your network & covariates.
# In it we can include our sienaNet objects as well as the constant/changing
# covariates and behavior variables we want to use in a model.

dat.1 <- sienaDataCreate(friend.net, drink.vc, smoke.cc)

# Here's another set of variables that you could use in a model.
# This one includes behavioral dependent variables too:

dat.2 <- sienaDataCreate(friend.net, drink.beh, smoke.vc)

class(dat.2) 

# Just remember not to include covariates in a model if you've also added 
# them as behavioral variables. For example, this would be a bad idea:
# sienaDataCreate(friend.net, drink.beh, drink.vc)


# MODEL SPECIFICATION: Effects included in the Siena model
#=================================================================#

# The next ting you have to specify are the effects to be included in your model.
# Like the data, the effects are stored in an object - a sienaEffects data frame.

# First we create the effects object using getEffects(some.siena.data). 
# You need to have the siena data object as a parameter, because based on the
# number and type of variables in it, you will be able to include different
# effects that are relevant in your case.  

eff.1 <- getEffects(dat.1)

# What you get is a sienaEffects data frame - each row is an effect, and the
# columns are different properties for that effect.
 
eff.1

# eff.1 includes all possible effects that can be used on your data, dat.1
# Effects for which "include" is set to TRUE will be included in the model.
# Some simple effects (rates of change, density, reciprocity, etc.) will
# automatically be included in the model. All others will be FALSE and 
# you'll have to  switch them to TRUE if you want them in your model.

# Let's look at the descriptions of the effects you can use with dat.1:

eff.1$effectName

# Each effect also has a "shortname" used in different functions:

eff.1$shortName

# Which of the effects are included (TRUE) or not included (FALSE)?

eff.1$include
 
# There are a number of different ways to include, remove, or change the effects
# specifying your model. Here we'll use the includeEffects function.

# Let's include an effect for transitivity.
# Its name is transTrip (transitive triplets) - you can read a detailed description 
# of different effects and their meaning in the RSiena manual, section 12.

eff.1 <- includeEffects(eff.1, transTrip)

eff.1

# How do we remove the effect that we've added to the model?
# We can do that adding an include=FALSE parameter to includeEffects.
# By default, if we skip the include attribute, it's assumed to be TRUE.

eff.1 <- includeEffects(eff.1, transTrip, include=FALSE)

eff.1


# Structural effects
#======================================#

# A few structural effects you could include (just a few examples,
# many more are available - look them up in the Siena manual):

# recip     - reciprocity
# cycle3    - 3-cycles (i->j, j->h, h->i)
# transTrip - transitive triplets (i->j, j->h, i->h)
# nbrDist2  - number of actors at distance 2 (i->j, j->h)
# inPop     - in-degree related popularity effect
#             nodes with high in-degree will get more incoming links.
#             (think preferential attachment, etc.)
# outPop    - out-degree related popularity effect
#             nodes with high out-degree will get more incoming links.
# inAct     - in-degree related activity effect
#             nodes with high in-degree will send more outgoing links.
# outAct    - out-degree related activity effect
#             nodes with high out-degree will send more outgoing links.
#...and so on - check out the RSiena manual, section 12.


# Covariate-related effects
#======================================#

# What covariates did we have in our data?
dat.1

# For effects that are premised on node attributes, you have to include
# not only the name of the effect, but also the name of the covariate that
# it applies to. You indicate that with "interaction1="

# Among the covariate effects, we have:
# egoX - actors with high scores will have more outgoing links
# altX - actors with high scores will have more incoming links
# sameX - actors with the same levels on the covariate will be more
#         likely to be connected (typically for categorical vars)
# simX - actors with similar scores on the covariate will be more
#         likely to be connected (typically for continuous vars)
#... and others (see manual)

# Do people who smoke more tend to form more friendship ties?
eff.1 <- includeEffects(eff.1, egoX, interaction1 = "smoke.cc")
eff.1

# Do people who smoke more tend to be more popular?
eff.1 <- includeEffects(eff.1, altX, interaction1 = "smoke.cc")
eff.1

# Are people more likely to form ties with others who have similar smoking level?
eff.1 <- includeEffects(eff.1, simX, interaction1 = "smoke.cc")
eff.1

# Are people more likely to form ties with others who have the same smoking level?
eff.1 <- includeEffects(eff.1, sameX, interaction1 = "smoke.cc")
eff.1

# We can also add effects for more complex interactions between two covariates
# using the function "includeInteraction".  Instead of a single covariate 
# in "interaction1", we will now have a vector of covariates to be included.
# When one of the interacted effects is a structural effect, instead of a 
# covariate name we will include just "" in the interaction1 vector.

# Here we're looking for an interaction effect between smoking and reciprocity:
# (do smokers have a greater tendency to reciprocate friendship ties than non-smokers?)
eff.1 <- includeInteraction(eff.1, egoX, recip, interaction1 = c("smoke.cc",""))
eff.1

# How about an interaction between smoking and drinking? Do people who do both form more ties?
eff.1 <- includeInteraction(eff.1, egoX, egoX, interaction1 = c("smoke.cc","drink.vc"))
eff.1


# Behavior-related effects
#======================================#

# Above we looked at effects for siena data dat.1 that had no behavioral variables.
# However, the dat.2 object we created contained a behavioral/dependent variable:

dat.2

# Let's get the effects for the variables in dat.2:

eff.2 <- getEffects(dat.2)
eff.2

# You'll notice that in the initial model, already included are effects for the rate of 
# change of smoking, as well as linear & quadratic shape effects.
# The last two try to capture trends in the behavior over time that are not related to
# network properties. For instance, as students grow older, they may drink more - which
# could look like a linear increase in our drinking variable.

# Just like we did  for other covariates, we can examine sender, receiver,
# and homophily effects of the behavioral var on the network structure:

eff.2 <- includeEffects(eff.2, egoX, altX, simX, interaction1="drink.beh")

# Behavioral variables, however, can also be *explained* by network structure.
# So that we can test influence mechanisms - do people become more similar to 
# their friends over time?
#
# When we are using a behavioral var as a DV, we add a parameter name="beh.variable"
# to the includeEffects function. In interaction1, as before we include the explanatory 
# variable. Here that would be our network variable friend.net (we're using network 
# characteristics to explain behavior)
#
# Some possible effects to add here include:
# avSim  - average similiarity in behavior between the actor & friends
# totSim - total (summed) similarity in behavior between the actor & friends
# avAlt  - average alter behavior scores (alters are the nodes an actor is tied to)
# indeg  - the effect of indegree on the behavior (e.g. does being popular make you drink?)
# outdeg - the effect of outdegree on behavior    (e.g. does being friendly make you smoke?)
#... and so on - check out the RSiena manual.

eff.2 <- includeEffects(eff.2, name="drink.beh", avAlt,indeg,outdeg, interaction1 = "friend.net")
eff.2 


# MODEL PARAMETER ESTIMATION in RSiena
#=============================================================#

# Now that we know what effects we can use, let's run a model.

# First we'll specify effects as we did above, using the data in dat.2:

eff.3 <- getEffects(dat.2)

# Let's include some structural effects:

eff.3 <- includeEffects(eff.3, transTrip, cycle3)

# A homophily effect for the smoking constant covariate:

eff.3 <- includeEffects(eff.3, simX, interaction1="smoke.vc")

# And an influence effect for the drinking behavior.
# We'll assume a student's drinking behavior is influenced 
# by the average drinking level of their friends:


eff.3 <- includeEffects(eff.3, name="drink.beh", avAlt, interaction1 = "friend.net")

# Take a look at our effects:

eff.3 


# Create a siena project. Your results will be saved in a file with the name you 
# gave as "projname" and an extension .out - so here Student_Behavior_Model.out
# If you run multiple estimations for the same model, the resutls from them will
# be appended at the end of that file. 

mod.1 <- sienaModelCreate(projname='Student_Behavior_Model')

# Finally, let's do the estimation and assign the results to result.1
# To do that we use function siena07() which takes as parameters your
# project, your data object, and your effects object.
# If you don't want the Siena pop-up to show up while you're running
# a model, you can add a parameter batch=TRUE.

result.1 <- siena07( mod.1, data=dat.2, effects=eff.3)
result.1
summary(result.1)
 
# As in previous labs, parameters here are considered significant if they are 
# larger than ~twice (1.96 times) their standard error.
#
# Note that the t-ratio here is not a t statistic for the parameter. It's a convergence
# t-ratio that refers to the difference between simulated and observed values.
# As a rule of thumb, model convergence is considered good if all convergence t-ratios
# are less than .1 in absolute value.

# In this model, we have significant parameters for all rates, density, reciprocity, 
# and transitivity. Our hypotheses about network influence on drinking and homophily in 
# smoking were not supported. The model convergence is good, with all t-ratios < 0.1

# If you do not get good model convergence, you could re-estimate the model using
# the parameters you currently have as the initial values for the new estimation.
# This process can be repeated multiple times until you get good convergence.
# You tell Siena to use the previous estimates using prevAns=results.from.prev.model

result.2 <- siena07( mod.1, data=dat.2, effects=eff.3, prevAns=results.1)
result.2 
 


# GOODNESS OF FIT for RSiena Models
#=============================================================#

# To test goodness of fit for a model with auxiliary statistics, we'll need another
# package called  RSienaTest. It's not yet available on CRAN, but you can find it in
# the R-Forge repository. You can install it here using:

install.packages("RSienaTest", repos="http://R-Forge.R-project.org")

# For some reason (don't even ask), you may get errors if you 
# load the RSiena package before you have loaded RSienaTest. 
# We'll detach RSiena and then load the two packages in the right order.

detach("package:RSiena")

library(RSienaTest)
library(RSiena)

# In order to test goodness of fit, we need to include returnDeps=TRUE
# in the siena07 function. This tells siena to return the simulated networks
# in the results it produces, so we can use them to test goodness of fit.

result.2 <- siena07( mod.1, data=dat.2, effects=eff.3, returnDeps=TRUE )
result.2  

# Check goodnes of fit using the sienaGOF function
# We'll use sienaGOF(some.siena.estimation.results, some.statistic.to.check )
# The statistics to check could be IndegreeDistribution, OutdegreeDistribution, 
# GeodesicDistribution, or TriadCensus.

# How similar is the indegree distribution of the observed data to 
# those of networks simulated based on our model?

ideg.gof <- sienaGOF(result.2, IndegreeDistribution)
ideg.gof
plot(ideg.gof)

# The p value is high, so the observed and simulated stats are similar
# (no statistically significant difference between them)
# This tells us we have good fit for in-degree distribution.

# Similarly, let's check out-degree and geodesic distances:
# (how similar our observed values to those for simulated nets?)

odeg.gof <- sienaGOF(result.2, OutdegreeDistribution)
odeg.gof
plot(odeg.gof)

geo.gof <- sienaGOF(result.2, GeodesicDistribution)
geo.gof
plot(geo.gof)

# The out-degree distribution also looks good, but not the geodesics. 


#=============================================================#


 