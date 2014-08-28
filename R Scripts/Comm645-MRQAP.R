
#=============================================================#
#                                                             #
# Network permutaton tests, correlation & multiple regression #
#  (QAP & MRQAP)                                              #
#  COMM 645 - Communication Networks                          #
#  Katya Ognyanova, 10/10/2012                                #
#                                                             #
#=============================================================#


# This lab will use R packages sna & network.
# Make sure you have those packages installed:

install.packages("sna")
install.packages("network")

# ... and loaded:

library(sna)
library(network) 


# The emon dataset - interorganizational networks
#=============================================================#

# We'll use the emon dataset: interorganizational Search and Rescue
# Networks (Drabek et al.), included in the "network" package.
# The dataset contains 7 networks, each node has 8 attributes

?emon
data(emon)
emon  # a list of 7 networks 

# Cheyenne network (the first one in the list)
ch.net <- emon$Cheyenne 
plot(ch.net)

# Extract the node attributes of the Cheyenne network into a data frame

ch.attr <- data.frame(matrix (0,14,8))
colnames(ch.attr) <-   c("vertex.names", "Command.Rank.Score", "Decision.Rank.Score", 
                       "Formalization", "Location", "Paid.Staff", "Sponsorship", 
                       "Volunteer.Staff")

# Copy each of the 8 vertex attributes to a variable in the ch.attr data frame
for (i in 1:8) { ch.attr[,i] <- (ch.net %v% colnames(ch.attr)[i]) }

ch.attr


# Correlation & Linear Regression in R
#=============================================================#

# Check the correlation between command rank and decision rank scores:
# (remember that ch.attr[[2]] is the same as ch.attr$Command.Rank.Score, etc.)

cor(ch.attr[[2]], ch.attr[[3]])      # Calculate the correlation btw command and decision rank
cor.test(ch.attr[[2]], ch.attr[[3]]) # Examine the significance

# Linear regression: lm(DV ~ IV1 + IV2 + IV3 + ...) 
# DV = dependent variable, IV = independent variables
# As is often the case in R, you can store the regression results in an object.
# Let's call it ch.lm.1 (results from Cheyenne linear regression model #1)

ch.lm.1 <- lm(ch.attr[[2]] ~ ch.attr[[3]] + ch.attr[[6]] + ch.attr[[8]])
summary(ch.lm.1) 

# Calculate centrality measures and store them in the ch.attr data frame:

ch.attr$IndegCent  <- degree(ch.net, gmode="digraph", cmode="indegree")  # indegree centralities
ch.attr$OutdegCent <- degree(ch.net, gmode="digraph", cmode="outdegree") # outdegree centralities 
ch.attr$BetweenCent <- betweenness(ch.net, gmode="digraph") # betweenness centralities

# We can use the centralities in a linear regression model:

ch.lm.2 <- lm(ch.attr[[2]] ~ ch.attr[[6]] + ch.attr[[11]]) # Model with betweenness centrality
summary(ch.lm.2)



# Conditional Uniform Graph (CUG) tests for Network Measures
#=============================================================#

# Let's examine a network measure computed for our observed data.
# (e.g. reciprocity, transitivity, in-degree centralization, etc.)
# How can we tell whether the score we're looking at is relatively high/low?
# Well, we can compute the same measure for all graphs with certain proprerties 
# (for instance all graphs that have the same size/density as our observed network).
# Then we can see how the the measure computed for our observed network compares.

# Of course, in many cases it's not very practical to look at all possible graphs
# of a particular size and density. Instead, we randomly draw graphs from a uniform
# distribution (i.e. each graph has an equal probability of being selected) and compare
# the metrics computed for those graphs with the one we got from our network.
# This is the basic idea behind the conditional uniform graph (CUG) test.

# The function we'll use is cug.test(OurNetwork, Some.network.measure.function, cmode="condition" )


# Conditioning on size
#====================================#
# We're drawing from all possible graphs with the same number of nodes as ours.
# How does the density of our network compare?
# (remember gden is the graph density function)

ch.cug.den <- cug.test(ch.net, gden, cmode="size")
ch.cug.den

 
# Conditioning on number of edges
#====================================#

# Is there more reciprocity in our observed network than we can expect based on its density?

ch.cug.recip <- cug.test(ch.net, grecip, cmode="edges")  
ch.cug.recip
plot(ch.cug.recip)

# Looks like the answer is yes - very few of those 1000 graphs with the same density as
# our network had reciprocity quite as high as we do. 

# In order to do this for centralization, we also have to tell cug.test what the arguments for 
# the "centralization" function are. Remember it takes argument that indicate what type of 
# centralization should be calculated - in-degree, out-degree, betweenness, etc.
# We supply those arguments using FUN.arg= list of arguments for the function:

ch.cug.cent <- cug.test(ch.net, centralization, cmode="edges", FUN.arg=list(FUN=degree, cmode="outdegree"))
ch.cug.cent
plot(ch.cug.cent)


# Conditioning on dyad census
#====================================#

# Here we're comparing our network to graphs that have the same number of dyads that are:
# (1) null - i.e. dyads (i,j)  with no links from i to j or j to i.
# (2) asymmetric - dyads (i,j) with only one directed link ( i -> j or j <- i)
# (3) mutual - dyads (i,j) with two directed links (i <-> j)
#
# Note that this fixes size, density, and reciprocity - all graphs will have the same number
# of reciprocated edges as our observed network.
  
# Given the density and reciprocity of our network, how extreme is its level of transitivity?
 
ch.cug.tr <- cug.test(ch.net, gtrans, cmode="dyad") # Conditioning on dyad census
ch.cug.tr
plot(ch.cug.tr)

# The answer is: very ;)


# Network correlations and the Quadratic Assignment Procedure (QAP)
#=============================================================#

# Permutation tests
# We can restrict the graphs we compare our network to by looking at ones that are
# somewhat similar in structure to our observed data. We generate those by taking the 
# matrix for our network and reshuffling (permutating) its rows and columns.

# This is the idea behind the Quadratic Assignment Procedure (QAP) used to test
# the significance of correlations between networks. 

# We'll use the Padgett florentine marriage & business ties dataset from the ergm package:

install.packages("ergm")
library(ergm)
data(florentine)
detach("package:ergm")

# The data contains two network objects - one with marital and another one
# with business relations between Florentine families.

flobusiness; flomarriage 


# Compute the network correlation:

gcor(flomarriage,flobusiness)
 
# Test the correlation using qaptest:
# qaptest(A.List.Of.Networks, Some.Network.Function, reps=100)
# where reps is the number of graphs to be generated.

# Out of a 100 permutations, how many had >= or <= correlation than the observed value?
# (the parameters g1=1 g2=2 tell the test to compare the first & second element of
# the list of networks that is the first argument of qaptest)

flo.qap <- qaptest(list(flomarriage,flobusiness), gcor, g1=1, g2=2, reps=100)
flo.qap
plot(flo.qap)

# Instead of correlation, we could test another measure of similarity/distance 
# between two networks.


# Multiple regression for network variables: MRQAP
#=============================================================#


# We can similarly use permutaton tests for a linear regression
# on network variables - known as MRQAP (multiple regression QAP)

# Note that the matrices used in the regression can contain different relations, the same
# relation at different points in time, or even dyadic covariate (attribute) data. For instance, 
# you can have a matrix of social ties as a DV and a same-gender attrubute matrix as an IV.
# (i.e. a matrix where (i,j)=1 if i and j are both female(male), and 0 if they have different genders)

# Let's create an array containing 3 random networks, each with 10 nodes:

x <- rgraph(10, 3) 
x[1,,] # matrix 1
x[2,,] # matrix 2
x[3,,] # matrix 3

# A matrix constructed as a linear combination of the first two networks in x:

y <- 3*x[1,,] + 5*x[2,,] + 7

# Now let's run a network regression of y on x and see what happens. 
# netlm takes a single matrix or network object as its first argument (DV) 
# and a stack of matrices (array containing all the IVs) or a list of 
# network objects as its second argument.

net.mod.1 <- netlm(y, x, reps=100)
summary(net.mod.1)

# Oh look - the first two parameters are significant and close to 3 and 5, and the intercept is 7
# Just as you'd expect based on the way we constructed y.



#=============================================================#










