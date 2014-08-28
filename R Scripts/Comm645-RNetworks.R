
#=============================================================#
#                                                             #
#  Networks in R                                              #
#  COMM 645 - Communication Networks                          #
#  Katya Ognyanova, 10/03/2012                                #
#                                                             #
#=============================================================#



# This lab will introduce you to R packages sna & network.
# Make sure you have those packages installed:

install.packages("sna")
install.packages("network") 


#=============================================================#
#                     Network Objects                         # 
#=============================================================# 

# Load the packages.

library(sna)
library(network)
 


# Create an empty network: network.initialize 
#=============================================================#

# 10 nodes, directed adjacency network, no self-loops (0 diagonal)
# multiple edges between the same pair of nodes are not allowed,
# and the network is not a hypergraph

net.1 <- network.initialize(10, directed=TRUE, loops=FALSE, bipartite=FALSE
                            hyper=FALSE,  multiple=FALSE)
net.1

# By default, networks are directed, with no loops, and not bipartite.
# They also don't have multiple edges and are not hypergrpahs.
# So we could get the same network this way:

net.1 <- network.initialize(10)
net.1


# Create a network from a sociomatrix
#=============================================================#

# rgraph creates m random graphs, each of which has n nodes
# drawing from a Bernoulli graph distribution, tie probability = tprob

m <- rgraph(n=10, m=1, tprob=.25) # create 1 random graph, 10 nodes, tie prob =.25
m  # an adjacency matrix

# Turn the matrix m into a network object
# each of the following rows will produce the same result:

net.2 <- network (m, directed=TRUE, loops=FALSE, bipartite=FALSE, matrix.type="adjacency")
net.2 <- network(m)              # construct a network object
net.2 <- as.network(m)           # coerce an object to network
net.2 <- as.network.matrix(m)    # coerce a matrix object to network


summary(net.2) # examine your newly constructed network
plot(net.2)    # look at a simple visualization


# Create a network from an edge list
#=============================================================#

# rgraph can also return a matrix containing an edge list for the graph.

el.1 <- rgraph(n=10, m=1, tprob=.25, return.as.edgelist=T)
el.1 # a matrix containing an edge list 

# The following two lines create identical network objects from the edge list:

net.3 <- network(el.1, directed=TRUE, loops=FALSE, bipartite=FALSE, matrix.type="edgelist")
net.3 <- network(el.1)

summary(net.3)
plot(net.3)


# What are the properties of your network?
#=============================================================#

is.network(net.3)       # Is it a network?
is.bipartite(net.3)     # Does it contain adjacency or affiliation data?
is.directed(net.3)      # Is it directed?
has.loops(net.3)        # Are self-loops allowed?  
network.size(net.3)     # What is the network size?
network.edgecount(net.3) # What is the number of edges?

# Access and change network objects & their elements
#=============================================================#


# You can index network objects much like you do with matrices:

summary(net.3)

net.3[3,5]      # the edge from node 3 to node 5
net.3[1:3,1:3]  # a submatrix of the adjacency matrix
net.3[1,2] <- 1 # changing the edge value
net.3[1,]  <- 1 # change the entire first row

summary(net.3)

# Note that [1,1] remained 0. This is because the network has loops=F.
# So the diagonal will stay 0, even if you try to change it.
# Similarly, networks with directed=F will be symmetrized
# and networks with multiple=F will take only 0 or 1 for edges.

# You can also add and delete edges and vertices using
# add.edges, delete.edges, add.vertices, delete.vertices
 


# Network, node, and edge attributes
#=============================================================#

# The network format has some built-in attributes 
# (e.g. "directed" and "bipartite" for networks, 
# vertex.names for nodes, etc.)
# But you can also add your own.

list.network.attributes(net.3) 
list.vertex.attributes(net.3)
list.edge.attributes(net.3)

# The built-in "na" attribute tells you if a vertex/edge is missing or not.
# Note that this refers to missing data, not to present/absent edges.

# Some edge attribute, set to 77 for all edges:
# (You could of course use a vector with different
#  values for each edge instead)

set.edge.attribute(net.3, "some-edge-attr", 77)
list.edge.attributes(net.3)
get.edge.value(net.3, "some-edge-attr")

# Add a "gender" attribute for each node

node.gender <- c(2,1,2,2,1,2,1,1,2,2)
set.vertex.attribute(net.3, "gender", node.gender)
list.vertex.attributes(net.3)
get.vertex.attribute(net.3, "gender")

# Change the node names via vertex.names attr:

my.node.names <- c("John", "Mary", "Tom", "Bob", "Ann", "Pete", "Kate", "Jill", "George", "Tim")
set.vertex.attribute(net.3, "vertex.names", my.node.names)
get.vertex.attribute(net.3, "vertex.names")
summary(net.3)

# Add a "year" attribute to the network

set.network.attribute(net.3, "year", 2008)
list.network.attributes(net.3)
get.network.attribute(net.3, "year")

# Change the network to symmetric by changing "directed" attr:

set.network.attribute(net.3, "directed", F)
summary(net.3)
 
# We can also get and set network attributes using operators
# %n% for networks, %v% for vertices (nodes), %e% for edges

net.3 %n% "year" <- 2009
net.3 %n% "year"

net.3 %e% "some-edge-attr" <- 88 # This could be a vector assigning different values
net.3 %e% "some-edge-attr"       # to each edge. Here we just assign 88 to all of them.

net.3 %v% "gender" <- 1 # Set gender to female for all nodes
net.3 %v% "gender"
 
# You can also delete attributes:

delete.network.attribute(net.3,"year")
delete.edge.attribute(net.3,"some-edge-attr")
delete.vertex.attribute(net.3,"gender")

summary(net.3)


#=============================================================#
#               Network & Node Measures                       # 
#=============================================================#


# Note that a lot of those will work on a single network, on a list
# of networks, a matrix, or an array with multiple matrices.

# Let's read the Analytic Technologies CAMPNET network we've come to know and love:

campnet <- read.csv(file="campnet.csv", header=T, row.names=1, as.is=T)

campnet.attr <- read.csv(file="campnet-attr.csv", header=T, as.is=T)

campnet; campnet.attr; # Display the files we just imported.

# Create a network object.
# Remember that read.csv gives you a data frame, but we can coerce it into a matrix.

camp.net <- network(as.matrix(campnet))
summary(camp.net)
plot(camp.net)


set.vertex.attribute(camp.net, "vertex.names", campnet.attr$Name) # Add name as a node attribute
set.vertex.attribute(camp.net, "Gender", campnet.attr$Gender)     # Add gender as a node attribute
set.vertex.attribute(camp.net, "Role", campnet.attr$Role)         # Add role as a node attribute

list.vertex.attributes(camp.net)

summary(camp.net)

# Check out the network as an edge list and as a matrix

as.edgelist.sna(camp.net)
as.sociomatrix.sna(camp.net)


# Density
#=============================================================#
# In this and other functions, mode is "digraph" for directed
# and "graph" for symmetric networks.

gden(camp.net, mode="digraph")


# Reciprocity and Mutuality
#=============================================================#

# The number of dyads where i -> j and j <- i.
mutuality(camp.net)

# Reciprocity - proportion of symmetric dyads
# dyadic - ratio of dyads where (i,j)==(j,i) to all dyads
# dyadic.nonnull - ratio of dyads where (i,j)==1 AND (j,i)==1 to all dyads where (i,j)==1
# edgewise - ratio of edges that are reciprocated to all edges where (i,j)==1


grecip(camp.net, measure="dyadic")
grecip(camp.net, measure="dyadic.nonnull")
grecip(camp.net, measure="edgewise")



# Geodesic Distances
#=============================================================#

# By default, nodes that cannot reach each other have a geodesic distance of Inf
# remember, Inf is the constant for infinity.
# Here we'll replace it with the longest theoretically possible path length
# in the network + 1, so 18

camp.geo <- geodist(camp.net, inf.replace=18)
camp.geo

camp.geo$gdist   # The length of the shortest path for all pairs of nodes.
camp.geo$counts  # The number of shortest path for all pairs of nodes.


# Centrality & Centralization
#=============================================================#

# Let's add node centrality measures to the list of attributes in campnet.attr:
# We'll calculate the measures treating the network as symmetric: gmode=graph

campnet.attr$DegreeCent <- degree(camp.net, gmode="graph")       # Degree centrality
campnet.attr$BetweenCent <- betweenness(camp.net, gmode="graph") # Betweenness centrality    
campnet.attr$CloseCent <- closeness(camp.net, gmode="graph")     # closeness centrality

campnet.attr

# Now let's take a look at in-degree, out-degree and total node centralities
# for the directed campnet network

degree(camp.net, gmode="digraph", cmode="indegree")  # indegree centralities
degree(camp.net, gmode="digraph", cmode="outdegree") # outdegree centralities
degree(camp.net, gmode="digraph", cmode="freeman")   # total freeman degree centralities

# Centralization is calculated as a single function with FUN indicating which
# type of measure is to be used (e.g. closeness, degree, betweenness, etc.)

# Degree centralization, network treated as symmetric:
centralization(camp.net, FUN=degree, mode="graph")

# In-degree and out-degree centralization, directed network:
centralization(camp.net, FUN=degree, mode="digraph", cmode="indegree")
centralization(camp.net, FUN=degree, mode="digraph", cmode="outdegree")

# Closeness and betweenness centralization, graph treated as symmetric:
centralization(camp.net, FUN=betweenness, mode="graph")
centralization(camp.net, FUN=closeness, mode="graph")



# Transitivity
#=============================================================#

# Compute transitivity:
# measure = "weak" returns proportion of triads where a->b->c => a->c
# measure ="weakcensus" returns the number of weakly transitive triads.
# measure = "strong" returns proportion of triads where a->b->c <=> a->c
# measure ="strongcensus" returns the number of strongly transitive triads.

gtrans(camp.net, mode="digraph", measure="weak")
gtrans(camp.net, mode="digraph", measure="weakcensus")
gtrans(camp.net, mode="digraph", measure="strong")
gtrans(camp.net, mode="digraph", measure="strongcensus")


# Triad census
#=============================================================#

# Count how many triads of different types exist in your network:

triad.census(camp.net, mode="digraph")

# Triad types (per Davis & Leinhardt):
# 
# 003  A, B, C, empty triad.
# 012  A->B, C 
# 102  A<->B, C  
# 021D A<-B->C 
# 021U A->B<-C 
# 021C A->B->C
# 111D A<->B<-C
# 111U A<->B->C
# 030T A->B<-C, A->C
# 030C A<-B<-C, A->C.
# 201  A<->B<->C.
# 120D A<-B->C, A<->C.
# 120U A->B<-C, A<->C.
# 120C A->B->C, A<->C.
# 210  A->B<->C, A<->C.
# 300  A<->B<->C, A<->C, completely connected.



# Components
#=============================================================#

# How many connected components do we have in the network? components()

# connected ="strong" means v1 and v2 are connected if there is a directed path 
# from v1 to v2 and from v2 to v1.
# connected = "weak" means v1 and v2 are connected if there is a semi-path
# (i.e. path ignoring the link direction) from v1 to v2 and v2 to v1.

# Number of components:
components(camp.net, connected="strong")
components(camp.net, connected="weak")


# Which node belongs to which component? component.dist()

camp.comp <- component.dist(camp.net, connected="strong")
camp.comp

camp.comp$membership # The component each node belongs to
camp.comp$csize      # The size of each component
camp.comp$cdist      # The distribution of component sizes

# Which nodes in the network, if removed, will increase the
# number of components we find?

cutpoints(camp.net, connected="strong")

# Let's remove one of the cutpoints and count components again.

camp.net.cut <- camp.net[-1,-1] # remember "-1" selects all bit the first row/column.
                                # So camp.net.cut will be camp.net with node 1 removed.

components(camp.net.cut, connected="strong")



# Cliques
#=============================================================#

# Let's go ahead and symmetrize camp.net by setting the attribute
# "directed" to F: this will add links (j,i) wherever (i,j) exists.
# You could symmetrize by minimum using symmmetrize() with rule="strong"

set.network.attribute(camp.net, "directed", F) 

# The clique census returns a list with several important elements - let's assign 
# that list to an object we'll call camp.net.cliques.

# The clique.comembership parameter takes values "none" (no co-membership is computed),
#  "sum" (the total number of shared cliques for each pair of nodes is computed),
#  or "bysize" (separate clique co-membership is computed for each clique size)

camp.net.cliques <- clique.census(camp.net, mode = "graph", clique.comembership="sum")

camp.net.cliques # an object that now contains the results of the clique census

# The first element of the result list is clique.count: a matrix containing the  
# number of cliques of different sizes (size = number of nodes in the clique).
# The first column (named Agg) gives you the total  number of cliqies of each size,
# the rest of the columns show the number of cliques each node participates in.
#
# Note that this includes cliques of sizes 1 & 2. We have those when the largest
# fully connected structure includes just 1 or 2 nodes.

camp.net.cliques$clique.count

# The second element is the clique co-membership matrix:

camp.net.cliques$clique.comemb


# The third element of the clique census result is a list of all found cliques:
# (Remember that a list can have another list as its element)

camp.net.cliques$cliques # a full list of cliques, all sizes

camp.net.cliques$cliques[[1]] # cliques size 1
camp.net.cliques$cliques[[2]] # cliques of size 2
camp.net.cliques$cliques[[3]] # cliques of size 3
camp.net.cliques$cliques[[4]] # cliques of size 4



# Structural Equivalence
#=============================================================#

# Structural equivalence - similarity/distance measures include:
# correlation, euclidean distance, Hamming distance, and gamma correlation

sedist(camp.net, mode="digraph", method="correlation" )



#=============================================================#
#                  Affiliation Networks                       # 
#=============================================================#

# Let's take a look at the Davis Southern Women dataset
# a classic affiliation network of women & social clubs

# It's included in the latentnet package, so we can load it from there.

install.packages("latentnet")

library("latentnet")
data("davis")
detach("package:latentnet")

# davis is a network object. Note the attribute "bipartite" equals 18.
# That indicates that we're looking at a 2-mode network: 
# the first 18 nodes of the network are actors and the other 14 are events.

summary(davis)
as.matrix.network(davis)

dav.mat <- as.matrix.network(davis)

# Let's project the affiliation network onto one mode.
# Get a one-mode network of women and their number of shared events
# And a one-mode network of events and their number of shared visitors
# Remember that we get that when we multiply the matrix by its transposed.

dav.women <- dav.mat %*% t(dav.mat)
dav.events <- t(dav.mat) %*% dav.mat

# Dichotomize the valued network of shared events for the Southern women
# method="mean" is the grand mean of the network
# method="absolute", thresh = 0.5 would dichotomize with a cutoff point thresh (here 0.5)
# method="quantile", thresh = .95 would dichotomize by thres quantile (here .95)


dav.women[,] <-  event2dichot(dav.women, method="mean") 
dav.women.net <- network(dav.women) # turn the dichotomized matrix back into a network


# Find node centralities in the dichotomized net. 

degree(dav.women.net, gmode="diraph", cmode="indegree")  # indegree centralities
degree(dav.women.net, gmode="diraph", cmode="outdegree") # outdegree centralities


# Find number of components and cliques.

components(dav.women.net, connected="strong")
dav.women.cs <- clique.census(dav.women.net, mode="digraph")
dav.women.cs
dav.women.cs$clique.count



#=============================================================#
#                 Basic network visualization                 # 
#=============================================================#


library(sna)

# Back to the CAMPNET data


campnet <- read.csv(file="campnet.csv", header=T, row.names=1, as.is=T)
campnet.attr <- read.csv(file="campnet-attr.csv", header=T, as.is=T)
camp.net <- network(as.matrix(campnet))
summary(camp.net)

campnet.attr$DegreeCent <- degree(camp.net, gmode="graph")
campnet.attr

### Visualize the network using the gplot (graph plotting) function from sna.

# The "mode" parameter contains the network layout:
# such as circle, fruchtermanreingold or kamadakawai 
#
# displaylabels = T displays node labels 
# label.cex = 0.5 sets label size to 0.5 of the default size
# label.col = orange sets the color of all labels to orange
# vertex.sides - number of sides for the node symbol (e.g 3=triangle, 4=tetragon, etc.)

gplot(camp.net, mode="fruchtermanreingold", displaylabels=T, label.cex=0.5, label.col="orange", vertex.sides=4)


# Node symbol size can be set using the vertex.cex parameter.
# Let's change node size based on its total degree centrality
# and also hide the arrows for edges:

gplot(camp.net, vertex.cex=campnet.attr$DegreeCent, vertex.sides=18, usearrows=FALSE)

# Network color can be set using the vertex.col parameter.
# This plot will use defaults to color nodes by Gender:

gplot(camp.net, vertex.cex=campnet.attr$DegreeCent, vertex.col=campnet.attr$Gender)

# Let's check out the color names that R knows about:

colors()

# Make a vector that will contain the right color for each node based on gender. 
# We'll use ifelse(condition, value-if-true, value-if-false)

NodeColVec <- ifelse(campnet.attr$Gender=="1","hotpink","lightblue")
NodeColVec

# Let's try to color nodes by gender again, this time using our vector:

gplot(camp.net, vertex.cex=campnet.attr$DegreeCent, vertex.col=NodeColVec)

# If you run the same network visualization more than once, nodes will typically change
# postition.  If you want to visualize the networks multiple times (or show multiple images
# with different relations between the same nodes) it helps to have the nodes keep their
# places across all plots.
#
# Luckily, the network plot returns the coordinates of each node, which we can store in an 
# object and later use again:


my.coordinates <- gplot(camp.net, vertex.cex=campnet.attr$DegreeCent, vertex.col=NodeColVec)


# We have preserved the positions of nodes in an object called my.coordinates.
# Let's create another plot with nodes fixed at the same coordinates using coord=

gplot(camp.net, vertex.cex=5, vertex.col=campnet.attr$Role, coord=my.coordinates)


# Examine the help file to find many more parameters controlling network visualization: 

?sna::gplot 

# Save the visualiziation to a file:

pdf("MyNetworkVis.pdf")  #Start a graphic device that creates PDFs 
gplot(camp.net, vertex.cex=campnet.attr$DegreeCent, vertex.col=NodeColVec) #Plot goes to pdf
dev.off() # shut down the PDF devices



#=============================================================#



