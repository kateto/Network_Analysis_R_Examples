
#=============================================================#
#                                                             #
#  Introduction to R                                          #
#  COMM 645 - Communication Networks                          #
#  Katya Ognyanova, 09/26/2012                                #
#                                                             #
#=============================================================#


# R software used in this course:
# R: http://cran.cnr.berkeley.edu 
# RStudio: http://www.rstudio.org


#===============================================#
# Installing R packages
#===============================================#

# Some of the functions you will need in R come in packages:
# topic-specific libraries that you can download from the R 
# project website. If you know the name of the package you need, 
# you can install it from the command line using install.package()


install.packages("sna")
install.packages("network")
install.packages("igraph")
install.packages("statnet")
 

#===============================================#
# Loading an installed R package
#===============================================#

# Before you can use a package, you need to load it using library(). 
# Packages can be detached using detach(package:).

library("MASS")

detach("package:MASS")


#===============================================#
# Loading built-in datasets
#===============================================#

# Many packages include sample datasets.
# You can see a list of those using data()

data()        # list datasets
data(Cars93)  # load this built-on dataset (from package MASS)
Cars93        # take a look at the data
head(Cars93)  # look at the first 6 rows of the data 
              # - helpful   for large files.

#===============================================#
# Help!
#===============================================#

# Find out more about specific functions: help() or ?.
# If you don't remember the exact name of the function you're looking for,
# use double question mark (??).

# Help on function called as.matrix.network
?as.matrix.network

# If you remember this function is in package "network":
?network::as.matrix.network  

# Use ?? if you're not quite sure about the name of the function
??"network matrix"

# You can get help on built in datasets too
?Cars93  # A description of the dataset.


#===============================================#
# Comments
#===============================================# 

# As you might have noticed, comments in R are marked with a "#"
# In RStudio, you can comment/uncomment blocks of text using Ctrl+Shift+C (Windows)


#===============================================#
#  Small things that will go wrong
#===============================================#


# Look out for the following things when using R:
#
#  - R doesn't care about spacing - you can write "a+b", "a + b" or "a    +    b".
#    However it is CASE SENSITIVE. If you get your capitalization wrong, things won't work.
#
#  - As you work with multiple packages, keep in mind that some of them will contain 
#    functions that have the same name. So while you may be trying to use function blah 
#    from package X, what you may actually be using is the blah from package Y. 
#    Detach Y or call the function as X::blah rather than just blah.


#===============================================#
#  Set your working directory
#===============================================#


# The working directory is the place where R will look for and save files, 
# unless otherwise instructed. Note that the directory path is Unix-style -
# it uses slash (/) rather than backslash (\).
# In RStudio you can change the working directory from Tools -> Set Working Directory

# Check the current working directory
getwd() 

# Set the working directory
setwd("C:/USC/Classes/11-2012-Fall/COMM 645/_LABS_/R Folder")


#===============================================#
# What's going on?
#===============================================#

# The data you load, variables you create, etc. are stored in a commong workspace

ls()      # get a list of the objects you have created in the current R session.
search()  # get a list of all currently attached packages.
 

#===============================================#
# Simple calculations
#===============================================#

2 + 2

2^2+2*2

sqrt(abs(-4))


#===============================================#
# Assigning and comparing values
#===============================================#


# You can assign a value to an object using assign() or "<-".  
# You could also use the equal sign "=". Since there are specific 
# cases were <- and = lead to different outcomes, it's a good idea 
# to stick with "<-" for assignment. You can remove objects with rm().

x <- 3         # Assignment
x              # R evaluates the expression and prints the result

y <- 4         # Assignment
y + 5          # Evaluation, note that y remains 4

z <- x + 17*y  # Assignment
z              # Evaluation

rm(z)          # Remove z, this deletes the object.
z              # Try to evaluate z... and fail.


# You can check equality with "==" and inequality with "!=".  
# Traditional comparison operators like less than (<) 
# and less than or equal (<=), etc. will also work.

2==2  

2!=2

x <= y 


#===============================================#
# Special constants: NA, NULL, Inf, Nan
#===============================================#

# NA
# In R, missing or undefined data has the special value NA
# You can check whether values are missing with is.na()
# If you use NA in an operation (e.g. 5+NA), the result 
# will generally also be NA

# NULL
# NULL represents an empty object - e.g. a null/empty list
# You can check whether values are NULL with is.null()

# Inf anf -Inf represent positive and negative infinity
# which can be returned by certain mathematical operations
# like division of a number by zero:

5/0
is.inf(5/0)

# NaN
# NaN or Not a Number is assigned when the result of an operation
# can not be reasonably defined - for instance 0/0

0/0
is.nan(0/0)


#===============================================#
# Vectors 
#===============================================#


# You can create a vector object by combining elements with c()

v1 <- c(1, 7, 11, 22)       # Numeric vector, length 4
v2 <- c("hello","world")    # Character vector, length 2 (a vector of strings)
v3 <- c(TRUE, TRUE, FALSE)  # Logical vector, same as c(T, T, F)

(T+T)*(T+T+T)               # As you remember from Boolean algebra, T=1, F=0

# Vectors can only contain elements of a single type. If you 
# try to make R combine types, it will coerce your elements
# into the least restrictive type. 

v4 <- c(v1,v2,v3,"boo") 	#All elements magically turn into strings


# Other ways to create vectors:

v <- 1:7         # same as c(1,2,3,4,5,6,7)  

v <- rep(0, 77)  # repeat zero 77 times: v is a vector of 77 zeroes
length(v)        # check the length of the vector

v <- rep(1:3, times=2) # Repeat 1,2,3 twice  

v <- rep(1:10, each=2) # Repeat each element twice  

v <- seq(10,20,2) # sequence: numbers between 10 and 20, in jumps of 2  


#===============================================#
# Vector operations 
#===============================================#

# NOTE that the operations listed in this section can also be applied 
# to the rows or columns of matrices, arrays and data frames.

# Element-wise operations

v1 <- 1:5         # 1,2,3,4,5. This assignment erases the old value of v1.
v2 <- rep(1,5)    # 1,1,1,1,1. 

v1 + v2      # Element-wise addition
v1 + 1       # Add 1 to each element
v1 * 2       # Multiply each element by 2
v1 + c(1,7)  # This doesn't work: (1,7) is a vector of different length

# Standard mathematical operations can be used with numerical vectors:

sum(v1)      # The sum of all elements
mean(v1)     # The average of all elements
sd(v1)       # The standard deviation
cor(v1,v1*5) # Correlation between v1 and v1*5 

# Logical operations

v1 > 2       # Each element is compared to 2, returns logical vector (TRUE or FALSE for each element)
v1==v2       # Are corresponding elements equivalent? Returns logical vector.
v1!=v2       # Are corresponding elements *not* equivalent? Same as !(v1==v2). Returns logical vector.

(v1>2) | (v2>0)  # | is the boolean OR, returns a vector.
(v1>2) & (v2>0)  # & is the boolean AND, returns a vector.

(v1>2) || (v2>0)  # || is the boolean OR - returns a single value
(v1>2) && (v2>0)  # && is the boolean AND - ditto

# Vector elements

v1[3]             # third element of v1
v1[2:4]           # elements 2, 3, 4 of v1
v1[c(1,3)]        # elements 1 and 3 - note that your indexes are a vector
v1[c(T,T,F,F,F)]  # elements 1 and 2 - only the ones that are TRUE
v1[v1>3]          # v1>3 is a logical vector TRUE for elements >3

# Note that all comparisons return logical vectors that can be assigned to an object:
LogVec <- v1 > 3  # Logical vector, same length as v1, TRUE (or T) where v1 elements > 3
v1[LogVec]        # Returns only those elements from v1 that are > 3

# To add more elements to a vector, simply assign them values.
# Elements 5 to 10 are assigned values 5,6,7,8,9,10:
v1[6:10] <- 6:10

# We can also directly assign the vector a length:
length(v1) <- 15 # the last 5 elements are added as missing data: NA



#===============================================#
# Matrices & Arrays
#===============================================#


# A matrix is a vector with dimensions:

m <- rep(1, 20)   # A vector of 20 elements, all 1
dim(m) <- c(5,4)  # Dimensions set to 5 & 4, so m is now a 5x4 matrix

# Create a matrix using matrix():

m <- matrix(data=1, nrow=5, ncol=4)  # same matrix as above, 5x4, full of 1s
m <- matrix(1,5,4) 			         # same matrix as above
dim(m)                               # What are the dimensions of m?

# Create a matrix by combining vectors:

m <- cbind(1:5, 5:1, 5:9)  # Bind 3 vectors as columns, 5x3 matrix
m <- rbind(1:5, 5:1, 5:9)  # Bind 3 vectors as rows, 3x5 matrix

m <- matrix(1:10,10,10)

# Select matrix elements: 

m[2,3]  # Matrix m, row 2, column 3 - a single cell
m[2,]   # The whole second row of m as a vector
m[,2]   # The whole second column of m as a vector
m[1:2,4:6] # submatrix: rows 1 and 2, columns 4, 5 and 6
m[-1,]     # all rows *except* the first one

m[1,]==m[,1]  # Are elements in row 1 equivalent to corresponding elements from column 1? 
m>3           # A logical matrix: TRUE for m elements >3, FALSE otherwise
m[m>3]        # Selects only TRUE elements - that is ones greater than 3

# Transpose m:
t(m)          # Again, note that this will not change the matrix, 
# it only evaluates t(m) and prints the result
m <- t(m)     # Now m is transposed. 

m %*% t(m)    # The operator %*% does matrix multiplication
m * m         # Note that * is element-wise multiplication

# When one or two dimensions (vectors & matrices) are not enough,
# arrays are the way to go.
# Created with the array() function:

a <- array(data=1:18,dim=c(3,3,2)) # 3d with dimensions 3x3x2
a <- array(1:18,c(3,3,2))          # a shorter way to get the same array


#===============================================#
# Factors
#===============================================#

# Factors are used for categorical data.

eye.col.v <- c("brown", "green", "brown", "blue", "blue", "blue") #vector

eye.col.f <- factor(c("brown", "green", "brown", "blue", "blue", "blue")) #factor

# R finds the different levels of the factor - e.g. all distinct values. The data is stored
# internally as integers - each number corresponding to a factor level.

levels(eye.col.f)  # The levels (distinct values) of the factor (categorical variable)

as.numeric(eye.col.f)  # The factor as numeric values: 1 is  blue, 2 is brown, 3 is green
as.numeric(eye.col.v)  # The character vector, however, can not be coerced to numeric


#===============================================#
# Lists
#===============================================#

# Lists are collections of objects (e.g. of strings, vectors, matrices, other lists, etc.)
# In contrast with vectors, matrices, and arrays, lists can hold different type of objects
# that also don't have to have the same length

l1 <- list(boo=v1,foo=v2,moo=v3,zoo="Animals!")  # A list with four components
l2 <- list(v1,v2,v3,"Animals!")

l1["boo"]      # Access boo using square brackets - this returns a list.

l1[["boo"]]    # Access boo using double square brackets - this returns the numeric vector

l1[[1]]        # Returns the first component of the list, equivalent to above.

l1$boo         # Named elements can also be accessed using the $ operator - equivalent to [[]]

# Add more elements to a list:

l1[[5]] <- "More elements!" # The list l1 had 4 elements, we're adding a 5th here.
l1[[8]] <- 1:11 # We added an 8th element, but not 6th or 7th. Those will be created empty (NULL)
l1$Something <- "A thing"  # Adds a ninth element - "A thing", called "Something"


#===============================================#
# Data Frames
#===============================================#

# The data frame is a special kind of list used for datasets.
# Think of rows as cases, columns as variables. Each column is a vector or factor.

head(Cars93)    # Let's take another look at the cars93 dataset
class(Cars93)   # What is the data format? Look, it's a data frame!

# Creating our own dataframe:

dfr1 <- data.frame(ID=1:4,FirstName=c("John","Jim","Jane","Jill"), Female=c(F,F,T,T), Age=c(22,33,44,55))
dfr1$FirstName   # This is the second column of dfr1. Notice that R thinks this is a categorical variable
# and so it's treating it like a factor, not a character vector.

# Let's get rid of the factor by telling R to treat FirstName as vector:
dfr1$FirstName <- as.vector(dfr1$FirstName)

# Alternatively, you can tell R you don't like factors from the start using stringsAsFactors=FALSE:
dfr2 <- data.frame(FirstName=c("John","Jim","Jane","Jill"),stringsAsFactors=FALSE)
dfr2$FirstName   # Success: not a factor.


# Access elements of the data frame
dfr1[1,]   # First row, all columns
dfr1[,1]   # First column, all rows
dfr1$Age   # Age column, all rows
dfr1[1:2,3:4] # Rows 1 and 2, colimns 3 and 4 - the gender and age of John & Jim
dfr1[c(1,3),] # Rows 1 and 3, all columns

dfr1$Age <- c(15,25,35,45) # Change the Age column to this vector

# Let's find the age of Jill.

dfr1[4,4]

# What if we don't know exactly which row is Jill?
# Find all rows for which FirstName == Jill. For those rows, select column 4 (Age):

dfr1[dfr1$FirstName=="Jill",4] 

# Find the names of everyone over the age of 30 in the data

dfr1[dfr1$Age>30,2]

# Find the average age of all females in the data:

mean ( dfr1[dfr1$Female==TRUE,4] )


# Operations on datasets work just like those for matrices.

dfr1$Age * 10            # Age multiplied by 10 (Evaluation, Age does not change)
dfr1$ID <- dfr1$ID * 10  # Replace ID with ID*10 (Assignment, ID does change)
dfr1$ID > 20             # A boolean vector, TRUE for cases with ID > 20

# Add a new column/variable to the data frame
dfr1$Human <- c(TRUE, TRUE, TRUE, TRUE)
dfr1[,6]   <- 1:4

# Add a new case/row to the data frame
dfr1[5,2] <- "Pete" # Note that we just set the name for this case.
# All other values are missing: NA

# We can also turn matrix objects into data frames:
m1 <- matrix(1:100,10,10)
m1
d1 <- as.data.frame(m) 
d1

class(m1); class(d1) 

# m1 is a matrix, it can contain only one kind of elements.
# If we replace one element with a string, all other elements
# of the matrix will turn into strings.

m1[1,1] <- "Boo!" 

# The same is not true of a data frame, which can store 
# different types of objects.

d1[1,1] <- "Boo!"



#===============================================#
# Reading and writing files
#===============================================#

# You can read multiple data formats into R, including text, excel, spss, and SAS files.
# Here we'll look at reading .csv files. Note that the data is read as a data frame
# CSV is a useful format - many programs (e.g. Excel) can save data as csv.

# Read MyFile.csv, header should be TRUE if you have column names on top, FALSE otherwise.
# stringsAsFactors=FALSE if you have strings that are not categorical variables (e.g. names)

MyData <- read.csv(file="MyFile.csv", header=TRUE, stringsAsFactors=FALSE)
class (MyData) # Note that the data is read as a data frame

# Write R data to csv. Cars93 to be saved to MyFile.csv. 
# na = what to use for missing values.  

write.csv(Cars93, file="MyFile.csv", na="999")


# The package "foreign" contains functions that will allow you to read many other formats.
# You can use it, for instance, to read files saved in SPSS or SAS.

library(foreign) # load the library
?read.spss   # look up the function reading SPSS files       
?read.ssd    # look up the function reading SAS files
?read.dta    # look up the function reading Stata files


#===============================================#
# Control Flow: Conditions and Loops
#===============================================#


##### Conditional statements

# if (condition) expr1 else expr2

x <- 5; y <- 10

if (x!=0) y <- y/x else y <- 0  # If x is not 0, assign y/x to y, otherwise y is 0

y
 
##### Loops 

# For loop
# for (variable in sequence) expr

ASum <- 0; AProd <- 1

for (i in 1:x)    # Repeat the following for each i between 1 and x
{
  ASum <- ASum + i
  AProd <- AProd * i
}

ASum; AProd

# While loop
# while (condintion) expr

while (x > 1) {print(x); x <- x-1;} # while x>1 is true, repeat the following.

# Repeat loop
# repeat expr, use break to exit the loop


repeat { print(x); x <- x+1; if (x>10) break} # repeat until x>10, then break
 

#===============================================#
 







