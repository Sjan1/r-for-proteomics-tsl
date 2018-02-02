## variable
x <- 1

## functions
x <- sqrt(4)

?sqrt
?mean

round(x = 3.1415, digits = 2)

round(digits = 2, x = 3.1415)
round(2, 3.1415)

round(3.1415, digits = 2)

## vector
# numeric
x <- c(1, 2, 3, 4, 5, 10, 11)
y <- 1

length(x)
length(y)

## characters
pep <- c("AGVHTW", "AJGYEYY", "GAFDJSL")
length(pep)

class(x)
class(pep)

x2 <- c("1", 2, 3)
class(x2)
x2

x3 <- c(1, "2")
class(x3)

as.numeric("2")

## logicals

l1 <- c(TRUE, FALSE, FALSE)
class(l1)
l1

## 1. create a vector of numeric 
## from 1 to 15. See ?seq

?seq
x <- seq(1, 15)

## 2. create a vector of logicals 
## of length 15, where 10 first elements
## are TRUE and the last 5 are FALSE

length(x)
class(x)
x
l <- x < 11
class(l)
length(l)

## subsetting
x <- rnorm(15)
x
## by indices
x[5]  ## 5th
x[1]  ##first
x[15] ## 15th

x[length(x)] ## last

x[1:3]
x[seq(1, 3)]
x[c(10, 3, 6, 1)]

## by logicals
pep
pep[1:2]
pep[c(TRUE, TRUE, FALSE)]
pep[c(1, 2)]

pep[TRUE]


pep[c(TRUE, FALSE)]

## missing

x <- c(2, 2, 1, NA, 1, NA)
x
sum(x)
sum(x, na.rm = TRUE)
