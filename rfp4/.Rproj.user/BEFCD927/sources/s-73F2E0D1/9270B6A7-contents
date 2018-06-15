e <- as(msnid, "MSnSet")
nrow(e)

set.seed(1)

x1 <- sample(featureNames(e), 500)
x2 <- sample(featureNames(e), 500)
x3 <- sample(featureNames(e), 500)
x4 <- sample(featureNames(e), 123)
x5 <- sample(featureNames(e), 400)
x6 <- sample(featureNames(e), 300)

library("Vennerable")

## library("devtools")
## install_github("js229/Vennerable")

v <- Venn(list(x1, x2, x3, x4))
v
plot(v)
v@IntersectionSets[["1111"]]
v@IntersectionSets[["1010"]]

## install.packages("UpSetR")
library("UpSetR")
x <- fromList(list(x1, x2, x3, x4, x5, x6))
names(x) <- LETTERS[1:6]
upset(x, nsets = 6)
?upset
upset(x, sets = c("A", "E", "C", "D"),
      keep.order = TRUE)

