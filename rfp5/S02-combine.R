x1 = msnset[1:10, ]
x2 = msnset[11:30, ]
x3 = msnset[31:55, ]

l = list(x1, x2, x3)

res = do.call(combine, l)

## combine(combine(combine(x1), x2), x3)

set.seed(11)
m <- matrix(sample(0:5, 12, replace = TRUE), nrow = 4)
m[4, 1] <- 2
m
colnames(m) <- LETTERS[1:3]
rownames(m) <- letters[1:4]
a <- readMSnSet2(as.data.frame(m), eco = 1:3)
fData(a)$notZero <- rowSums(exprs(a) != 0)
fData(a)$rws <- rowSums(exprs(a))
fData(a)$fac <- fData(a)$rws * fData(a)$notZero

set.seed(1)
m <- matrix(sample(0:5, 12, replace = TRUE), nrow = 4)
colnames(m) <- LETTERS[10:12]
rownames(m) <- letters[2:5]
b <- readMSnSet2(as.data.frame(m), eco = 1:3)
fData(b)$notZero <- rowSums(exprs(b) != 0)
fData(b)$rws <- rowSums(exprs(b))
fData(b)$fac <- fData(b)$rws * fData(b)$notZero

exprs(a) <- exprs(a) * fData(a)$fac
exprs(b) <- exprs(b) * fData(b)$fac

a <- updateFvarLabels(a, "A")
b <- updateFvarLabels(b, "B")

x <- combine(a, b)

fData(x)$acc <- c("P1", "P1", "P2", "P2", "P3")

x <- impute(x, method = "zero")

## combineFeatures(x, fcol = "acc", method = "sum")
y <- combineFeatures(x, groupBy = fData(x)$acc, 
                     fun = "sum")

exprs(combineFeatures(x, fcol = "acc", fun = "median"))
exprs(combineFeatures(x, fcol = "acc", fun = "mean"))

?iPQF
