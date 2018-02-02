## data.frame
mulvey <- read.csv("data/mulvey2015.csv", 
                   stringsAsFactors = FALSE)
?read.csv

dim(mulvey)
ncol(mulvey)
nrow(mulvey)
str(mulvey)
names(mulvey)
colnames(mulvey)
rownames(mulvey)

head(mulvey) ## 6 first obs
tail(mulvey) ## 6 last obs

## subsetting
mulvey[ , 1] ## first col/var
class(mulvey[1,  ]) ## first row/obs
mulvey[1, 1]

## emulate (1) head and (2) tail using [, ]
mulvey[1:6, ]

mulvey[, 2332:2337]

from <- nrow(mulvey) - 5
to <- nrow(mulvey)
mulvey[from:to, ]
mulvey[(nrow(mulvey) - 5):nrow(mulvey), ]

head(names(mulvey))

## subsetting by names
head(mulvey[, 1])
head(mulvey[, "FeatureNames"])
head(mulvey$FeatureNames)


head(rownames(mulvey))

rownames(mulvey) <- mulvey$FeatureNames 

## Question: what data do we have for 
##           protein P17809
## Tip - see line 36

mulvey[4, ]
mulvey["P17809", ]

sel <- mulvey$FeatureNames == "P17809"
sum(sel)
table(sel)
mulvey[sel, ]

which(sel)

names(mulvey)

res <- mulvey[mulvey$pv < 0.01, ]


names(mulvey)

grep("rep", names(mulvey))
