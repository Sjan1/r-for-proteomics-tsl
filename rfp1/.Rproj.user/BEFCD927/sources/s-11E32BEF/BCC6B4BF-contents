## test id file
library("msdata")

biocLite("msdata")

idf <- ident(full.names = TRUE)

iddf <- readMzIdData(idf)
class(iddf)
nrow(iddf)
colnames(iddf)

table(iddf$isDecoy)
table(iddf$chargeState)

table( iddf$chargeState[ !iddf$isDecoy ] )

table(iddf$chargeState, 
      iddf$isDecoy)

sel <- !iddf$isDecoy

iddf2 <- iddf[sel, ]

length(unique(iddf2$DatabaseAccess))
table(iddf2$DatabaseAccess)
table(table(iddf2$DatabaseAccess))

hist(iddf2$MS.GF.QValue)
boxplot(iddf2$MS.GF.QValue)
mean(iddf2$MS.GF.QValue)

boxplot(iddf$MS.GF.QValue ~ iddf$isDecoy)

mean(iddf[sel, "MS.GF.QValue"])
mean(iddf[!sel, "MS.GF.QValue"])

fiddf <- filterIdentificationDataFrame(iddf, 
                                       verbose = TRUE)