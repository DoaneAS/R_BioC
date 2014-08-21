#Five Ways to subset x
 #1 with positive integers, which return elements at the specified positions.##
###

x[c(3, 1)] #> [1] 3.3 2.1
x[order(x)] #> [1] 2.1 3.3 4.2 5.4

# Duplicated indices yield duplicated values
x[c(1, 1)] #> [1] 2.1 2.1

# Real numbers are silently truncated to integers
x[c(2.1, 2.9)] #> [1] 4.2 4.2


#2 with negative integers, which omit elements at the specified positions
x[-c(3, 1)] #> [1] 4.2 5.4
#It's an error to mix positive and negative integers in a single subset:

x[c(-1, 2)] #> Error: only 0's may be mixed with negative subscripts


#3 with a LOGICAL vector, which selects elements where the corresponding logical value is TRUE. Because you write the expression that creates the logical vector, this is probably the MOST USEFUL  type of subsetting.

x[c(TRUE, TRUE, FALSE, FALSE)]#> [1] 2.1 4.2
x[x > 3]#> [1] 4.2 3.3 5.4

##If the logical vector is shorter than the vector being subsetted, it will be recycled to be the same length.
x[c(TRUE, FALSE)]#> [1] 2.1 3.3
# Equivalent to
x[c(TRUE, FALSE, TRUE, FALSE)]#> [1] 2.1 3.3

#A missing value in the index always yields a missing value in the output:
  x[c(TRUE, TRUE, NA, FALSE)] #> [1] 2.1 4.2  NA

#with NOTHING, which returns the original vector unchanged. This is not useful in 1d but it's very useful in 2d. It can also be useful in conjunction with assignment.
x[]#> [1] 2.1 4.2 3.3 5.4


#with ZERO, which returns a zero-length vector. This is not something you'd usually do on purpose, but it can be helpful for generating test data.

x[0]
#> numeric(0)

#If the vector is named, you can also subset with: a CHARACTER vector, which returns elements with matching names.

(y <- setNames(x, letters[1:4]))#>   a   b   c   d  #> 2.1 4.2 3.3 5.4

y[c("d", "c", "a")]
#>   d   c   a 
#> 5.4 3.3 2.1

# Like integer indices, you can repeat indices
y[c("a", "a", "a")]
#>   a   a   a 
#> 2.1 2.1 2.1

# When subsetting with [ names are always matched exactly
z <- c(abc = 1, def = 2)
z[c("a", "d")]
#> <NA> <NA> 
#>   NA   NA



####MATRICES
a <- matrix(1:9, nrow = 3)
a
colnames(a) <- c("A", "B", "C")
a[1:2, ]
#>      A B C
#> [1,] 1 4 7
#> [2,] 2 5 8
a[c(T, F, T), c("B", "A")]
#>      B A
#> [1,] 4 1
#> [2,] 6 3
a[0, -2]
#>      A C




df <- data.frame(x = 1:3, y = 3:1, z = letters[1:3])
df
#selecting rows
df[df$x == 2, ]

#>   x y z
#> 2 2 2 b
df[c(1, 3), ]
#>   x y z
#> 1 1 3 a
#> 3 3 1 c

# There are two ways to select columns from a data frame
# Like a list:
df[c("x", "z")]
#>   x z
#> 1 1 a
#> 2 2 b
#> 3 3 c
# Like a matrix
df[, c("x", "z")]
#>   x z
#> 1 1 a
#> 2 2 b
#> 3 3 c

# There's an important difference if you select a single column:
# matrix subsetting simplifies by default, list subsetting does not.
str(df["x"])
#> 'data.frame':    3 obs. of  1 variable:
#>  $ x: int  1 2 3
str(df[, "x"])
#>  int [1:3] 1 2 3

#subsetting genes
fData(smallset) <- annot

Gene.symbol <- fData(smallset)$Gene.symbol
ESR1 <- grep("ESR1", Gene.symbol)
#or
ESR1v <-grep("ESR1", fData(smallset)$Gene.symbol)
ESR1 == ESR1v #TRUE
ESR1ex <- smallset[ESR1,]

#SUBSETTING EXPRESSIONSETS#
###########################
PAM50 <- phenoData(VDXset)$PAM50
stime <- phenoData(VDXset)$t.dmfs
Basal <- grep("Basal", PAM50)
LumA <- grep("LumA", PAM50)
BasalEx <- VDXset[,Basal]
LumAEx <- VDXset[,LumA]
LumAs <- LumAEx[,1:10]
Basals <- BasalEx[,1:10]
smallset <- combine(LumAs, Basals)
##########################



boxplot(exprs(smallset))
boxplotplus2(t(exprs(smallset)))

y <- (exprs(smallset))
boxplot(y, col=as.numeric(pData(smallset)$PAM50)+1)
Index <- as.numeric(pData(smallset)$PAM50)
d <- rowMeans(y[,Index==3]) - rowMeans(y[, Index==1])
a <- rowMeans(y)
smoothScatter(a, d, main="MA plot", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")


boxplotplus2(t(exprs(ESR1ex)), pt.pch=16, .las=1)


g <- y[4752,]
m <- mean(g[Index==1])-mean(g[Index==3])
plot(jitter(Index), g, col=Index, axis(1, labels=c("RES", "SENS"), at=1:2))





library(help=survival) # see the version number, list of available functions and data sets.
data(aml)



########Matching and replacing######
#
#if x1 values are unique and row numbers are the same (use %in% first), define a new variable x2 is df1 that is equal to df2x2:
df1$x2[match(df2$x1,df1$x1)] <- df2$x2
#if values are not unique:
for(id in 1:nrow(df2)){
  df1$x2[df1$x1 %in% df2$x1[id]] <- df2$x2[id]
}
