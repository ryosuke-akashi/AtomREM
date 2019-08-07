

Struct_3_94 <- c(-0.35193254E+0001,   -0.35193136E+0001,   -0.35188119E+0001,   -0.35183444E+0001,   -0.35186672E+0001,   -0.49936536E+0001,   -0.49915951E+0001)
Struct_3_94 <- sum(Struct_3_94) / length(Struct_3_94)
Struct_3_8 <- c(-0.26937777E+0001,   -0.43757668E+0001,   -0.36024604E+0001,   -0.36018418E+0001,   -0.43741389E+0001,   -0.43758833E+0001,   -0.36030443E+0001)
Struct_3_8 <- sum(Struct_3_8) / length(Struct_3_8)
Struct_3_722 <- c(-0.27522539E+0001,   -0.42641142E+0001,   -0.27520963E+0001,   -0.27531912E+0001,   -0.42640257E+0001,   -0.42641902E+0001,   -0.50060063E+0001)
Struct_3_722 <- sum(Struct_3_722) / length(Struct_3_722)
Struct_3_707 <- c(-0.27010642E+0001,   -0.35103386E+0001,   -0.27015282E+0001,   -0.42627563E+0001,   -0.35092772E+0001,   -0.42615648E+0001,   -0.50081590E+0001)
Struct_3_707 <- sum(Struct_3_707) / length(Struct_3_707)


d <- read.table(file = "SADDLES_GOOD_STR.dat")

d$V3 <- round(d$V3,4)

d$V4 <- round(d$V4,2)
d$V5 <- round(d$V5,2)
d$V6 <- round(d$V6,2)

sad <- matrix(nrow = length(d$V1),ncol = 8)

sad[,1] <- d$V3

sad[,2] <- d$V4
sad[,3] <- d$V5
sad[,4] <- d$V6
sad[,5] <- d$V7
sad[,6] <- d$V8
sad[,7] <- d$V9
sad[,8] <- d$V10

un_str <-  unique(sad)

un_11 <- un_str[which(un_str[,6] == 1),]
un_another <- un_str[which(un_str[,6] == 0),]

# Delete the duplication by hands! 
un_str <- un_str[-9,]

yr <- c(-4.0,-3.0) # max(un_str[,1])+0.01)
#yr <- c(-3.8,-3.0)
xr <- c(0,length(un_str[,1])+1)


plot(1:length(un_str[,1]),un_str[,1],type = "p",col = "blue", ylim = yr,xlim = xr, 
     xlab="",ylab = expression(paste("Energy, ","10"^ "-12","  erg")),axes = FALSE)

axis(2)

for (i in 1:length(un_str[,1])) {
  lines(c(i,i),c(-3.95,un_str[i,1]-0.016))
  points(i,-3.95*un_str[i,5],pch = 19, cex=1.1, add = TRUE,col = "green")
  points(i,-3.8*un_str[i,6],pch = 18, cex=1.5, add = TRUE,col = "red")
  points(i,-3.722*un_str[i,7],pch = 25, cex=1, bg=1, add = TRUE,col = "blue")
  points(i,-3.7*un_str[i,8],pch = 17, cex=1, add = TRUE,col = "red")
  points(i,un_str[i,1]+0.05,pch = letters[i], cex=1, add = TRUE,col = "black") 
  }

par(xpd=TRUE)
n1 <- "pentagonal bipyramid"
n2 <- "capped octahedron"
n3 <- "tricapped tetrahedron"
n4 <- "bicapped trigonal bipyramid"

legend(xr[1], yr[1]+0.03, c(n1,n2), 
      cex=1.1,col=c("green","red"), pch = c(19,18), pt.cex = c(1.2,1.6), horiz = TRUE, box.lty = 0)


legend(xr[1], yr[1]-0.08, c(n3,n4), 
       cex=1.1,col=c("blue","red"), pch = c(25,17), pt.bg = 1, pt.cex = c(1,1.1), horiz = TRUE, box.lty = 0)

