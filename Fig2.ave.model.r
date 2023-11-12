rm(list=ls())

########### 7 gene interaction model

treat.c = 'y'
n.sub = 7
n.int = 7

# generate all possible genotypes
cap.m = c()
for (i in 1:n.sub){
  cap.m = cbind(cap.m, rep(c(1,0), each=2^(n.sub-i), times=2^(i-1)))
}
cap.m.c = c()
for (i in 1:ncol(cap.m)){
  let.c = sapply(cap.m[,i], function(x){
    if (x){
      return(LETTERS[i])
    } else {
      return(letters[i])
    }
  })
  cap.m.c = cbind(cap.m.c, let.c)
}
all.geno = apply(cap.m.c, 1, function(x) paste0(x, collapse=''))
rownames(cap.m) = all.geno
colnames(cap.m) = LETTERS[1:n.sub]

### generate the matrix - this was changed from network reconstitution
ord.int = apply(cap.m, 1, sum)
rec.m = cap.m
rec.mx = c()  ## rec.mx was introduced
if (n.int > 1) {
  for (i in 1:n.int){  ## 2:n.int in net reconst was changed to 1:n.int
    if (i < n.sub){
      comb.s = cap.m[ord.int == i,]
      rec.m1 = apply(comb.s, 1, function(x){
        apply(cap.m, 1, function(y){
          sum(x * y) == i
        })
      })
      rec.m1 = t(apply(rec.m1, 1, function(x) {
        if (sum(x)==0) {
          return(rep(0,length(x)))
        } else {
          return(x/sum(x))
        }
      }))
      colnames(rec.m1) = apply(comb.s, 1, function(x){
        paste0(LETTERS[which(x==1)],collapse=':')
      })
    } else {
      rec.m1 = rep(0, nrow(cap.m))
      rec.m1[ord.int==n.sub] = 1
    }
    rec.mx = cbind(rec.mx, rec.m1)   ## rec.m in net reconst was changed to rec.mx
  }
  colnames(rec.mx)[colnames(rec.mx)=='rec.m1'] = paste0(LETTERS[1:n.sub],collapse=':')  ## rec.m in net reconst was changed to rec.mx
}
if (treat.c == 'y'){
  rec.mx = cbind(rec.mx, remainder=1) ## rec.m in net reconst was changed to rec.mx
}

## rec.mx is the matrix for averaging model
av7.mat = rec.mx
nr7.mat = rec.mx
nr7.mat[,1:7] = apply(nr7.mat[,1:7],2, function(x) {
  x[x > 0] = 1
  return(x)
})
anov7.mat = rec.mx
anov7.mat = apply(anov7.mat,2, function(x) {
  x[x > 0] = 1
  return(x)
})

## anov7 random
set.seed(5)
anov7.rand = t(sapply(1:1e4, function(x) {
  z.v = runif(128, 1, 10)
  solve(anov7.mat, z.v)
}))

## nr7 random
set.seed(5)
nr7.rand = t(sapply(1:1e4, function(x) {
  z.v = runif(128, 1, 10)
  solve(nr7.mat, z.v)
}))

## av7 random
set.seed(5)
av7.rand = t(sapply(1:1e4, function(x) {
  z.v = runif(128, 1, 10)
  solve(av7.mat, z.v)
}))

### plot, Fig 2
eff.numb = c(7,choose(7,2),choose(7,3),choose(7,4),
             choose(7,5),choose(7,6),1)
col1 = rev(rainbow(7, start=0,end=0.8))
col2 = rep(col1, eff.numb)
eff.numb.b = cumsum(eff.numb)
eff.numb.b = c(0.5, eff.numb.b+0.5)
leg = paste(1:7, 'gene', sep='-' )
pdf('outputs/Fig2.random.7genes.anov.nr.av.pdf', height=10, width=6.5)
opar=par(mfrow=c(3,1), mar=c(1,7,1,1))
boxplot(anov7.rand[,-128], # remove the intercept
        xaxt ='n', outline=F,
        col=col2, ylim=c(-100,100)) 
abline(v=eff.numb.b,lty=3, col='gray70')
x.pos = (eff.numb.b[-1]+eff.numb.b[-7])/2
x.pos[7] = 129
x.pos = x.pos-2
text(x.pos, -100, leg, pos=4, srt=90, col=col1, cex=1.5)
text(120,93,'Additive Model', pos=2, cex=1.5)
boxplot(nr7.rand[,-128], xaxt ='n', outline=F,
        ylim=c(-15,15), col=col2)
abline(v=eff.numb.b,lty=3, col='gray70')
text(120,14,'NR Model', pos=2, cex=1.5)
boxplot(av7.rand[,-128], xaxt ='n', outline=F,
        ylim=c(-15,15), col=col2)
abline(v=eff.numb.b,lty=3, col='gray70')
text(120,14,'Averaging Model', pos=2, cex=1.5)
par(opar)
dev.off()


