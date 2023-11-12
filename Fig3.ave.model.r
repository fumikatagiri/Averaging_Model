############ Tsuda et al. data analysis
rm(list=ls())
### package
library(nlme)
library(RColorBrewer)
library(stringr)

### load data
input.d.n = 'AvrRpt2_ETI.txt'
input.dat = read.delim(paste0('Supplemental.Dataset1/', input.d.n), header=T, colClasses='character')
cat(paste('Data file:', input.d.n, '\n'))

summary(input.dat)

### fix the input.dat 1
input.dat = input.dat[input.dat$genotype != 'rpm1/rps2' & input.dat$genotype != 'npr1',]
colnames(input.dat)[colnames(input.dat) == 'colony'] = 'value'
input.dat$treatment[input.dat$treatment == 'avrRpt2'] = 't'
input.dat$treatment[input.dat$treatment == '_pLAFR'] = 'c'

geno.names = c()
for (geno.n in input.dat$genotype) {
  if (grepl('dde2',geno.n)) gn = 'a' else gn = 'A'
  if (grepl('ein2',geno.n)) gn = c(gn, 'b') else gn = c(gn, 'B')
  if (grepl('pad4',geno.n)) gn = c(gn, 'c') else gn = c(gn, 'C')
  if (grepl('sid2',geno.n)) gn = c(gn, 'd') else gn = c(gn, 'D')
  geno.names = c(geno.names, paste0(gn, collapse = ''))
}
input.dat$genotype = geno.names

## check columns for 'genotype' and 'value'
col.n = colnames(input.dat)
if (!('genotype' %in% col.n)){
  cat('No genotype column')
  stop()
}
if (!('value' %in% col.n)){
  cat('No value column')
  stop()
}
input.dat$value = as.numeric(input.dat$value)
factor.col = col.n[col.n!='value']
for (i in factor.col){
  col.val = input.dat[,i]
  col.val = trimws(col.val,which='both')
  input.dat[,i]=factor(col.val)
}

### input info 1. treatment
## check the data
treat.c = 'y'
if (treat.c == 'y'){
  if (!('treatment' %in% col.n)){
    cat('No treatment column')
    stop()
  }
  if (!setequal(levels(input.dat$treatment), c('c','t')) ) {
    cat('treatment levels must be c or t')
    stop()
  }
}

### input info 2. subnetwork numbers
n.sub = 4
cat(paste(n.sub, 'subnetworks\n'))
## check the genotype values
geno.v = input.dat$genotype
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
# separate the cases for treatment
if (treat.c == 'y'){
  dat.1 = input.dat[input.dat$treatment == 'c', 'genotype']
  dat.2 = input.dat[input.dat$treatment == 't', 'genotype']
  dat.1 = unique(as.character(dat.1))
  dat.2 = unique(as.character(dat.2))
  if ( prod(dat.1 %in% all.geno) == 0 | prod(dat.2 %in% all.geno) == 0 ){
    cat('The genotype value is not correct')
    stop()
  }
  if (!setequal(dat.1,all.geno) | !setequal(dat.2,all.geno)){
    cat('Warning: missing genotypes')
  }
} else {
  dat.1 = unique(as.character(input.dat$genotype))
  if (!setequal(dat.1,all.geno) ){
    cat('Warning: missing genotypes')
    stop()
  }
}

########################################################
out.dir = 'outputs.fig3/'
####### n.int = 4
### input info 3: up to what order of interactions
n.int = 4

### generate the matrix
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
av4.mat = rec.mx
nr4.mat = rec.mx
nr4.mat[,1:4] = apply(nr4.mat[,1:4],2, function(x) {
  x[x > 0] = 1
  return(x)
})
anov4.mat = rec.mx
anov4.mat = apply(anov4.mat,2, function(x) {
  x[x > 0] = 1
  return(x)
})

#### NR model first
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'nr.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'nr.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'nr.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'nr.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = nr4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### ANOVA additive model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'add.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'add.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'add.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'add.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = anov4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### averaging model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'Av.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'Av.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'Av.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'Av.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = rec.mx[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()

############################
######### n.int = 3
### input info 3: up to what order of interactions
n.int = 3

## matrices
rec.mx = rec.mx[,-15]
av4.mat = av4.mat[,-15]
nr4.mat = nr4.mat[,-15]
anov4.mat = anov4.mat[,-15]

#### NR model first
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'nr.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'nr.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'nr.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'nr.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = nr4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects

height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### ANOVA additive model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'add.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'add.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'add.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'add.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = anov4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### averaging model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'Av.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'Av.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'Av.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'Av.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = rec.mx[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


###########################################
########## n.int=2
### input info 3: up to what order of interactions
n.int = 2

### generate the matrix
rec.mx = rec.mx[,-(11:14)]
av4.mat = av4.mat[,-(11:14)]
nr4.mat = nr4.mat[,-(11:14)]
anov4.mat = anov4.mat[,-(11:14)]

#### NR model first
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'nr.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'nr.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'nr.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'nr.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = nr4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### ANOVA additive model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'add.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'add.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'add.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'add.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = anov4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### averaging model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'Av.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'Av.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'Av.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'Av.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = rec.mx[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


########################
##### n.int = 1
### input info 3: up to what order of interactions
n.int = 1

rec.mx = rec.mx[,-(5:10)]
av4.mat = av4.mat[,-(5:10)]
nr4.mat = nr4.mat[,-(5:10)]
anov4.mat = anov4.mat[,-(5:10)]

#### NR model first
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'nr.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'nr.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'nr.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'nr.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = nr4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### ANOVA additive model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'add.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'add.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'add.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'add.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = anov4.mat[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)

out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()


#### averaging model
##### mixed-effect linear model
### output file names
outp.n = paste(input.d.n, n.int, 'Av.output.txt', sep='.')
boxp.n = paste(input.d.n, n.int, 'Av.boxplot.pdf', sep='.')
diag.n = paste(input.d.n, n.int, 'Av.diagnostic.pdf', sep='.')
plot.n = paste(input.d.n, n.int, 'Av.reconstitution.pdf', sep='.')

outp.n = paste0(out.dir, outp.n)
boxp.n = paste0(out.dir, boxp.n)
diag.n = paste0(out.dir, diag.n)
plot.n = paste0(out.dir, plot.n)

### generate the matrix for the data
m. = rec.mx[as.character(input.dat$genotype),]
m.[input.dat$treatment=='c',] = 0

left.m = 5
if (n.sub > 5) left.m = left.m + (n.sub-5)/3
### model fitting with boxplot of the data for fixed effects
height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
if (height < 8) height = 8
pdf(boxp.n,width=7,height=height)
opar = par(mar=c(5, left.m, 2, 2))
boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
        xlab='value')
par(opar)
dev.off()

## the model with treatment
lme.rec = lme( value ~ m. -1 + genotype,
               random =~ 1|replicate/flat/pot,
               data = input.dat)
out.t = summary(lme.rec)$tTable
out.t = data.frame(Contribution=rownames(out.t),out.t)
write.table(out.t, file=outp.n, sep='\t', quote=F, row.names=F)

### diagnostic plots
pdf(diag.n, width=5, height=10)
opar=par(mfrow=c(2,1))

## Q-Q plot
qqnorm(resid(lme.rec),cex=0.7,
       main = 'Q-Q: Residual vs. Standard Normal')
qqline(resid(lme.rec), col='red',lty=2)

## Resid vs fitted plot 
plot(fitted(lme.rec),resid(lme.rec), cex=0.7,
     main = 'Residual vs. Fitted',
     xlab='Fitted', ylab='Residual')
abline(h=0, col='red',lty=2)
par(opar)
dev.off()

## barplot the output
out.t1 = out.t[grep('m.',rownames(out.t)),c('Value','Std.Error','p.value')]
r.name = rownames(out.t1)
r.name = sub('m.','',r.name)
mean.v = out.t1[,'Value']
names(mean.v) = r.name
error.m = out.t1[,'Value']-out.t1[,'Std.Error']
error.p = out.t1[,'Value']+out.t1[,'Std.Error']
max.range = range( c(error.m, error.p) )
max.range = max.range * 1.3
if (max.range[2]>0 & max.range[1]> -0.35*max.range[2]){
  max.range[1]=-0.35*max.range[2]
}
if (max.range[1]<0 & max.range[2] < -0.25*max.range[1]){
  max.range[2] = -0.25*max.range[1]
}

## color assignment
colon.n = str_count(r.name, pattern=':')
col.n = max(colon.n)+2
if (col.n > 12) {
  col.n = 12
}
col1 = brewer.pal(col.n, 'Set3')
col.x = col1[(colon.n %% 11) +2]
if (treat.c == 'y'){
  col.x[length(col.x)]=col1[1]
}

## p-value correction and formatting
p.vala = p.adjust(out.t1[,'p.value'])
p.valf = sprintf('%1.1E',p.vala)
p.valc = rep('blue', length(p.vala))
p.valc[p.vala < 0.05] = 'red'

## output the barplot
height = nrow(out.t1) /20 * 6.5 +2

pdf(plot.n, width=6.5,height=height)
opar = par( mar=c(5,6,2,2) )
y = barplot(mean.v, horiz=T, las=2, col=col.x,
            xlim = c(max.range[1],max.range[2]),
            xlab = 'Contribution level')
box()
text(max.range[1]*1.03,y,p.valf,cex=0.6,col=p.valc, pos=4)
arrows(error.m, y, error.p, y, angle=90,
       length = 0.08, code=3)
par(opar)
dev.off()

##############
##### Fig 3
rm(list=ls())

### packages
library(RColorBrewer)

### import data
out.dir = 'outputs.fig3/'
files.ind = dir(out.dir)
out.files = files.ind[grep('output.txt', files.ind)]

search.w = c('add','nr','Av\\.')

lm.out.list = list()
lm.out.list = lapply(search.w, function(s.word) {
  f.names = sort(out.files[grep(s.word, out.files)])
  n.int.dat = lapply(f.names, function(f.n) {
    read.delim(paste0(out.dir,f.n), header=T, row.names=1)
  })
  names(n.int.dat) = f.names
  return(n.int.dat)
})
search.w = c(search.w[1:2], 'Av')
names(lm.out.list) = search.w

coef.names = rownames(lm.out.list[[1]][[4]])
coef.names = coef.names[grep('^m\\.', coef.names)]
model.names = paste(rep(search.w, each=4), rep(1:4, times=3), sep='.')

lwr.b.mat = upr.b.mat = matrix(NA, nrow=16, ncol=3*4)
dimnames(lwr.b.mat) = dimnames(upr.b.mat) = 
  list(coef.names, model.names)

for (m.n in model.names) {
  m.n.split = unlist(strsplit(m.n, '\\.'))
  ce.tab = lm.out.list[[m.n.split[1]]][[as.integer(m.n.split[2])]]
  ce.tab.n = rownames(ce.tab)
  ce.tab.n = ce.tab.n[ce.tab.n %in% coef.names]
  ce.tab.sel = ce.tab[ce.tab.n,]
  t.95ci = qt(0.975, df=ce.tab[1,'DF'])
  c.int = t(apply(ce.tab.sel, 1, function(x) {
    mu = x['Value']
    se = x['Std.Error']
    ub = mu + se * t.95ci
    lb = mu - se * t.95ci
    return(c(lb,ub))
  }))
  lwr.b.mat[ce.tab.n, m.n] = c.int[ce.tab.n, 1]
  upr.b.mat[ce.tab.n, m.n] = c.int[ce.tab.n, 2]
}

### plotting
lwr.b.mat = -lwr.b.mat
upr.b.mat = -upr.b.mat
lwr.b.mat = rbind(lwr.b.mat[16,],lwr.b.mat[-16,])
upr.b.mat = rbind(upr.b.mat[16,],upr.b.mat[-16,])
y.min.max = range(c(lwr.b.mat, upr.b.mat), na.rm = T)
y.min.max
col1 = c('green3','blue','red','black')
mod.n = c('Additive model','NR model','Averaging model')

pdf(paste0('outputs/', 'Fig3.drop.high.orders.v3.pdf'), width = 7, height = 10)

opar = par(mfrow=c(3,1), mar=c(1,6,0.5,0.5), oma=c(4,2,0,0))
y.val.l = list()
y.val = c('Intercept','J','E','P','S','J:E','J:P','J:S','E:P','E:S','P:S',
          'J:E:P','J:E:S','J:P:S','E:P:S','J:E:P:S')
y.val.l[[1]] = factor(y.val, levels=rev(y.val))

y.val = c('Intercept','J','E','P','S','J:E','J:P','J:S','E:P','E:S','P:S',
          'J;E;P','J;E;S','J;P;S','E;P;S','J;E;P;S')
y.val.l[[2]] = factor(y.val, levels=rev(y.val))

y.val = c('Intercept','J','E','P','S','J;E','J;P','J;S','E;P','E;S','P;S',
          'J;E;P','J;E;S','J;P;S','E;P;S','J;E;P;S')
y.val.l[[3]] = factor(y.val, levels=rev(y.val))

for (s.w in 0:2) {
  stripchart(rep(0,16) ~ factor(y.val.l[[s.w+1]]), type='n',
             xlim=c(y.min.max[1], y.min.max[2]+1), ylim = c(0.6, 16.7),
             xaxt='n', ylab=NA, 
             las = 2,
             cex.axis=1.3, cex.lab=1.5)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray87")
  rect(par("usr")[1], 1.5,par("usr")[2], 5.5,col = "gray80" , border=NA)
  rect(par("usr")[1], 11.4,par("usr")[2], 15.36,col = "gray80" , border=NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
  abline(v=0, lty=2, lwd=2, col='gray45')
  abline(v=c(-2,-1,1,2,3), lty=3, col='gray60')
  for (cn in 1:4) {
    col.numb = s.w * 4 + cn
    y.offset = 0.2*(cn-1)-0.5
    segments(lwr.b.mat[,col.numb], 16:1 +y.offset,  
             upr.b.mat[,col.numb], 16:1 +y.offset, 
             col=col1[cn], lwd=2)
    points((lwr.b.mat[,col.numb]+upr.b.mat[,col.numb])/2, 16:1 +y.offset, pch=16,
           col=col1[cn])
  }
  text(y.min.max[1],16,mod.n[s.w+1],pos=4, cex=1.4)
  if (s.w == 0) {
    legend('topright',c('Full model','4-gene interaction omitted',
                        '3,4-gene interaction omitted','2,3,4-gene interaction omitted'),
           col=rev(col1), lwd=2, bg='gray90')
  }
  if (s.w==2) {
    axis(side=1, cex.axis=1.4)
    mtext(side = 1, text = expression('Contribution to Immunity [' ~ log[10] ~ '(cfu /' ~ cm^2 ~ ') ]'), line = 3)
    
  } else {
    axis(side=1,labels=F)
  }
}

par(opar)

dev.off()

