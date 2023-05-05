
################################
##### Fig 5
##### redo analysis with averaging model
#### one by one
############ Tsuda et al. data analysis
rm(list=ls())
### package
library(nlme)
library(RColorBrewer)
library(stringr)
library(R.utils)

##### function to make a matrix for averaging model
make.rec.mx = function(all.g = LETTERS[1:4], sub.n = LETTERS[1:4], treat.c='y') {
  dec.n = 2:(2^length(sub.n)) -1
  bin.n = intToBin(dec.n)
  bin.n = strsplit(bin.n,'')
  coef.ei = sapply(bin.n, function(x) {
    paste(sub.n[as.logical(as.integer(x))], collapse = ':')
  })
  coef.ei = sort(coef.ei)
  coef.ei = coef.ei[order(nchar(coef.ei))]  # variables in the model
  
  add.m = sapply(coef.ei, function(y) {
    y1 = unlist(strsplit(y,':'))
    if (length(y1) == 1) {
      y1 = as.integer(grepl(y1, all.g))
    } else {
      y1 = sapply(y1, function(z) {
        as.integer(grepl(z, all.g))
      })
      y1 = apply(y1, 1, prod)
    }
    return(y1)
  })
  rownames(add.m) = all.g
  
  rec.mx = add.m
  co.count = sapply(colnames(rec.mx), str_count, pattern=':')
  for (c.c in 0:max(co.count)) {
    ap.cols = co.count == c.c
    if (sum(ap.cols) == 1) next
    rec.mx1 = rec.mx[,ap.cols]
    rec.mx[,ap.cols] = t(apply(rec.mx1, 1, function(x) {
      one.c = sum(x)
      row.v = x
      if (one.c > 1) {
        non.0 = 1/sum(x)
        row.v[x==1] = non.0
      } 
      return(row.v)
    }))
  }
  if (treat.c == 'y'){
    rec.mx = cbind(rec.mx, remainder=1) ## rec.m in net reconst was changed to rec.mx
  }
  return(rec.mx)  ## rec.mx is the matrix for averaging model
}

############# fit averaging model
### input file names
i.files = c('AvrRpt2_ETI.txt',
            'AvrRpm1_ETI.txt',
            'flg22_PTI.txt',
            'elf18_PTI.txt')
dir.name = './Supplemental.Dataset1/'
out.name = './outputs/'
out.file.name = c()

keep.gene.l = list()

for (file.numb in 1:length(i.files)) {
  ### load data
  input.fn = i.files[file.numb]
  input.d.n = paste0(dir.name, input.fn)
  input.dat = read.delim(input.d.n, header=T, colClasses='character')
  
  ### fix the input.dat 1
  if (file.numb == 1) {
    input.dat = input.dat[input.dat$genotype != 'rpm1/rps2' & input.dat$genotype != 'npr1',]
    colnames(input.dat)[colnames(input.dat) == 'colony'] = 'value'
    input.dat$treatment[input.dat$treatment == 'avrRpt2'] = 't'
    input.dat$treatment[input.dat$treatment == '_pLAFR'] = 'c'
  }
  if (file.numb == 2) {
    input.dat = input.dat[input.dat$genotype != 'rpm1/rps2' & input.dat$genotype != 'npr1',]
    colnames(input.dat)[colnames(input.dat) == 'colony'] = 'value'
    input.dat$treatment[input.dat$treatment == 'avrRpm1'] = 't'
    input.dat$treatment[input.dat$treatment == '_pLAFR'] = 'c'
  }
  if (file.numb == 3) {
    input.dat = input.dat[input.dat$genotype != 'fls2',]
    colnames(input.dat)[colnames(input.dat) == 'colony'] = 'value'
    input.dat$treatment[input.dat$treatment == 'flg22'] = 't'
    input.dat$treatment[input.dat$treatment == '_mock'] = 'c'
    colnames(input.dat)[1] = 'replicate'
  }
  if (file.numb == 4) {
    input.dat = input.dat[input.dat$genotype != 'fls2' & input.dat$genotype != 'efr',]
    colnames(input.dat)[colnames(input.dat) == 'colony'] = 'value'
    input.dat$treatment[input.dat$treatment == 'elf18'] = 't'
    input.dat$treatment[input.dat$treatment == '_mock'] = 'c'
    colnames(input.dat)[3] = 'pot'
  }
  
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
  
  ### v.3, check the significance of each sector
  ## the model with treatment
  if (file.numb == 5) {
    lme.pheno = lme( value ~ genotype -1,
                   random =~ 1|replicate,
                   data = input.dat)
  } else {
    lme.pheno = lme( value ~ genotype/treatment -1 + genotype,
                   random =~ 1|replicate/flat/pot,
                   data = input.dat)
  }
  
  if (file.numb == 5) subt = 16 else subt = 0
  coefs = lme.pheno$coefficients[['fixed']]
  coefs = coefs[17:32-subt]
  
  geno.names = names(coefs)
  geno.names = substr(geno.names, 9, 12)
  
  vcov = lme.pheno$varFix
  vcov = vcov[17:32-subt, 17:32-subt]
  
  names(coefs) = rownames(vcov) = colnames(vcov) = geno.names
  df = summary(lme.pheno)$tTable[17-subt, 3]
  
  ## for each gene, p.val calculation for the differences A - a
  gene.n = LETTERS[1:4]
  
  p.val.mat = c()
  for (g.n in gene.n) {
    geno.nA = geno.names[grep(g.n, geno.names)]
    str1 = paste0('^(.*)',g.n,'(.*)$', collapse = '')
    str2 = paste0('\\1',tolower(g.n),'\\2', collapse = '')
    geno.na = sub(str1, str2, geno.nA, perl=T)
    geno.Aa = rbind(geno.nA, geno.na)
    p.vals = apply(geno.Aa, 2, function(x) {
      mean.diff = coefs[x[1]] - coefs[x[2]]
      semd = sqrt(vcov[x[1],x[1]] + vcov[x[2],x[2]] - 2* vcov[x[1], x[2]]) # corrected v.4
      t.val = abs(mean.diff)/semd
      p.val = 2*pt(t.val, df, lower.tail = F)
      return(p.val)
    })
    p.val.mat = rbind(p.val.mat, p.vals)
  }
  # BH-FDR correction  # changed in v.4 from no correction 
  p.val.mat = matrix(p.adjust(p.val.mat, method = 'BH'), nrow=4)
  rownames(p.val.mat) = gene.n
  keep.genes = apply(p.val.mat, 1, min) < 0.05
  
  keep.gene.l[[file.numb]] = keep.genes # for record
  act.gene.numb = sum(keep.genes)

  all.g = levels(input.dat$genotype)
  sub.n = names(keep.genes)[keep.genes]
  treat.c='y'; if (file.numb ==5) treat.c = 'n'
  rec.mx = make.rec.mx(all.g=all.g, sub.n=sub.n, treat.c=treat.c)

  #### averaging model
  ##### mixed-effect linear model
  ### output file names
  var.kept = paste0(sub.n, collapse = '')
  input.o.n = paste0(out.name,input.fn)
  outp.n = paste(input.o.n, var.kept, 'Ave.v4.output.txt', sep='.')
  out.file.name = c(out.file.name, outp.n)
  boxp.n = paste(input.o.n, var.kept, 'Ave.v4.boxplot.pdf', sep='.')
  diag.n = paste(input.o.n, var.kept, 'Ave.v4.diagnostic.pdf', sep='.')
  plot.n = paste(input.o.n, var.kept, 'Ave.v4.reconstitution.pdf', sep='.')
  ### generate the matrix for the data
  m. = rec.mx[as.character(input.dat$genotype),]
  m.[input.dat$treatment=='c',] = 0
  
  n.int = n.sub = 4
  left.m = 5
  if (n.sub > 5) left.m = left.m + (n.sub-5)/3
  ### model fitting with boxplot of the data for fixed effects
  height = length(levels(input.dat$genotype)) * length(levels(input.dat$treatment)) /4
  if (height < 8) height = 8
  pdf(boxp.n,width=7,height=height)
  opar = par(mar=c(5, left.m, 2, 2))
  if (file.numb == 5) {
    boxplot(value~genotype,data=input.dat,horizontal=T,las=1,
            xlab='value', ylab=NA)
  } else {
    boxplot(value~genotype:treatment,data=input.dat,horizontal=T,las=1,
            xlab='value', ylab=NA)
  }
  par(opar)
  dev.off()
  
  ## the model with treatment
  if (file.numb == 5) {
    lme.rec = lme( value ~ m. ,
                   random =~ 1|replicate,
                   data = input.dat)
  } else {
    lme.rec = lme( value ~ m. -1 + genotype,
                   random =~ 1|replicate/flat/pot,
                   data = input.dat)
  }
  
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
  text(max.range[1]*0.95,y,p.valf,cex=0.6,col=p.valc, pos=4)
  arrows(error.m, y, error.p, y, angle=90,
         length = 0.08, code=3)
  par(opar)
  dev.off()
  
}

#################
###### plotting, Av model results
plot.dat = list()

### import output file
i.files = out.file.name
set.n = c('AvrRpt2-ETI','AvrRpm1-ETI', 'flg22-PTI', 'elf18-PTI')
names(i.files) = set.n

####
dec.n = 2:(2^4) -1
bin.n = intToBin(dec.n)
bin.n = strsplit(bin.n,'')
sub.n = LETTERS[1:4]
coef.ei = sapply(bin.n, function(x) {
  paste(sub.n[as.logical(as.integer(x))], collapse = ':')
})
coef.ei = sort(coef.ei)
coef.ei = coef.ei[order(nchar(coef.ei))]  # variables in the model
coef.ei = c('remainder', coef.ei)

coef.mat = matrix(0, nrow=length(set.n), ncol=length(coef.ei))
dimnames(coef.mat) = list(set.n, coef.ei)
confint95.mat = sem.mat = coef.mat
pval.mat = matrix(NA, nrow=length(set.n), ncol=length(coef.ei))
dimnames(pval.mat) = list(set.n, coef.ei)
  
for (d.set in set.n) {
  lm.out.file = i.files[d.set]
  lm.out.coef = read.delim(lm.out.file, header=T, row.names=1)
  
  coef.names = rownames(lm.out.coef)
  coef.names = coef.names[grep('^m\\.', coef.names)]
  coef.n2 = sub('^m\\.(.+)$', '\\1', coef.names)
  
  coef.mat[d.set, coef.n2] = lm.out.coef[paste0('m.', coef.n2), 'Value']
  sem.mat[d.set, coef.n2] = lm.out.coef[paste0('m.', coef.n2), 'Std.Error']
  pval.mat[d.set, coef.n2] = lm.out.coef[paste0('m.', coef.n2), 'p.value']
  df = lm.out.coef[paste0('m.', coef.n2[1]), 'DF']
  t.95ci = qt(0.975, df=df)
  confint95.mat[d.set, ] = sem.mat[d.set, ] * t.95ci
}
confint95upr.mat = coef.mat + confint95.mat
confint95lwr.mat = coef.mat - confint95.mat

for (d.set in set.n) {
  lwr.b = -confint95lwr.mat[d.set,]
  upr.b = -confint95upr.mat[d.set,]
  m.est = -coef.mat[d.set,]
  p.val = pval.mat[d.set,]
  
  ## p-value correction and formatting
  p.vala = p.adjust(p.val)
  p.valf = sprintf('%1.1E',p.vala)
  p.valf[p.valf=="NA"] = '-'
  p.valc = rep('blue', length(p.vala))
  p.valc[p.vala < 0.05] = 'red'
  
  plot.dat[[d.set]] = list(lwr.b=lwr.b, upr.b=upr.b, est=m.est, 
                                p.vala = p.vala, p.valf = p.valf, p.valc = p.valc)
  
}

### plotting
y.val = c('Intercept','J','E','P','S','J;E','J;P','J;S','E;P','E;S','P;S',
          'J;E;P','J;E;S','J;P;S','E;P;S','J;E;P;S')
y.val = factor(y.val, levels=rev(y.val))
y.min.max = c(-2, 2)
col1 = c('green3', 'red', 'blue','orange', 'black')
col1 = rep(col1, c(1,4,6,4,1))
imm.t = c('AvrRpt2-ETI', 'AvrRpm1-ETI', 'flg22-PTI', 'elf18-PTI')

pdf('./outputs/reanalyzed.w.average.model.wremove.v4.pdf', width = 12, height = 8.5)

opar = par(mfcol=c(2,2), mar=c(1,6,0.5,0.5), oma=c(4,2,0,0))
for (file.numb in 1:4) {
  file.n = names(plot.dat)[file.numb]
  plot.d = plot.dat[[file.n]]
  
  stripchart(rep(0,16) ~ factor(y.val), type='n',
             xlim=c(y.min.max[1], y.min.max[2]), ylim = c(0.6, 17.7),
             xaxt='n', ylab=NA,
             las = 2,
             cex.axis=1.3, cex.lab=1.5)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray87")
  rect(par("usr")[1], 1.5,par("usr")[2], 5.5,col = "gray80" , border=NA)
  rect(par("usr")[1], 11.4,par("usr")[2], 15.36,col = "gray80" , border=NA)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = NA)
  abline(v=0, lty=2, lwd=2, col='gray45')
  abline(v=c(-2,-1,1,2), lty=3, col='gray60')
  if (file.numb ==5) {
    segments(plot.d$lwr.b[-1] * log10(2), 15:1,  
             plot.d$upr.b[-1] * log10(2), 15:1, lwd=2)
    points(plot.d$est[-1] * log10(2), 15:1, pch=16)
    axis(side=1, cex.axis=1.4)
    mtext(side = 1, text = expression('Contribution to Immunity [' ~ log[10] ~ '(plant DNA/fungal DNA) ]'), line = 3,
          cex=0.8)
    text(y.min.max[1]*0.98,15:1, plot.d$p.valf[-1], cex=0.7, col=plot.d$p.valc[-1], pos=4)
    text(y.min.max[1]*0.98,16,'NA', cex=0.7, pos=4)
    
  } else {
    segments(plot.d$lwr.b, 16:1,  
             plot.d$upr.b, 16:1, lwd=2)
    points(plot.d$est, 16:1, pch=16)
    text(y.min.max[1]*0.98,16:1, plot.d$p.valf, cex=0.7, col=plot.d$p.valc, pos=4)
  }
  text(y.min.max[1],17.2,imm.t[file.numb],pos=4, cex=1.4)
  if (file.numb == 2|file.numb == 4) {
    axis(side=1, cex.axis=1.4)
    mtext(side = 1, text = expression('Contribution to Immunity [' ~ log[10] ~ '(cfu /' ~ cm^2 ~ ') ]'), line = 3,
          cex=0.8)
    
  } else {
    axis(side=1,labels=F)
  }
}

par(opar)

dev.off()
