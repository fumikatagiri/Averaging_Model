#### Fig 1 additive model
z.vals = c(8,3,11,2,6,2,2,2)
names(z.vals) = c('ABC','aBC','AbC','ABc','abC','aBc','Abc','abc')

### ANOVA model, comp to ABC
anov.m.mat = matrix(c(1,0,0,0,0,0,0,0,
                      1,1,0,0,0,0,0,0,
                      1,0,1,0,0,0,0,0,
                      1,0,0,1,0,0,0,0,
                      1,1,1,0,1,0,0,0,
                      1,1,0,1,0,1,0,0,
                      1,0,1,1,0,0,1,0,
                      1,1,1,1,1,1,1,1 ), nrow=8, byrow=T)
dimnames(anov.m.mat) = list(c(names(z.vals)), c('intercept','a','b','c','a:b','a:c','b:c','a:b:c'))
anov.m.sol = solve(anov.m.mat, z.vals)
anov.m.sol


### ANOVA model, comp to abc
anov.mat = matrix(c(1,1,1,1,1,1,1,1,
                    1,0,1,1,0,0,1,0,
                    1,1,0,1,0,1,0,0,
                    1,1,1,0,1,0,0,0,
                    1,0,0,1,0,0,0,0,
                    1,0,1,0,0,0,0,0,
                    1,1,0,0,0,0,0,0,
                    1,0,0,0,0,0,0,0 ), nrow=8, byrow=T)
dimnames(anov.mat) = list(c(names(z.vals)), c('intercept','A','B','C','A:B','A:C','B:C','A:B:C'))
anov.sol = solve(anov.mat, z.vals)
anov.sol

