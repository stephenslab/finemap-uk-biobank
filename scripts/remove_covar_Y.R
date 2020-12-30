library(data.table)
pheno = fread('/project2/mstephens/yuxin/ukb-bloodcells/bloodcells.pheno.txt')

covar = fread('/project2/mstephens/yuxin/ukb-bloodcells/bloodcells.covar.txt')

Y = as.matrix(pheno[,3:18])
Z = as.matrix(cbind(1, covar[,3:140]))

A   <- crossprod(Z) # Z'Z
# chol decomposition for (Z'Z)^(-1)
R = chol(solve(A)) # R'R = (Z'Z)^(-1)
RZt = tcrossprod(R, Z)
RZtY = R %*% crossprod(Z, Y) # RZ'y

Y_s = Y - crossprod(RZt, RZtY)

Y_s = cbind(pheno[,1:2], Y_s)

fwrite(Y_s, '/project2/mstephens/yuxin/ukb-bloodcells/bloodcells.pheno.resid.txt', quote=FALSE, row.names = FALSE)
