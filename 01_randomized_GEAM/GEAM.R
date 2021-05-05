# ============================================================================
# GEAM SCRIPT ================================================================
# ============================================================================
# GEAM script generates the intra and inter randomized datasets and apply the 
# GEAM (Genotype Environment Association Model) as a Stepwise Regression.


# DATA -----------------------------------------------------------------------

# Table of the 21 SNPs under selection by individual
snps <- read.table(".../SNPs_outliers_coded.txt", header = T, sep = "\t")
# Scores of the first axis of the environmental PCA by geographic cell
sc <- read.table(".../scores_pc1_all.txt", header = T, sep = "\t")


# RANDOMIZATIONS --------------------------------------------------------------

# Random_dist assigns geographic cells to individuals by shuffling the cells of
# all populations (type= Inter) or the cells of its own populations 
# (type= Intra). Finally, it adds the 21 SNPs under selection as columns. 

Random_dist <- function(table_snps, table_sc, type = "Intra"){ 
  
  l_dist <- list()
  cells <- table_sc[table_sc$POP != 0, c(1, 2)]
  
  for(i in 1:1000){
    l_dist[[i]] <- as.data.frame(matrix(ncol = 4, nrow = nrow(table_snps)))
    l_dist[[i]][, c(1, 2)] <- table_snps[, c(1, 2)]
    colnames(l_dist[[i]]) <- c("Taxa", "POP", "Cell", "PC1")
    
    for(y in 1:length(unique(table_snps$POP))){
      if(type == "Intra"){
        pop <- unique(table_snps$POP)[y]
        cells <- table_sc[table_sc$POP == pop, c(1, 2)] 
      }
      length_rcells <- sum(table_snps$POP == pop)
      random_cells <- cells[sample(1:nrow(cells), length_rcells, replace = T), ]
      l_dist[[i]][l_dist[[i]]$POP == pop, c(3, 4)] <- random_cells
    }
    l_dist[[i]] <- l_dist[[i]][,-2]
    
    l_dist[[i]] <- merge(l_dist[[i]], table_snps[, -2], by = "Taxa") 
    l_dist[[i]][, 2:24] <- apply(l_dist[[i]][, 2:24], 2, 
                                 function(x){x <- as.numeric(as.character(x))})
  }
  return(l_dist)
}

dist_intra <- Random_dist(snps, sc, type="Intra")
dist_inter <- Random_dist(snps, sc, type="Inter")


# STEPWISE REGRESSION ----------------------------------------------------------
# lmp is a function to obtain the p-values, obtained from: 
# https://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html

lmp <- function (modelobject) { 
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

GEAM <- function(n, dist_data, label){
  
  # res is the result-table to save the p-value and R squared of the model. 
  # Each row contains the regression results of an iteration of the loop,
  # i.e. of a randomized distribution table.
  res <- as.data.frame(matrix(ncol = 3, nrow=n)) 
  colnames(res) <- c("p-value", "adj.R.squared", "Label")
  res$Label <- c(rep(label, n))
  
  # res tables are included in the list l_res, each element
  # contains results from one step of the while.
  list_res <- list()
  
  # coeffs is the result-table to save coefficients of the model for each SNP.
  coeffs <- list()  
  
  y = 0
  bad = 1
  
  # The while is ON while there are still SNPs with p value > 0.05 to remove
  while(length(bad) > 0){ 
    y = y + 1
    
    list_res[[y]] <- res
    coeffs[[y]] <- list()
    
    for(x in 1:4){
      # The nº of columns of coeffs is the nº of SNPs plus the intercept.
      num_col <- ncol(dist_data[[1]])-2
      coeffs[[y]][[x]] <- as.data.frame(matrix(nrow = 1000, ncol = num_col))
      names_col <- colnames(dist_data[[1]])[4:ncol(dist_data[[1]])]
      colnames(coeffs[[y]][[x]]) <- c("Interc", names_col)
    }
    
    # Loop that applies a regression analysis to each dataset.
    for(i in 1:length(dist_data)){ 
      
      # we leave only the predictor and the dependent variables
      data <- dist_data[[i]][, c(3, 4:ncol(dist_data[[i]]))] 
      
      # Multiple regression
      pc <- lm(formula = PC1 ~ ., data=data)
      
      # Save results
      list_res[[y]][0 + i, "adj.R.squared"] <- summary(pc)$adj.r.squared 
      list_res[[y]][0 + i, "p-value"] <- lmp(pc) 
      coeffs[[y]][[1]][i, ] <- t(as.data.frame(summary(pc)$coefficients))[1, ]
      names(coeffs[[y]])[1] <- "Estimate" 
      coeffs[[y]][[2]][i, ] <- t(as.data.frame(summary(pc)$coefficients))[2, ]
      names(coeffs[[y]])[2] <- "Std Error" 
      coeffs[[y]][[3]][i, ] <- t(as.data.frame(summary(pc)$coefficients))[3, ]
      names(coeffs[[y]])[3] <- "t value"
      coeffs[[y]][[4]][i, ] <- t(as.data.frame(summary(pc)$coefficients))[4, ]
      names(coeffs[[y]])[4] <- "p value"
      
      print(c(y, i))
    }
    
    mean_pval_SNPs <- apply(coeffs[[y]][[4]], 2, mean)
    
    # We sequentially remove the least-fitted SNPs:
    # First, the ones with p-value > 0.5
    bad <- names(mean_pval_SNPs[mean_pval_SNPs > 0.5]) 
    # If there is none, the ones with p-values > 0.15
    if(length(bad) == 0){ 
      bad <- names(mean_pval_SNPs[mean_pval_SNPs > 0.15])
    }
    # If there is none, the ones with p-values > 0.05
    if(length(bad) == 0){ 
      bad <- names(mean_pval_SNPs[mean_pval_SNPs > 0.05])
    }
    if(length(bad) > 0){ 
      bad <- which(colnames(dist_data[[i]]) %in% bad) 
      dist_data <- lapply(dist_data, function(x){
        x <- x[, -c(bad)]
      })
    }
    
  }
  return(list("p&R_model" = list_res, "Coeffs" = coeffs))
}

resu.intra <- GEAM(1000, dist_intra, "pc1_intra")
resu.inter <- GEAM(1000, dist_inter, "pc1_inter")


## Results --------------------------------------------------------------------

# Mean of R squared and p value of the best model among iterations 
# (the one of the last step in the while)
mean(resu.intra[[1]][[6]][, "adj.R.squared"]) 
mean(resu.intra[[1]][[6]][, "p-value"])
mean(resu.inter[[1]][[6]][, "adj.R.squared"])
mean(resu.inter[[1]][[6]][, "p-value"])

# Means of coefficients of the best intra-population model
means <- as.data.frame(rbind(apply(resu.intra[[2]][[6]][[1]], 2, mean), 
                             apply(resu.intra[[2]][[6]][[2]], 2, mean),
                             apply(resu.intra[[2]][[6]][[3]], 2, mean), 
                             apply(resu.intra[[2]][[6]][[4]], 2, mean))) 

write.table(resu.intra[[2]][[6]][[1]], "Estimate_bestmodel_pc1_intra.txt",
            sep="\t", row.names=F)
write.table(resu.inter[[2]][[6]][[1]], "Estimate_bestmodel_pc1_inter.txt",
            sep="\t", row.names=F)
rownames(means) <- c("Estimates", "Std. Error", "t value", "p value")
write.table(means, "...", row.names=F, sep="\t")



