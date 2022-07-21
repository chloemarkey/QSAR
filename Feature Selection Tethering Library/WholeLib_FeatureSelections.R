# LOAD ALL NECESSARY LIBRARIES -----------------------------------------------------------------
library(caret)
library(glmnet)
library(ModelMetrics)
library(pls)
library(randomForest)
library(kernlab)
library(earth)
library(prospectr)
library(olsrr)

# LOAD WHOLE LIBRARY DATA --------------------------------------------------------------------
setwd('/Volumes/GoogleDrive/My Drive/UROP/UROP - Chloe/Paper Writing/FinalCommentedCodes/Feature Selection Tethering Library')
dat <- read.csv('/Volumes/GoogleDrive/My Drive/UROP/UROP - Chloe/Paper Writing/FinalCommentedCodes/Feature Selection Tethering Library/WholeLibrary_AllDescUPDATED.csv')
colnames(dat)[1] <- "ID"
colnames(dat)[3] <- "y"
set.seed(1234)

# IDENTIFY AND REMOVE INVARIANT DESCRIPTORS ----------------------------------------------
set.seed(1234)

# Create a new data set with invariant descriptors removed
all_desc <- colnames(dat)
invar_desc_removed_dat <- dat[vapply(dat, function(z) length(unique(z))>1, logical(1L))]

# Find column indices corresponding to invariant descriptors
reduced_desc <- colnames(invar_desc_removed_dat)
invariant_desc <- setdiff(all_desc, reduced_desc)
invariant_inds <- integer()

for (i in 1:length(invariant_desc)){
  index <- which(colnames(dat) == invariant_desc[i])
  invariant_inds <- c(invariant_inds, index)
}

# Save remaining descriptors as a model matrix
x_desc = model.matrix(y ~.-ID-SMILES, data = invar_desc_removed_dat)
x_desc = x_desc[,-1]

# Save molecule identifier information and Tethering (y) values
ID = invar_desc_removed_dat$ID
SMILES = invar_desc_removed_dat$SMILES
y = as.numeric(invar_desc_removed_dat$y)

# IDENTIFY AND REMOVE PERFECTLY CORRELATED DESCRIPTORS -----------------------------------------------------

# Generate matrix containing correlations between each descriptors and Tethering
wholeLib_desc_tethering_cormat = cor(x_desc, y, method = "kendall")

# Generate matrix containing correlations between all descriptors and each other
wholeLib_descriptor_cormat = cor(x_desc, method = "kendall")
wholeLib_descriptor_cormat = round(wholeLib_descriptor_cormat, 3)

# Set diagonal and upper triangle of matrix to 0 to prevent repeat analysis
wholeLib_descriptor_cormat[upper.tri(wholeLib_descriptor_cormat)] <- 0
diag(wholeLib_descriptor_cormat) <- 0

# Initialize empty list to store indices of descriptors to be removed
ind_to_remove <- list()

# Iterate through matrix to identify perfectly correlated descriptors pairs, removing descriptor 
# of pair with lower correlation to Tethering
for (i in 1:ncol(wholeLib_descriptor_cormat)){
  for (j in 1:nrow(wholeLib_descriptor_cormat)){
    if (abs(wholeLib_descriptor_cormat[j,i]) == 1.000){
      if (wholeLib_desc_tethering_cormat[i] > wholeLib_desc_tethering_cormat[j]){
        ind_to_remove <- c(ind_to_remove, j)
      }
      if (wholeLib_desc_tethering_cormat[j] > wholeLib_desc_tethering_cormat[i]){
        ind_to_remove <- c(ind_to_remove, i)
      }
      if (wholeLib_desc_tethering_cormat[j] == wholeLib_desc_tethering_cormat[i]){
        ind_to_remove <- c(ind_to_remove,j)
      }
    }
  }
}

# Generate ordered list of descriptor indicies removed via analysis
new_ind_to_remove <- unique(ind_to_remove)
new_ind_to_remove <- unlist(new_ind_to_remove)
new_ind_to_remove <- sort(new_ind_to_remove)

# Save list of descriptors removed via analysis
removed_descriptors <- reduced_desc[new_ind_to_remove]

# IDENTIFY AND REMOVE DESCRIPTORS CORRELATED > |0.8| -------------------------------------------------

# Update correlation matrices with new descriptors removed
x_desc = x_desc[,-c(new_ind_to_remove)]
wholeLib_desc_tethering_cormat = cor(x_desc, y, method = "kendall")
wholeLib_descriptor_cormat = cor(x_desc, method = "kendall")
wholeLib_descriptor_cormat = round(wholeLib_descriptor_cormat, 3)
wholeLib_descriptor_cormat[upper.tri(wholeLib_descriptor_cormat)] <- 0
diag(wholeLib_descriptor_cormat) <- 0

# Update list of descriptor names with new descriptors removed
descriptor_list <- colnames(wholeLib_descriptor_cormat)

# Generate a vector of descriptors to remove to reduce pair-wise correlations
built_in_inds_to_remove_80 <- findCorrelation(wholeLib_descriptor_cormat, cutoff = 0.800)

# Initialize empty list to store indices of descriptors to be removed
ind_to_remove_80 <- list()

# Initialize variable to save correlated pair information
d = NULL

# Iterate through matrix to identify correlated descriptors pairs, removing descriptor 
# of pair with lower correlation to Tethering
for (i in 1:ncol(wholeLib_descriptor_cormat)){
  for (j in 1:nrow(wholeLib_descriptor_cormat)){
    if (abs(wholeLib_descriptor_cormat[j,i]) >= 0.800 & !(j %in% ind_to_remove_80) & !(i %in% ind_to_remove_80)){
      if (abs(wholeLib_desc_tethering_cormat[i]) > abs(wholeLib_desc_tethering_cormat[j]) & !(j %in% ind_to_remove_80)){
        ind_to_remove_80 <- c(ind_to_remove_80, j)
      }
      if (abs(wholeLib_desc_tethering_cormat[j]) > abs(wholeLib_desc_tethering_cormat[i]) & !(i %in% ind_to_remove_80)){
        ind_to_remove_80 <- c(ind_to_remove_80, i)
      }
      if (abs(wholeLib_desc_tethering_cormat[j]) == abs(wholeLib_desc_tethering_cormat[i]) & !(i %in% ind_to_remove_80) & !(j %in% ind_to_remove_80)){
        if (j %in% built_in_inds_to_remove_80 & !(i %in% built_in_inds_to_remove_80)){
          ind_to_remove_80 <- c(ind_to_remove_80, j)
        }
        if (i %in% built_in_inds_to_remove_80 & !(j %in% built_in_inds_to_remove_80)){
          ind_to_remove_80 <- c(ind_to_remove_80, i)
        }
        else{
          ind_to_remove_80 <- c(ind_to_remove_80, i, j)
        }
      }
      d = rbind(d, data.frame(colnames(wholeLib_descriptor_cormat)[i], wholeLib_desc_tethering_cormat[i], rownames(wholeLib_descriptor_cormat)[j], wholeLib_desc_tethering_cormat[j]))
    }
  }
}

# Generate ordered list of descriptor indicies removed via analysis
new_ind_to_remove_80 <- unique(ind_to_remove_80)
new_ind_to_remove_80 <- unlist(new_ind_to_remove_80)
new_ind_to_remove_80 <- sort(new_ind_to_remove_80)

# Save list of descriptors removed via analysis
removed_descriptors_80 <- descriptor_list[new_ind_to_remove_80]

# Save list of descriptors remaining post analysis and generate updated data set
remaining_descriptors_80 <- x_desc[,-c(new_ind_to_remove_80)]
updated_data <- dat[,c("ID", "SMILES", "y", colnames(remaining_descriptors_80))]

# Save and export csv file of highly correlated descriptor pairs
names(d) <- c('Var1', 'V1 Cor to Tethering', 'Var2', 'V2 Cor to Tethering')
write.csv(d, file = "WholeLibCorrelatedDescriptorPairs.csv")

# SCALE DESCRIPTORS -------------------------------------------------------------------------

# Generate model matrix containing remaining descriptors
x_scaled = model.matrix(y ~ .-ID-SMILES, data = updated_data)

# Scale all descriptor values in model matrix
x_means <- colMeans(x_scaled)
x_sd <- apply(x_scaled, 2, sd)
for (i in 1:(ncol(x_scaled))){
  for (j in 1:(nrow(x_scaled))){
    scaled_val <- (x_scaled[j,i]-x_means[i])/(x_sd[i])
    x_scaled[j,i] <- scaled_val
  }
}
x_scaled[,1] = 1

# LASSO FEATURE SELECTION -------------------------------------------------------------------------------------------

# Fit a cross validated lasso model between scaled descriptors and Tethering
fit = glmnet(x_scaled, y, alpha = 1, standardize = FALSE)
cvfit = cv.glmnet(x_scaled, y, alpha = 1, keep = TRUE, standardize = FALSE, nfolds = length(y))

# Identify descriptors associated with the lasso model with a lambda value 1 standard error from the minimum 
beta <- coef(cvfit, s = "lambda.1se")
ii <- which(beta!=0)
lambda_1se_descriptors <- beta[ii, 1]
lambda_1se_chosen_desc <- names(lambda_1se_descriptors)
lambda_1se_chosen_desc <- lambda_1se_chosen_desc[-1]

# Export descriptor names to csv file
write.csv(lambda_1se_chosen_desc, file = "WholeLib_Lasso1seChosenDesc.csv")

# MARS FEATURE SELECTION --------------------------------------------------------------------

# Fit a MARS model between scaled descriptors and Tethering 
earth.mod <- earth(x_scaled, y)

# Function to identify names of descriptors used in MARS model
get.used.pred.names <- function(obj) # obj is an earth object
{
  any1 <- function(x) any(x != 0)    # like any but no warning if x is double
  names(which(apply(obj$dirs[obj$selected.terms, , drop=FALSE], 2, any1)))
}

# Identify descriptors used in MARS model
mars_chosen_desc <- get.used.pred.names(earth.mod)

# Export descriptor names to csv file
write.csv(mars_chosen_desc, file = "WholeLib_MARSchosenDesc.csv")

# STEPWISE FEATURE SELECTION ----------------------------------------------------

# Fit a linear model between scaled descriptors and Tethering 
new_linear_dat <- cbind.data.frame(y, x_scaled[,-1])
model <- lm(y ~ ., data = new_linear_dat)

# Apply stepwise selection to linear model
WholeLib_STEPfit <-ols_step_both_p(model, penter = 0.05, details = TRUE)

# Identify descriptors selected by stepwise
WholeLib_Step_ChosenDesc <- WholeLib_STEPfit$predictors

# Export descriptor names to csv file
write.csv(WholeLib_Step_ChosenDesc, file = "WholeLib_StepwiseSelectionDesc.csv")