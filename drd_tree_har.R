# DRD-HAR tree model (not recursive)
# Author: Daniele Ballinari
# R-version: 4.0.1 (2020-06-06)
# Notes: functions starting with "." are for internal use
# Todo: print summary of a tree, e.g. splits and estimates in each leaf



#'  Main function to build the tree
#'
#' @param data data frame with columns of all variables needed for the model and the ids to identify the correlation pairs
#' @param formula defines the model in the leafs
#' @param split_variables vector of variable names used in splitting the nodes
#' @param id character name of the variable that identifies the correlation pairs
#' @param miter number of iterations when growing the tree
#' @param min_obs minimum number of observations in each leaf
#' @param min_ids minimum number of observations in each leaf for each correlation pair
#' @param mtry how many split variables should be considered when making a split (if NULL, all variables are used)
#' @param mesh number of grids when determining a split 
#' @param fit logical defining if the final models in the leafs should be estimated 
#'
#' @return a list defining the tree where the active tree is the most complex tree
#' 
build_tree <- function(data, formula, split_variables, id, miter=7, min_obs=10, min_ids=10, mtry=NULL, mesh=8, fit=TRUE) {
  # ensure that the data is in a data.table object
  data <- data.table::as.data.table(data)
  # define grid of split variables
  if (!is.null(split_variables)) {
    alphas <- (1:(mesh-1))/mesh
    data_x <- data[,..split_variables]
    data_split <- data_x[, lapply(.SD, quantile, probs = alphas, type = 1)]
    data_x[,idx___ := 1:nrow(data_x)]
  } else {
    data_x <- NULL
    data_split <- NULL
  }
  # get the variables in the formula
  formula <- as.formula(formula)
  vars <- all.vars(formula)
  # change the names of the identifier and the variables of the leaf-models
  data.table::setnames(data, id, "id___")
  data.table::setnames(data, vars[1], "y___")
  vars_select <- c("id___", "y___")
  for (i in 2:length(vars)) {
    var_i <- sprintf("x%i___", i-1)
    data.table::setnames(data, vars[i], var_i)
    vars_select <- append(vars_select, var_i)
  }
  data <- data[,..vars_select]
  data[,idx___:=1:nrow(data)]
  
  # define number of ids in the dataset
  n_ids <- length(unique(data[,id___]))
  # number of observations (for each id)
  n_obs <- nrow(data)/n_ids
  # number of parameters in the leaf model
  n_param <- length(vars) - 1
  
  # grow the tree:
  history_tree <- .grow_tree(data = data, data_x = data_x, data_split = data_split, 
                             n_ids = n_ids, min_obs = min_obs, min_ids = min_ids, fit = fit, miter = miter, mtry = mtry)
  
  
  tree <- list("history_tree"=history_tree, 
               "active_tree"=history_tree[[length(history_tree)]], 
               "n_obs"=n_obs, 
               "n_ids"=n_ids, 
               "n_param"=n_param, 
               "vars"=vars, 
               "id"=id)
  
  return(tree)
}

#' Build tree buy determining the complexity through blocked-cross-validation
#'
#' @param data data frame with columns of all variables needed for the model and the ids to identify the correlation pairs
#' @param formula defines the model in the leafs
#' @param split_variables vector of variable names used in splitting the nodes
#' @param id character name of the variable that identifies the correlation pairs
#' @param time character name of the variable that identifies the time
#' @param cv number of folds for the cross-validation
#' @param h prediction horizon (we remove h-observations from the boundaries of the validation set)
#' @param miter number of iterations when growing the tree
#' @param min_obs minimum number of observations in each leaf
#' @param min_ids minimum number of observations in each leaf for each correlation pair
#' @param mtry how many split variables should be considered when making a split (if NULL, all variables are used)
#' @param mesh number of grids when determining a split 
#' @param cluster parallel cluster to run the cross-validation in parallel (can be NULL)
#'
#' @return a list defining the tree where the active tree is the one minimizing the euclidean CV-loss
#' 
cv_build_tree <- function(data, formula, split_variables, id, time, cv=5, h=1, miter=7, min_obs=10, min_ids=10, mtry=NULL, mesh=8, cluster = NULL) {
  require(foreach)
  # ensure that the data is in a data.table object
  data <- data.table::as.data.table(data)
  # define grid of split variables
  if (!is.null(split_variables)) {
    alphas <- (1:(mesh-1))/mesh
    data_x <- data[,..split_variables]
    data_split <- data_x[, lapply(.SD, quantile, probs = alphas, type = 1)]
    data_x[,idx___ := 1:nrow(data_x)]
  } else {
    data_x <- NULL
    data_split <- NULL
  }
  # get the variables in the formula
  formula <- as.formula(formula)
  vars <- all.vars(formula)
  # change the names of the identifier and the variables of the leaf-models
  data.table::setnames(data, time, "time___")
  data.table::setnames(data, id, "id___")
  data.table::setnames(data, vars[1], "y___")
  vars_select <- c("time___", "id___", "y___")
  for (i in 2:length(vars)) {
    var_i <- sprintf("x%i___", i-1)
    data.table::setnames(data, vars[i], var_i)
    vars_select <- append(vars_select, var_i)
  }
  data <- data[,..vars_select]
  # add an index variable: this is used when splitting the data set at a specific node (data_x has the same indexing)
  data[,idx___:=1:nrow(data)]
  # define time index 
  time_index <- sort(unique(data[,time___]))
  # define number of ids in the dataset
  n_ids <- length(unique(data[,id___]))
  # number of observations (for each id)
  n_obs <- nrow(data)/n_ids
  # number of parameters in the leaf model
  n_param <- length(vars) - 1
  # define block size:
  block_size <- ceiling(length(time_index)/cv)
  
  # grow the tree using all data:
  history_tree <- .grow_tree(data = data, data_x = data_x, data_split = data_split, miter = miter,
                             n_ids = n_ids, min_obs = min_obs, min_ids = min_ids, fit = TRUE, mtry = mtry)
  # define tree
  tree <- list("history_tree"=history_tree,
               "n_obs"=n_obs, 
               "n_ids"=n_ids, 
               "n_param"=n_param, 
               "vars"=vars, 
               "id"=id)
  # if there are no split variables, we can directly return the tree
  if (is.null(split_variables)) {
    tree[["active_tree"]] <- tree$history_tree[[1]]
    return(tree)
  }
  # define lambda range: at this lambdas the tree changes complexity
  lambdas_theo <- -diff(Reduce(f = c, x = lapply(history_tree, function(x) x$tree_loss)))/n_obs/3
  # we consider then a range of tuning parameters that goes from 0 (maximal complexity) up to a 
  # parameter which is 20% larger than the tuning parameter that reduces the tree to a simple root
  lambdas <- seq(0, max(lambdas_theo)*1.2, length.out = 100)
  # define CV tree:
  tree_cv <- list("n_obs"=n_obs, 
                  "n_ids"=n_ids, 
                  "n_param"=n_param, 
                  "vars"=vars, 
                  "id"=id)
  
  # use blocked cross-validation to determine the optimal complexity 
  doParallel::registerDoParallel(cluster)
  mse <- foreach(cv_i = seq(1, cv), .combine = "cbind", 
                 .packages = c("data.table"), 
                 .export = c(".grow_tree", "prune_tree", "predict_tree", ".lambda_cost_tree", 
                             ".model_loss_drd_har", ".fit_model_drd_har", ".predict_model_drd_har", 
                             ".get_split", ".fit_tree")) %dopar% {
                               # define the dates of the validation and train datasets
                               time_index_valid <- time_index[(cv_i + (cv_i-1)*block_size):min(block_size*cv_i, n_obs)]
                               time_index_train <- time_index[!time_index %in% time_index_valid]
                               # remove h-observations from the validation set to ensure independence from the train set
                               if (cv_i > 1) time_index_valid <- tail(time_index_valid, -h)
                               if (cv_i < cv) time_index_valid <- head(time_index_valid, -h)
                               # define train and validation in term of the indexing variables
                               index_data_train <- data[time___ %in% time_index_train, idx___] 
                               index_data_valid <- data[time___ %in% time_index_valid, idx___] 
                               # grow the tree on the training blocks
                               tree_cv[["history_tree"]] <- .grow_tree(data = data[idx___ %in% index_data_train, -"time___"], 
                                                                       data_x = data_x[idx___ %in% index_data_train, ], 
                                                                       data_split = data_split, miter = miter,
                                                                       n_ids = n_ids, min_obs = min_obs, min_ids = min_ids, 
                                                                       fit = TRUE, mtry = mtry)
                               tree_cv[["n_obs"]] <- length(time_index_train)
                               # define the validation dataset
                               data_valid <- cbind(data[idx___ %in% index_data_valid,-"idx___"], data_x[idx___ %in% index_data_valid,-"idx___"])
                               mse_i <- vector(mode = "numeric", length = length(lambdas))
                               # for each tuning parameter lambda, prune the tree and forecast the validation data
                               for (j in seq_along(lambdas)) {
                                 tree_cv <- prune_tree(tree_cv, cost_function = "lambda", lambda = lambdas[j])
                                 yhat_ij <- predict_tree(tree_cv, data_valid, .data_ready=TRUE)
                                 mse_i[j] <- mean((data_valid[,y___] - yhat_ij)^2)
                               }
                               return(mse_i)
                             }
  # find optimal complexity:
  avg_mse <- rowMeans(mse)
  opt_lambda <- lambdas[which.min(avg_mse)]
  
  # prune the tree with the optimal lambda parameter obtained with the CV approach
  tree <- prune_tree(tree, cost_function = "lambda", lambda = opt_lambda)
  
  return(tree)
  
  
}

#' Function that estimates a random forest DRD-HAR model
#'
#' @param data  data frame with columns of all variables needed for the model, the ids to identify the correlation pairs and a time index
#' @param formula defines the model in the leafs
#' @param split_variables vector of variable names used in splitting the nodes
#' @param id character name of the variable that identifies the correlation pairs
#' @param time character name of the variable that identifies the time
#' @param miter number of iterations when growing the tree
#' @param mtree number of trees in the random forest
#' @param min_obs minimum number of observations in each leaf
#' @param min_ids minimum number of observations in each leaf for each correlation pair
#' @param mtry how many split variables should be considered when making a split (if NULL, all variables are used)
#' @param mesh number of grids when determining a split 
#' @param block_length integer that defines the block length used for bootstrapping the data when constructing the tree
#' @param cluster optional cluster to run the code in parallel (can be NULL)
#'
#' @return a list defining the random forest
#' 
build_random_forest <- function(data, formula, split_variables, id, time, 
                                miter=7, mtree=10, min_obs=10, min_ids=10, mtry=NULL, mesh=8, block_length=22, cluster = NULL) {
  require(foreach)
  # ensure that the data is in a data.table object
  data <- data.table::as.data.table(data)
  # define grid of split variables
  alphas <- (1:(mesh-1))/mesh
  data_x <- data[,..split_variables]
  data_split <- data_x[, lapply(.SD, quantile, probs = alphas, type = 1)]
  data_x[,idx___ := 1:nrow(data_x)]
  # get the variables in the formula
  formula <- as.formula(formula)
  vars <- all.vars(formula)
  # change the names of the identifier and the variables of the leaf-models
  data.table::setnames(data, time, "time___")
  data.table::setnames(data, id, "id___")
  data.table::setnames(data, vars[1], "y___")
  vars_select <- c("time___", "id___", "y___")
  for (i in 2:length(vars)) {
    var_i <- sprintf("x%i___", i-1)
    data.table::setnames(data, vars[i], var_i)
    vars_select <- append(vars_select, var_i)
  }
  data <- data[,..vars_select]
  # add an index variable: this is used when splitting the data set at a specific node (data_x has the same indexing)
  data[,idx___:=1:nrow(data)]
  # define time index 
  time_index <- sort(unique(data[,time___]))
  # define number of ids in the dataset
  n_ids <- length(unique(data[,id___]))
  # number of observations (for each id)
  n_obs <- nrow(data)/n_ids
  # number of parameters in the leaf model
  n_param <- length(vars) - 1
  # create random samples of dates  
  samples <- boot::tsboot(1:n_obs, function(x) x, R = n_obs, sim = "fixed", l = block_length, n.sim = mtree)$t
  
  # use blocked cross-validation to determine the optimal complexity 
  doParallel::registerDoParallel(cluster)
  forest <- foreach(tree_i = seq(1, mtree),
                    .packages = c("data.table"), 
                    .export = c(".grow_tree", "prune_tree", "predict_tree", ".lambda_cost_tree", 
                                ".model_loss_drd_har", ".fit_model_drd_har", ".predict_model_drd_har", 
                                ".get_split", ".fit_tree")) %dopar% {
                                  index_i <- table(samples[,tree_i])
                                  
                                  data_i <- NULL
                                  data_x_i <- NULL
                                  for (m in unique(index_i)) {
                                    idx_m <- data[time___ %in% time_index[as.numeric(names(index_i[index_i==m]))], idx___]
                                    for (j in 1:m) {
                                      data_i <- rbind(data_i, data[idx___ %in% idx_m, ])
                                      data_x_i <- rbind(data_x_i, data_x[idx___ %in% idx_m, ])
                                    }
                                    
                                  }
                                  
                                  # grow a tree on the bootstrapped data
                                  tree_i <- .grow_tree(data = data_i[, -"time___"], 
                                                       data_x = data_x_i, 
                                                       data_split = data_split, miter = miter,
                                                       n_ids = n_ids, min_obs = min_obs, min_ids = min_ids, 
                                                       fit = FALSE, mtry = mtry)
                                  
                                  # only keep most complex tree
                                  tree_i <- tree_i[[length(tree_i)]]
                                  
                                  # fit the final tree
                                  tree_i$tree_fits <- .fit_tree(tree_i$tree_splits, data = data_i, data_x = data_x_i)
                                  
                                  
                                  return(list("active_tree" = tree_i))
                                }
  
  forest <- list("forest"=forest, 
                 "n_obs"=n_obs, 
                 "n_ids"=n_ids, 
                 "n_param"=n_param, 
                 "vars"=vars, 
                 "id"=id)
  return(forest)
}

#' Function that makes predictions based on a fitted tree HAR-DRD model
#'
#' @param tree list defining a tree HAR-DRD model
#' @param newdata data frame with columns of all variables needed for the model
#' @param .data_ready internal variable, when equal to true, 'newdata' is assumed to be already prepared (used when cross-validating the tree-complexity)
#'
#' @return a vector of predicted correlations
#' 
predict_tree <- function(tree, newdata, .data_ready = FALSE) {
  # check that the active tree has been fitted
  if (is.null(tree$active_tree$tree_fits)) stop("Tree has not been fitted!")
  # for internal use: the data might already be prepared for the prediction 
  if (!.data_ready) {
    # ensure that the data is in a data.table object
    newdata <- data.table::as.data.table(newdata)
    # change the names of the identifier and the variables of the leaf-models
    data.table::setnames(newdata, tree$id, "id___")
    data.table::setnames(newdata, tree$vars[1], "y___")
    vars_select <- c("id___", "y___")
    for (i in 2:length(tree$vars)) {
      var_i <- sprintf("x%i___", i-1)
      data.table::setnames(newdata, tree$vars[i], var_i)
      vars_select <- append(vars_select, var_i)
    }
  }
  # Add an identifier to the data: this helps in putting the predictions in the correct order
  newdata[,"__identifier__"] <- 1:nrow(newdata)
  # define the activate tree as the tree
  tree <- tree$active_tree
  # Check if the splits defining the active tree have length 1, i.e. root tree
  if (length(tree$tree_splits) == 1) {
    predictions <- .predict_model_drd_har(newdata = newdata, param = tree$tree_fits[[1]])
  } else {
    # Initialize the vector of predictions
    predictions <- vector(mode="numeric", length = nrow(newdata))
    # For each leaf, make the predictions
    for (i in seq_along(tree$tree_splits)) {
      # Get path to leaf
      node <- tree$tree_splits[i]
      # Get parameters of the leaf
      param <- tree$tree_fits[[i]]
      # Filter out all observations that are in the current leaf
      newdata_i <- newdata[eval(parse(text=node)), ]
      # Make the predictions
      predictions[newdata_i[["__identifier__"]]] <- .predict_model_drd_har(newdata = newdata_i, param = param)
    }
  }
  
  return(predictions)
}

#' Function that makes predictions based on a fitted random forest HAR-DRD model
#'
#' @param forest a list defining the random forest HAR-DRD model 
#' @param newdata data frame with columns of all variables needed for the model
#' @param cluster optional cluster to run the code in parallel (can be NULL)
#'
#' @return a vector of predicted correlations
#'
predict_forest <- function(forest, newdata, cluster = NULL) {
  require(foreach)
  # ensure that the data is in a data.table object
  newdata <- data.table::as.data.table(newdata)
  # change the names of the identifier and the variables of the leaf-models
  data.table::setnames(newdata, forest$id, "id___")
  data.table::setnames(newdata, forest$vars[1], "y___")
  vars_select <- c("id___", "y___")
  for (i in 2:length(forest$vars)) {
    var_i <- sprintf("x%i___", i-1)
    data.table::setnames(newdata, forest$vars[i], var_i)
    vars_select <- append(vars_select, var_i)
  }
  forest <- forest$forest
  
  doParallel::registerDoParallel(cluster)
  yhat <- foreach(forest_i = forest, .combine = "+", .packages = "data.table", .export = c("predict_tree", ".predict_model_drd_har")) %dopar% {
    return(predict_tree(tree = forest_i, newdata = newdata, .data_ready = TRUE))
  }
  yhat <- yhat/length(forest)
  return(yhat)
}

#' Function that prunes a tree HAR-DRD model
#'
#' @param tree list defining a tree HAR-DRD model
#' @param cost_function string that specifies which complexity cost should be used when pruning the tree ("aic", "bic", or "lambda")
#' @param lambda complexity cost parameter, higher values lead to less complex trees
#'
#' @return list defining a tree HAR-DRD model where the active tree is the pruned tree
#' 
prune_tree <- function(tree, cost_function = "aic", lambda = NULL) {
  if (cost_function == "aic") {
    tree_cost <- Reduce(f = c, x = lapply(tree$history_tree, .ic_cost_tree, k=2, n_obs = tree$n_obs, n_param = tree$n_param))
  } else if (cost_function == "bic") {
    tree_cost <- Reduce(f = c, x = lapply(tree$history_tree, .ic_cost_tree, k=log(tree$n_obs), n_obs = tree$n_obs, n_param = tree$n_param))
  } else if (cost_function == "lambda") {
    tree_cost <- Reduce(f = c, x = lapply(tree$history_tree, .lambda_cost_tree, lambda=lambda, n_obs = tree$n_obs, n_param = tree$n_param))
  } else {
    stop("Not supported cost function for pruning!")
  }
  
  tree[["active_tree"]] <- tree[["history_tree"]][[which.min(tree_cost)]]
  
  return(tree)
}


#' Function that computes the variable importance for a random forest HAR-DRD
#'
#' @param forest  a list defining the random forest HAR-DRD model 
#'
#' @return a named vector with the average variable importance across all trees in the forest
#'
variable_importance_forest <- function(forest) {
  var_imp <- Reduce(f = "rbind", x = lapply(forest$forest, function(x) x$active_tree$variable_importance))
  return(colMeans(var_imp))
}



# Internal functions ------------------------------------------------------

# Evaluate the loss of the model
.model_loss_drd_har <- function(data) {
  data <- data.table::copy(data)
  data[, m:=mean(x1___), by=id___]
  data[, `:=`(y___ = y___ - m, x2___ = x2___ - m, x3___ = x3___ - m, x1___ = x1___ - m) ]
  data <- as.matrix(data[,.(y___, x1___, x2___, x3___)])
  
  model <- speedglm::speedlm.fit(y = data[,1], X = data[,-1], intercept = FALSE, sparse = FALSE)
  
  sse <- model$RSS
  
  if (any(model$coef<0) | sum(model$coef) >= 1 ) {
    sse <- Inf
  }
  
  return(sse)
}

# Fit the HAR-DRD model on a data frame
.fit_model_drd_har <- function(data) {
  data <- data.table::copy(data)
  data[, m:=mean(x1___), by=id___]
  avg <- base::unique(data[,.(id___, m)])
  data[, `:=`(y___ = y___ - m, x2___ = x2___ - m, x3___ = x3___ - m, x1___ = x1___ - m) ]
  data <- as.matrix(data[,.(y___, x1___, x2___, x3___)])
  
  reg <- speedglm::speedlm.fit(y = data[,1], X = data[,-1], intercept = FALSE, sparse = FALSE)
  
  model <- list("coef"=reg$coef, "m"=avg)
  
  return(model)
}

# Make predictions based on a fitted HAR-DRD model
.predict_model_drd_har <- function(newdata, param) {
  means <- param$m
  coef <- param$coef
  model_matrix <- dplyr::left_join(newdata, means, by = "id___")
  model_matrix <- as.matrix(model_matrix[,.(m, x1___, x2___, x3___)])
  coef <- as.matrix(c(1-sum(coef), coef), ncol=1)
  y_hat <- model_matrix %*% coef
  return(as.numeric(y_hat))
}

# Function that grows the tree and returns the final tree and the history of trees for each growing-step
.grow_tree <- function(data, data_x, data_split, n_ids, min_obs, min_ids, fit, miter, mtry) {
  if (is.null(data_split)) {
    if (fit) tree_fits <- list(.fit_tree(tree_splits = "", data = data, data_x = NULL))
    history_tree <- list(list("tree_splits" = "", "tree_losses" = .model_loss_drd_har(data), 
                              "tree_fits" = tree_fits))
    return(history_tree)
  }
  variable_importance <- rep(0, ncol(data_split))
  names(variable_importance) <- colnames(data_split)
  tree_losses <- .model_loss_drd_har(data)
  tree_splits <- ""
  candidate_splits <- history_tree <- tree_fits <- vector(mode = "list", length = 1)
  if (fit) tree_fits[[1]] <- .fit_tree(tree_splits = tree_splits, data = data, data_x = NULL)
  history_tree[[1]] <- list("tree_splits" = tree_splits, "tree_losses" = tree_losses, 
                            "tree_fits" = tree_fits, "variable_importance" = variable_importance)
  
  for (i in 1:miter) {
    best_node_to_split <- NULL
    best_loss_reduction <- 0
    for (j in seq_along(tree_splits)) {
      
      if (is.null(candidate_splits[[j]])) {
        node <- tree_splits[j]
        if (node == "") {
          idx_j <- data_x[, idx___]
        } else {
          idx_j <- data_x[eval(parse(text=node)), idx___]
        }
        split_j <- .get_split(data = data[idx___ %in% idx_j, ], data_x = data_x[idx___ %in% idx_j, ], data_split = data_split, 
                              n_ids = n_ids, min_obs = min_obs, min_ids = min_ids, mtry = mtry)
        candidate_splits[[j]] <- split_j
      } else {
        split_j <- candidate_splits[[j]]
      }
      
      loss_reduction_j <- tree_losses[j] - split_j$loss_left - split_j$loss_right
      
      if (loss_reduction_j>0 & loss_reduction_j>best_loss_reduction) {
        best_node_to_split <- j
        best_loss_reduction <- loss_reduction_j
      }
    }
    
    if (is.null(best_node_to_split)) break
    
    path_to_split <- tree_splits[best_node_to_split]
    if (path_to_split!="") {path_to_split <- paste(path_to_split, "&")}
    best_split <- candidate_splits[[best_node_to_split]]
    new_tree_split <- c(paste(path_to_split, 
                              best_split$variable, "<", best_split$value), 
                        paste(path_to_split, 
                              best_split$variable, ">=", best_split$value))
    new_tree_loss <- c(best_split$loss_left, best_split$loss_right)
    # add the reduction in loss to the variable importance of the variable used for the split
    variable_importance[best_split$variable] <- variable_importance[best_split$variable] + best_loss_reduction
    # update the lists defining the tree: paths to leafs, losses in each leaf, candidate splits of the leafs (we keep them so that we don't 
    # have to compute them at every iteration)
    candidate_splits <- append(candidate_splits, vector(mode = "list", length = 2), after = best_node_to_split)[-best_node_to_split]
    tree_splits <- append(tree_splits, new_tree_split, after = best_node_to_split)[-best_node_to_split]
    tree_losses <- append(tree_losses, new_tree_loss, after = best_node_to_split)[-best_node_to_split]
    # if requested, fit the leaf model and save the parameters
    if (fit) {
      new_tree_fit <- .fit_tree(new_tree_split, data = data, data_x = data_x)
      tree_fits <- append(tree_fits, new_tree_fit, after = best_node_to_split)[-best_node_to_split]
    }
    # save the current tree structure: we can easely access all steps of the tree-growing process
    history_tree[[i+1]] <- list("tree_splits" = tree_splits, "tree_loss" = sum(tree_losses), 
                                "tree_fits"=tree_fits, "variable_importance"=variable_importance)
  }
  
  return(history_tree)
}

# Get the best split for a dataset 
.get_split <- function(data, data_x, data_split, n_ids, min_obs, min_ids, mtry = NULL) {
  
  best_loss <- best_loss_left <- best_loss_right <- Inf
  best_index <- NULL
  best_value <- NULL
  best_variable <- NULL
  
  split_variables <- colnames(data_split)
  if (!is.null(mtry)) split_variables <- split_variables[sample.int(length(split_variables), mtry)]
  
  for (variable in split_variables) {
    values <- data_split[[variable]]
    # for discrete split variables, e.g. binary variables, we only check for the unique values
    values <- unique(values)
    for (value in values) {
      
      # Determine the split
      left_index <- which(data_x[,..variable]<value)
      right_index <- which(data_x[,..variable]>=value)
      if (length(left_index) == 0 | length(right_index)==0) {
        next
      } else {
        # Divide dataset
        left <- data[left_index,]
        right <- data[right_index, ]
        # Check the occurrences for each id in each split
        left_ids <- table(left[,id___])
        right_ids <- table(right[,id___])
        
        # Check if split is valid:
        if (length(left_ids) < n_ids |
            nrow(left) < min_obs |
            min(left_ids) < min_ids) next
        if (length(right_ids) < n_ids |
            nrow(right) < min_obs |
            min(right_ids) < min_ids) next
      }
      
      loss_left <- .model_loss_drd_har(left)
      loss_right <- .model_loss_drd_har(right)
      loss <- loss_left + loss_right
      
      
      
      if (loss < best_loss) {
        best_loss <- loss
        best_loss_left <- loss_left
        best_loss_right <- loss_right 
        best_variable <- variable
        best_value <- value
      }
    }
  }
  return(list("variable"=best_variable, 
              "value"=best_value,
              "loss_left"=best_loss_left,
              "loss_right"=best_loss_right))
}

# Fit the model to a leaf
.fit_tree <- function(tree_splits, data, data_x) {
  
  if (length(tree_splits) == 1) {
    fitted_tree <-  .fit_model_drd_har(data = data)
    return(fitted_tree)
  }
  
  fitted_tree <- vector(mode="list", length = length(tree_splits))
  for (i in seq_along(tree_splits)) {
    node <- tree_splits[i]
    idx_i <- data_x[eval(parse(text=node)), idx___]
    fitted_tree[[i]] <- .fit_model_drd_har(data = data[idx___ %in% idx_i, ])
  }
  
  return(fitted_tree)
}

# Compute the complexity cost of the tree based on information criteria (e.g. AIC, BIC)
.ic_cost_tree <- function(historical_tree, k=2, n_obs, n_param) {
  ic <-  k*length(historical_tree$tree_splits)*n_param + n_obs*log(historical_tree$tree_loss)
  return(ic)
}

# Compute the complexity cost of the tree based on a regression tree complexity function (see Hastie et. al. 2009, page 308)
.lambda_cost_tree <- function(historical_tree, lambda, n_obs, n_param) {
  complexity_cost <- historical_tree$tree_loss/n_obs + lambda*length(historical_tree$tree_splits)*n_param
  return(complexity_cost)
}