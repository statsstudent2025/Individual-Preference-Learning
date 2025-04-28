library(abind)
library(PLMIX)
library(gtools)
library(stats)

########### Mallows (Regular-No Clusters) sampling functions ##########################################################

random_permutation <-function(n){
  permutation <- sample(1:n)
  return(permutation)
}

validate_permutation <- function(rho) {
  return(is.numeric(rho) && all(sort(rho) == seq_along(rho)))
}

sample_mallows_custom <- function (rho0, alpha0, n_samples, leap_size = max(1L, floor(length(rho0)/5)), 
                                   metric = "footrule", diagnostic = FALSE, burnin = ifelse(diagnostic, 
                                                                                            0, 1000), thinning = ifelse(diagnostic, 1, 1000), items_to_plot = NULL, 
                                   max_lag = 1000L) 
{
  if (!(validate_permutation(rho0) && sum(is.na(rho0)) == 0)) {
    stop("rho0 must be a proper ranking with no missing values.")
  }
  if (diagnostic && n_samples == 1) {
    stop("Must have more than one sample to create diagnostic plots.")
  } else if (n_samples <= 0) {
    stop("n_samples must be positive.")
  }
  n_items <- length(rho0)
  
  if (diagnostic) {
    internal_burnin <- 0
    internal_thinning <- 1
    internal_n_samples <- burnin + n_samples * thinning
  } else {
    internal_burnin <- burnin
    internal_thinning <- thinning
    internal_n_samples <- n_samples
  }
  
  samples <- t(BayesMallows:::rmallows(rho0 = rho0, alpha0 = alpha0, n_samples = internal_n_samples, 
                                       burnin = internal_burnin, thinning = internal_thinning, 
                                       leap_size = leap_size, metric = metric))
  
  if (diagnostic) {
    # Fix: Only use the default selection when items_to_plot is NULL
    if (is.null(items_to_plot)) {
      if (n_items > 5) {
        message("Items not provided by user. Picking 5 at random.")
        items_to_plot <- sample.int(n_items, 5)
      } else {
        items_to_plot <- seq(from = 1, to = n_items, by = 1)
      }
    } 
    
    # Compute autocorrelation
    autocorr <- apply(samples[, items_to_plot, drop = FALSE], 
                      2, stats::acf, lag.max = max_lag, plot = FALSE, demean = TRUE)
    names(autocorr) <- items_to_plot
    autocorr <- do.call(rbind, Map(function(x, xnm) {
      data.frame(item = xnm, acf = x$acf[, 1, 1], lag = x$lag[, 1, 1])
    }, x = autocorr, xnm = names(autocorr)))
    autocorr$item <- as.factor(as.integer(autocorr$item))
    
    ac_plot <- ggplot2::ggplot(autocorr, ggplot2::aes(x = .data$lag, 
                                                      y = .data$acf, color = .data$item)) + ggplot2::geom_line() + 
      ggplot2::theme(legend.title = ggplot2::element_blank()) + 
      ggplot2::xlab("Lag") + ggplot2::ylab("Autocorrelation") + 
      ggplot2::ggtitle("Autocorrelation of Rank Values")
    print(ac_plot)
    
    con <- getOption("ask_opts.con", stdin())
    print("Press [enter] to see the next plot")
    response <- readLines(con = con, n = 1)
    
    colnames(samples) <- seq(from = 1, to = n_items, by = 1)
    diagnostic <- as.data.frame(samples)
    diagnostic$iteration <- seq_len(nrow(diagnostic))
    diagnostic <- stats::reshape(diagnostic, direction = "long", 
                                 varying = setdiff(names(diagnostic), "iteration"), 
                                 v.names = "value", timevar = "item", times = setdiff(names(diagnostic), 
                                                                                      "iteration"), idvar = "iteration", ids = diagnostic$iteration)
    
    # Fix: Ensure we only keep the specified items
    diagnostic <- diagnostic[diagnostic$item %in% items_to_plot, , drop = FALSE]
    diagnostic$item <- as.factor(as.integer(diagnostic$item))
    
    rho_plot <- ggplot2::ggplot(diagnostic, ggplot2::aes(x = .data$iteration, 
                                                         y = .data$value, color = .data$item)) + ggplot2::geom_line() + 
      ggplot2::theme(legend.title = ggplot2::element_blank()) + 
      ggplot2::xlab("Iteration") + ggplot2::ylab("Rank value") + 
      ggplot2::ggtitle("Trace Plot of Rank Values")
    print(rho_plot)
    
    samples <- samples[seq(from = burnin + 1, by = thinning, 
                           length.out = n_samples), ]
  }
  
  return(samples)
}

########## Plackett-Luce Sampling function-Updated ( Now includes latent cluster labels)#######################################
gibbsPLMIX_updated<-function (pi_inv, K, G, init = list(z = NULL, p = NULL), n_iter = 1000, 
                              n_burn = 500, hyper = list(shape0 = matrix(1, nrow = G, ncol = K), 
                                                         rate0 = rep(0.001, G), alpha0 = rep(1, G)), centered_start = FALSE) 
{
  cl = match.call()
  if (class(pi_inv)[1] != "top_ordering") {
    if (class(pi_inv)[1] == "RankData") {
      pi_inv = as.top_ordering(data = pi_inv)
    }
    if (class(pi_inv)[1] == "rankings") {
      pi_inv = as.top_ordering(data = pi_inv)
    }
    if (class(pi_inv)[1] == "matrix" | class(pi_inv)[1] == 
        "data.frame") {
      pi_inv = as.top_ordering(data = pi_inv, format_input = "ordering", 
                               aggr = FALSE)
    }
  }
  pi_inv <- PLMIX:::fill_single_entries(data = pi_inv)
  N <- nrow(pi_inv)
  n_rank <- PLMIX:::howmanyranked(pi_inv)
  rho <- matrix(1:K, nrow = G, ncol = K, byrow = TRUE)
  if (is.null(init$z)) {
    z <- binary_group_ind(class = sample(1:G, size = N, replace = TRUE), 
                          G = G)
  }
  else {
    z <- init$z
  }
  omega <- colMeans(z)
  if (is.null(init$p)) {
    if (centered_start) {
      print("CENTERED START !!")
      mle1comp <- matrix(prop.table(table(factor(pi_inv[, 
                                                        1], levels = 1:K))), nrow = 1)
      p <- random_start(mlesupp = mle1comp, givenweights = omega)
    }
    else {
      print("COMPLETELY RANDOM (uniform support, rescaled) START")
      p <- matrix(rgamma(n = G * K, shape = 1, rate = 1), 
                  nrow = G, ncol = K)
    }
  }
  else {
    p <- init$p
  }
  shape0 <- hyper$shape0
  rate0 <- hyper$rate0
  alpha0 <- hyper$alpha0
  u_bin <- PLMIX:::umat(pi_inv = pi_inv)
  log_lik <- c(loglikPLMIX(p = p, ref_order = rho, weights = omega, 
                           pi_inv = pi_inv), rep(NA, n_iter))
  log_prior <- c(log(ddirichlet(omega, alpha0)) + sum(dgamma(p, 
                                                             shape = shape0, rate = rate0, log = TRUE)), rep(NA, n_iter))
  objective <- log_lik + log_prior
  Pi <- array(NA, dim = c(G, K, n_iter + 1)) #for the utility vectors
  Pi[, , 1] <- p
  Zeta <- z
  Zeta_matrix <-  matrix(NA, nrow = n_iter + 1, ncol = N) # FOr the latent vectors
  Zeta_matrix[1,] <-max.col(z) 
  
  Omega <- matrix(NA, nrow = n_iter + 1, ncol = G) # for the weights
  Omega[1, ] <- omega
  for (l in 1:n_iter) {
    if (l%%500 == 0) {
      print(paste("GIBBS iteration", l))
    }
    Omega[l + 1, ] <- rdirichlet(n = 1, alpha = alpha0 + 
                                   colSums(Zeta)) #weights updated
    temprate <- PLMIX:::CompRateYpartial(p = adrop(Pi[, , l, drop = FALSE], 
                                           3), pi_inv = pi_inv, ref_order = rho, z = Zeta, n_rank = n_rank)
    Ypsilon <- PLMIX:::SimYpsilon(rate = temprate, n_rank = n_rank)
    Pi[, , l + 1] <- matrix(rgamma(n = G * K, shape = shape0 + 
                                     PLMIX:::gammamat(u_bin = u_bin, z_hat = Zeta), rate <- PLMIX:::CompRateP(pi_inv = pi_inv, 
                                                                                              Y = Ypsilon, z = Zeta, u_bin = u_bin, n_rank = n_rank, 
                                                                                              rate0 = rate0)), nrow = G, ncol = K)
    Zeta <- binary_group_ind(apply(PLMIX:::CompProbZpartial(p = adrop(Pi[, 
                                                                 , l + 1, drop = FALSE], 3), pi_inv = pi_inv, Y = Ypsilon, 
                                                    u_bin = u_bin, n_rank, omega = Omega[l + 1, ]), 1, 
                                   FUN = sample, x = 1:G, replace = TRUE, size = 1), 
                             G = G)
    Zeta_matrix[l+1,] <- max.col(Zeta) 
    log_lik[l + 1] <- loglikPLMIX(p = adrop(Pi[, , l + 1, 
                                               drop = FALSE], 3), ref_order = rho, weights = Omega[l + 
                                                                                                     1, ], pi_inv = pi_inv)
    log_prior[l + 1] <- log(ddirichlet(Omega[l + 1, ], alpha0)) + 
      sum(dgamma(adrop(Pi[, , l + 1, drop = FALSE], 3), 
                 shape = shape0, rate = rate0, log = TRUE))
    objective[l + 1] <- log_lik[l + 1] + log_prior[l + 1]
  }
  log_lik <- log_lik[-c(1:(n_burn + 1))]
  objective <- objective[-c(1:(n_burn + 1))]
  Omega <- Omega[-c(1:(n_burn + 1)), , drop = FALSE]
  colnames(Omega) <- paste0("w_", 1:G)
  
  Zeta_matrix <- Zeta_matrix[-c(1:(n_burn + 1)), , drop = FALSE]
  colnames(Zeta_matrix) <- paste0('z_', 1:ncol(Zeta_matrix))
  Pi <- array(apply(Pi, 3, FUN = function(x) x/rowSums(x)), 
              c(G, K, n_iter + 1))
  Pi = t(apply(Pi, 3, c))[-c(1:(n_burn + 1)), ]
  colnames(Pi) <- paste0("p_", rep(1:G, K), rep(1:K, each = G))
  out = list(W = Omega, P = Pi, log_lik = log_lik, deviance = -2 * 
               log_lik, Z = Zeta_matrix,objective = objective, call = cl)
  class(out) = "gsPLMIX"
  return(out)
}
############## Partial Data Function########################################################################

making_partial_new <- function(dataset,lambda_K = ncol(dataset)/2,Kmax = ncol(dataset),Poisson_way = TRUE){
  # Step 4: Mask unranked items with 0
  N = nrow(dataset)
  partial_rankings <- matrix(0, nrow = N, ncol = ncol(dataset) ) # Initialize with zeros
  if (Poisson_way){
    K_j <- rtpois(N, lambda = lambda_K, a = 1, b = Kmax)  # K_j values between 1 and n
  }
  
  else{
    K_j <- rep(lambda_K,N)
  }
  for (j in 1:N) {
    top_K <- K_j[j]  # Get the number of revealed ranks for this assessor
    
    # Find indices of the top-K ranked items
    top_K_indices <- which(dataset[j, ] %in% 1:top_K)  
    # print(top_K_indices)
    
    partial_rankings[j, top_K_indices] <- as.numeric(dataset[j,top_K_indices])  # Keep top K, mask others
    
    Mallows_partial <- data.frame(partial_rankings)
    Plackett_partial <- as.top_ordering((Mallows_partial),format_input = 'ranking', aggr = FALSE)
    
  }
  
  return(list(Mallows_partial,Plackett_partial,K_j))
  
  
  
}
############# Plackett-Luce Functions for computing probabilities################################



pl_top_m_probs_exact <- function(theta, m) {
  n <- length(theta)  # Number of items
  probs <- numeric(n)  # Vector to store probabilities
  top_m_sets <- permutations(n, m, v = 1:n)  # All orderings of top-m choices
  
  # Compute probability for each item
  for (i in 1:n) {
    
    prob_sum <- 0  # To store the total probability for item i in top-m
    
    # Iterate over all top-m rankings where item i appears
    for (row in 1:nrow(top_m_sets)) {
      top_m_ranking <- top_m_sets[row, ]
      
      if (i %in% top_m_ranking) {
        # Compute Plackett-Luce probability of this specific ranking
        prob <- 1
        remaining_theta <- sum(theta)  # Total sum of theta
        
        for (j in 1:m) {
          prob <- prob * (theta[top_m_ranking[j]] / remaining_theta)
          remaining_theta <- remaining_theta - theta[top_m_ranking[j]]
        }
        
        prob_sum <- prob_sum + prob  # Sum over all valid rankings
      }
    }
    
    # Store the exact probability
    probs[i] <- prob_sum
  }
  
  return(probs)
}

# Function to compute probabilities row-wise considering ranked and unranked items
# Computes probability of an unranked item being in the 'next top m'

compute_probabilities <- function(df, theta_vec, m = 2) {
  n <- ncol(df)  # Number of columns (also range 1 to n)
  result_probs <- matrix(0, nrow = nrow(df), ncol = n)  # Store probabilities
  
  for (row_idx in 1:nrow(df)) {
    row <- df[row_idx, ]  # Get current row
    ranked_items <- row[row != 0]  # Items already ranked
    unranked_items <- setdiff(1:n, ranked_items)  # Items missing in this row
    
    # Initialize probability vector (ranked items get 1)
    row_probs <- rep(0, n)
    row_probs[ranked_items] <- 1  # Ranked items assigned probability = 1
    
    # Compute probabilities only for unranked items
    if (length(unranked_items) > m) {
      theta_unranked <- theta_vec[unranked_items]  # Get theta values for unranked
      top_m_probs <- pl_top_m_probs_exact(theta_unranked, m)  # Compute exact probs
      
      # Assign computed probabilities to unranked items
      row_probs[unranked_items] <- top_m_probs
    }
    else {
      
      row_probs[unranked_items]<- 1
    }
    
    # Store results
    result_probs[row_idx, ] <- row_probs
  }
  
  return(result_probs)
}

###### Functions for rank prediction using Bayes-Risk Methodology , taking footrule distance as Loss######
compute_footrule <- function(candidate, reference) {
  pos1 <- setNames(seq_along(candidate), candidate)
  pos2 <- setNames(seq_along(reference), reference)
  common <- intersect(names(pos1), names(pos2))
  sum(abs(pos1[common] - pos2[common]))
}

compute_bayes_risk <- function(candidate_order, posterior_samples) {
  mean(sapply(posterior_samples, function(sample_order) {
    compute_footrule(candidate_order, sample_order)
  }))
}

generate_neighbors <- function(perm) {
  neighbors <- list()
  for (i in 1:(length(perm) - 1)) {
    swapped <- perm
    swapped[c(i, i+1)] <- swapped[c(i+1, i)]
    neighbors[[length(neighbors) + 1]] <- swapped
  }
  neighbors
}

greedy_bayes_optimizer <- function(posterior_samples, tol = 0.1, max_iter = 100) {
  current <- posterior_samples[[1]]  # Start with ranking from first iteration
  current_risk <- compute_bayes_risk(current, posterior_samples)
  
  risk_trace <- data.frame(iteration = 0, bayes_risk = current_risk)
  
  for (i in 1:max_iter) {
    neighbors <- generate_neighbors(current)
    risks <- map_dbl(neighbors, compute_bayes_risk, posterior_samples = posterior_samples)
    
    best_risk <- min(risks)
    best_neighbor <- neighbors[[which.min(risks)]]
    
    risk_trace <- rbind(risk_trace, data.frame(iteration = i, bayes_risk = best_risk))
    
    if ((current_risk - best_risk) < tol) {
      message("Stopping: no significant improvement.")
      break
    }
    
    current <- best_neighbor
    current_risk <- best_risk
    
  }
  
  return(list(optimal_ranking = current, bayes_risk = current_risk,  risk_trace = risk_trace))
}

#################### Plackett-Luce Sampling Rankings ########################

SamplingCompleteRankings_PL <- function(score_params){
  
  # Simulate exponential variables for each rate in score_function'
  simulated_values<- sapply(score_params,function(rate) rexp(1,rate))
  
  #Now I want to get the complete rankings for Plackett-Luce Data
  
  rankings <- rank(simulated_values)
  
  
  
  return(rankings)
}

MonteCarloPosteriorProbs <- function (rankings,item_names) {
  P <- matrix(0,nrow = length(item_names), ncol = length(item_names))
  n<- nrow(rankings)
  for (i in 1:nrow(rankings)){
    
    curr_row <- as.vector(rankings[i,])
    K<- length(item_names)
    
    # Ranking is a 1xK vector
    
    for(item in 1:K){
      ranking_place <- curr_row[item]
      P[ranking_place,item] <- P[ranking_place,item] + 1
    }
    
    
  }
  # Now we find the probabilities just by dividing by the amount of rows of data (empirically)- Monte Carlo
  P<- P/n
  
  colnames(P) <- item_names
  P_df <- as.data.frame(P)
  #Now we have the posterior probabilities in a dataframe
  return(P_df)
  
}

# Function to simulate N rankings
SimulateNRankings <- function(score_params, N) {
  simulations <- replicate(N, SamplingCompleteRankings_PL(score_params))
  df <- as.data.frame(t(simulations))
  colnames(df) <- names(score_params)  # name columns if items are named
  rownames(df) <- paste0("Sim_", 1:N)
  return(df)}


######### Finding Probability of next m (UPDATED the ones above take too long) ##############################################
# Function to compute exact probability of each item in the top m (Plackett-Luce)
Next_m_probs<- function(m,data_){
  start_time <- Sys.time() #Start time
  T<- nrow(scores_PL)
  N <- nrow(data_) # Number of assessors
  K <- ncol(data_)
  posterior_nextm_PL <- (matrix(0, nrow = N, ncol = K) ) # N = assessors, K =  no of items
  for (i in 1:T){
    
    utility_i <- scores_PL[i,] #ith iteration of utility scores
    for (j in 1:nrow(data_)){ #iterating by assessor
      curr_assessor <- PL_partial[j,]
      ranked_items <- curr_assessor[curr_assessor!= 0]  # Items already ranked
      unranked_items <- setdiff(1:K, ranked_items)  # Items missing in this row
      
      posterior_nextm_PL[j,ranked_items] = posterior_nextm_PL[j,ranked_items] + 1
      
      # Compute probabilities only for unranked items
      if (length(unranked_items) > m) {
        theta_unranked <- utility_i[unranked_items]  # Get theta values for unranked
        simulated_ranking <- SamplingCompleteRankings_PL(theta_unranked)  # Compute exact probs
        
        # Then add + 1 only to the places that got a rank <=m
        positions_to_update <- which(simulated_ranking <=m)
        posterior_nextm_PL[j,unranked_items[positions_to_update]]=  posterior_nextm_PL[j,unranked_items[positions_to_update]] + 1
      } 
      
    }
    
    
  }
  posterior_nextm_PL =  posterior_nextm_PL/T
  end_time<- Sys.time()
  # Calculate the difference
  execution_time <- end_time - start_time
  print(paste("Time taken: ", execution_time))
  
  return (posterior_nextm_PL)
}

#For Mallows top-m probabilities


Next_m_probs_MM<- function(m){
  start_time <- Sys.time() #Start time
  df_try <-result_BMM$augmented_data[result_BMM$augmented_data$iteration > burnin(result_BMM),]
  
  number_of_samples_per_assesor <- rowSums(PL_partial!=0)
  
  assessors <- unique(df_try$assessor)
  items <- unique(df_try$item)
  
  posterior_prob_BMM <- matrix(0,nrow = length(assessors), ncol = length(items), dimnames = list(assessors,items))
  
  n_samples <- length(unique(df_try$iteration))
  
  for (i in seq_along(assessors)) {
    k <- number_of_samples_per_assesor[i]
    for (j in seq_along(items)) {
      count <- sum(df_try$assessor == assessors[i] & df_try$item == items[j] & df_try$value <= k+m ) # Item in top k + m
      posterior_prob_BMM[i, j] <- count / n_samples
    }
  }
  end_time<- Sys.time()
  # Calculate the difference
  execution_time <- end_time - start_time
  print(paste("Time taken: ", execution_time))
  return(list(posterior_prob_BMM,assessors,items))
  
  
  
}
######################## Plotting Functions/Other functions##########################
library(scales)
plot_calibration_curve <- function(posterior_matrix, model_name, m,data_) {
  # Ensure posterior is in long form
  posterior_df <- as.data.frame(as.table(posterior_matrix))
  colnames(posterior_df) <- c("Assessor", "Item", "Probability")
  
  # Construct ranking reference
  total_true_rank <- as.data.frame(as.table(data_))
  colnames(total_true_rank) <- c("Assessor", "Item", "True_Rank")
  
  # Merge with ranking and filter for unranked
  df <- total_true_rank %>%
    left_join(Partial_Long, by = c('Assessor', 'Item')) %>%
    mutate(unranked = Rank == 0)
  
  df_merged <- df %>%
    left_join(posterior_df, by = c("Assessor", "Item")) %>%
    filter(unranked == TRUE) %>%
    filter(samples_per_assesor[as.integer(Assessor)] < ncol(data_)-2) %>%
    mutate(Success = samples_per_assesor[as.integer(Assessor)] + m >= True_Rank)
  
  # Bin probabilities and compute empirical success/failure rate
  df_binned <- df_merged %>%
    mutate(Prob_Bin = cut(Probability, breaks = seq(0, 1, by = 0.25), include.highest = TRUE)) %>%
    group_by(Prob_Bin, Success) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Prob_Bin) %>%
    mutate(Total = sum(Count),
           Empirical_Prob = Count / Total)
  
  # Plot calibration
  ggplot(df_binned, aes(x = Prob_Bin, y = Empirical_Prob, fill = Success)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_text(aes(label = paste0(scales::percent(Empirical_Prob, accuracy = 0.1), "\n(n = ", Total, ")")),
              position = position_dodge(width = 0.9),
              vjust = -0.3, size = 3)+
    scale_fill_manual(values = c("red", "green"), labels = c("Fail", "Success")) +
    ylab("Empirical Proportion") +
    xlab("Predicted Probability Bin") +
    ggtitle(paste("Empirical Success/Fail by Predicted Probability Bin\nModel:", model_name)) +
    expand_limits(y = 1.1 * max(df_binned$Empirical_Prob, na.rm = TRUE)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    )
}

# Convert to long format- this is for mallows
#Partial_Long <- MM_partial %>%
 # mutate(Assessor = factor(row_number())) %>%
 # pivot_longer(cols = starts_with('X'),
         #      names_to = 'Item',
          #     values_to = 'Rank')

library(dplyr)
library(ggplot2)
library(tidytext)  # for reorder_within() and scale_x_reordered()

plot_top_m_unranked_probs <- function(prob_data, model_name, top_m = 2) {
  data <- prob_data %>%
    left_join(Partial_Long, by = c('Assessor', 'Item')) %>%
    left_join(true_rank, by = c('Assessor', 'Item')) %>%
    mutate(unranked = Rank == 0) %>%
    filter(unranked == TRUE ) %>%
    mutate(Assessor = factor(Assessor, labels = paste("Individual", 1:length(unique(Assessor)))))
  
  if (nrow(data) == 0) {
    warning(paste("No unranked items with True Rank <=", top_m, "found for model:", model_name))
    return(NULL)
  }
  
  data %>%
    group_by(Assessor) %>%
    mutate(Item = reorder_within(Item, True_Rank, Assessor)) %>%
    ggplot(aes(x = Item, y = Probability, fill = Probability)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~ Assessor, scales = 'free_x', ncol = 3) +
    scale_x_reordered() +
    scale_fill_gradient(low = "blue", high = "red", name = "Probability") +
    labs(
      title = paste('Posterior Probabilities for Top', top_m, 'Unranked Items -', model_name),
      x = 'Item (sorted by True Rank)',
      y = 'Probability'
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "black", color = "black"),
      strip.text.x = element_text(color = "white"),
      panel.spacing = unit(1, "lines")
    )
}



