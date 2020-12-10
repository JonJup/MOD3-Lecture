## -- RCode3 Functions 

# Code S3.R
# R code function for estimating Gaussian copula correlation matrix parameters
# via Monte Carlo Expectation Maximisation (MCEM).
#
# Example Usage of the code can be found at the end of the file.
#


library(stats)
library(utils)
library(graphics)
library(MASS)

library(ZIM)
library(extraDistr)
library(mnormt)
library(dfoptim)

###################################################
# The following functions form the user interface #
###################################################

#' Find the approximate MLE for the correlation matrix of a gaussian copula with discrete marginals.
#'
#' @param data A matrix of non-negative integer valued observations. Each row is
#'   a single sampling unit and each column relates to a single species. Each column
#'   must be named and the name should match one of the entries in marginal_details.
#'
#' @param marginal_details A data.frame detailing the form and parameters of the
#'   marginal distribution for each species.
#'
#'   The data.frame must contain the following five named columns: "Species",
#'   "Distr", "mu", "theta" and "pi". Entries in the "Species" column must each
#'   match a species in colnames(data). "Distr" indicates what type of
#'   marginal is assigned to that species; currently supported options are:
#'   "POIS" (Poisson), "NB" (negative binomial), "ZIP" (zero-inflated Poisson)
#'   and "ZINB" (zero-inflated negative binomial). Entries in the "mu", "theta"
#'   and "pi" columns indicate the relevant parameter for the marginal
#'   distribution. Any entries not required for a given type of marginal should
#'   be set to "NA"
#'
#' @param max_iterations The maximum number of MCEM iterations that should be
#'   performed.
#'
#' @param mcem_samples The number of Monte-carlo samples that should be used
#'   when estimate the distribution of each latent variable. Higher values here
#'   will produce improved estimates, but will increase the computational time
#'   required by the algorithm.
#'
#' @param verbosity Non-negative integer. Higher values will output more
#'   information whilst the function is running. A value of 0 here will surpress
#'   all non-critical information from being output whilst the function is
#'   running.
#'
#' @param ML_estimator_func A function which returns the MLE for the covariance
#'   matrix of a MVN model given a set of weighted data. Changing this from the
#'   default value allows the estimation of other MVN based copula models
#'   besides the standard gaussian copula. The function, f, should have the
#'   signature: f(data, weights, optim_params) where optim_params is a list
#'   containing hints which may be used to optimise the function's performance.
#'
#' @param initial_cov The initial covariance matrix that should be used as the
#'   starting point of the MCEM process.
#'
#' @param fixed_latent_data A boolean which determines whether the monte-carlo
#'   samples representing the latent data should be changed between each
#'   iteration. Setting this value to TRUE will speed up the function, but may
#'   lead to over fitting.
#'
#' @param latent_data A matrix containing predetermined latent data. The number
#'   of columns should equal the number of species listed in marginal_details.
#'   The number of rows should be a multiple of the number of rows in data. If
#'   this parameter is not NULL, then mcem_samples must be NULL and
#'   fixed_jitters must be TRUE.
#'
#' @return A named list containing details of the MCEM execution. This list includes:
#'         "convergence" - A non-negative integer. A non-zero value indicates that the MCEM algorithm did NOT converge.
#'         "cov_final" - The best MLE of the copula's correlation matrix found by the function.
#'         "cov_final_log_lik" - An estimate of the log likelihood of "cov_final"
#'
#' @examples
#' data <- matrix(c(14, 11, 17, 6, 9, 9, 16, 11, 13, 7, 26, 14, 49, 1,
#'                  1, 2, 40, 11, 20, 2, 24, 22, 29, 0, 0, 0, 26, 20, 23, 0, 46,
#'                  31, 92, 0, 0, 0, 51, 23, 30, 0),
#'                nrow = 10,
#'                ncol = 4,
#'                dimnames = list(NULL, c("species_1",  "species_2", "species_3", "species_4")))
#' marginal_details <- data.frame("Species" = c("species_1", "species_2", "species_3", "species_4"),
#'                                "Distr" = c("POIS", "NB", "ZIP", "ZINB"),
#'                                "mu" = c(10, 15, 20, 30),
#'                                "theta" = c(NA, 1, NA, 2),
#'                                "pi" = c(NA, NA, 0.2, 0.3))
#' mcem_result <- estimate_copula(data, marginal_details)
#'
#' summary(mcem_result)
#'
#' @export
estimate_copula <- function(data,
                            marginal_details,
                            max_iterations = 100,
                            mcem_samples = 10000,
                            verbosity = 1,
                            ML_estimator_func = MLE_correlation_matrix,
                            initial_cov = NULL,
                            fixed_latent_data = FALSE,
                            latent_data = NULL) {
        data <- as.matrix(data)
        
        N <- nrow(data)
        p <- ncol(data)
        
        # Check that species listed in marginal_details matches thos present in data.
        if (sum(colnames(data) %in% marginal_details[["Species"]]) != p) {
                stop("Not all species listed in the data are present in the marginals_details")
        }
        
        # Find the order that the species in data are listed in marginal_details.
        species_indices <- match(colnames(data), marginal_details[["Species"]])
        
        marginals <- as.list(rep(NA, p))
        for (i in 1:p) {
                # Potentially reorder the marginals, so that the order of the species matches that in data.
                species_index <- species_indices[[i]]
                
                # Convert the marginal details into the internal format used by the mcem_copula function.
                marginals[[i]] <- generate_marginal(marginal_details[species_index, "Distr"],
                                                    marginal_details[species_index, c("mu", "theta", "pi")])
        }
        
        # Call the function which actually finds the MLE
        mcem_result <- mcem_copula(data,
                                   marginals,
                                   max_iterations = max_iterations,
                                   mcem_samples = mcem_samples,
                                   initial_cov = initial_cov,
                                   unobserved_data = latent_data,
                                   fixed_jitters = fixed_latent_data,
                                   verbosity = verbosity,
                                   cov_matrix_function = ML_estimator_func)
        return(mcem_result)
}

#' Estimate the log likelihood of gaussian copula given a set of discrete observations.
#'
#' @param data A matrix of non-negative integer valued observations. Each row is
#'   a single sampling unit and each column relates to a single species. Each column
#'   should be named and the name should match one of the entries in marginal_details.
#'
#' @param marginal_details A data.frame detailing the form and parameters of the
#'   marginal distribution for each species. The data.frame must contain the
#'   following five named columns: "Species", "Distr", "mu", "theta" and "pi".
#'   Entries in the "Species" column must each match a species listed in
#'   colnames(data). "Distr" indicates what type of marginal is assigned to
#'   that species; currently supported options are: "POIS" (Poisson), "NB"
#'   (negative binomial), "ZIP" (zero-inflated Poisson) and "ZINB"
#'   (zero-inflated negative binomial). Entries in the "mu", "theta" and "pi"
#'   columns indicate the relevant parameter for the marginal distribution. Any
#'   entries not required for a given type of marginal should be set to "NA"
#'
#' @param cov A matrix which species the covariance structure of the copula. For
#'   a gaussian copula, this will be it's correlation matrix. This should be a
#'   square matrix with dimension equal to the number of marginals provided in
#'   marginal_details.
#'
#' @param mc_samples The number of monte-carlo samples to take when estimating
#'   the distribution of the latent data. Increasing this value will improve the
#'   accuracy of the result at the expense of increased computational times.
#'
#' @param latent_data A matrix containing predetermined latent data. The number
#'   of columns should equal the number of species listed in marginal_details.
#'   The number of rows should be a multiple of the number of rows in data. If
#'   this parameter is not NULL, then mc_samples must be NULL.
#'
#' @return The log likelihood as a numeric value.
#'
#' @examples
#' data <- matrix(c(14, 11, 17, 6, 9, 9, 16, 11, 13, 7, 26, 14, 49, 1,
#'                  1, 2, 40, 11, 20, 2, 24, 22, 29, 0, 0, 0, 26, 20, 23, 0, 46,
#'                  31, 92, 0, 0, 0, 51, 23, 30, 0),
#'                nrow = 10,
#'                ncol = 4,
#'                dimnames = list(NULL, c("species_1",  "species_2", "species_3", "species_4")))
#' marginal_details <- data.frame("Species" = c("species_1", "species_2", "species_3", "species_4"),
#'                                "Distr" = c("POIS", "NB", "ZIP", "ZINB"),
#'                                "mu" = c(10, 15, 20, 30),
#'                                "theta" = c(NA, 1, NA, 2),
#'                                "pi" = c(NA, NA, 0.2, 0.3))
#' corr_guess <- cor(data)
#'
#' print("Estimate of the log Likelihood of corr_guess given data and marginal_details:")
#' estimate_log_likelihood(data, marginal_details, corr_guess)
#'
#' @export
estimate_log_likelihood <- function(data, marginal_details, cov, mc_samples = 10000, latent_data = NULL) {
        data <- as.matrix(data)
        
        N <- nrow(data)
        p <- ncol(data)
        
        # Check that species listed in marginal_details matches those present in data.
        if (sum(colnames(data) %in% marginal_details[["Species"]]) != p) {
                stop("Not all species listed in the data are present in the marginals_details")
        }
        
        # Find the order that the species in data are listed in marginal_details.
        species_indices <- match(colnames(data), marginal_details[["Species"]])
        
        marginals <- as.list(rep(NA, p))
        for (i in 1:p) {
                # Potentially reorder the marginals, so that the order of the species matches that in data.
                species_index <- species_indices[[i]]
                # Convert the marginal details into the internal format used by the mcem_copula function.
                marginals[[i]] <- generate_marginal(marginal_details[species_index, "Distr"],
                                                    marginal_details[species_index, c("mu", "theta", "pi")])
        }
        
        # Call the function which will find the log likelihood.
        if (! is.null(latent_data)) {
                llik <- log_likelihood_from_latent(data, marginals, cov, latent_data = latent_data)
        }
        else {
                llik <- log_likelihood(data, marginals, cov, mc_samples)
        }
        
        return(llik)
}

#' Generate a set of simulated discrete observations from a gaussian copula.
#'
#' @param N The number of observations to generate.
#'
#' @param marginal_details A data.frame detailing the form and parameters of the
#'   marginal distribution for each species. The data.frame must contain the
#'   following five named columns: "Species", "Distr", "mu", "theta" and "pi".
#'   Entries in the "Species" column must each match a species listed in
#'   colnames(data). "Distr" indicates what type of marginal is assigned to that
#'   species; currently supported options are: "POIS" (Poisson), "NB" (negative
#'   binomial), "ZIP" (zero-inflated Poisson) and "ZINB" (zero-inflated negative
#'   binomial). Entries in the "mu", "theta" and "pi" columns indicate the
#'   relevant parameter for the marginal distribution. Any entries not required
#'   for a given type of marginal should be set to "NA"
#'
#' @param cov The correlation matrix for the gaussian copula. A square matrix
#'   whose dimension should match the number of entries in marginal_details.
#'
#' @return A list, x, with three entries:
#'         x[["observed"]] - the simulated discrete valued observations;
#'         x[["censored"]] - the latent fractional components of each observation; and
#'         x[["cov_actual"]] - the covariance matrix of the full data (including the latent component) in the MVN space.
#'
#' @examples
#' marginal_details <- data.frame("Species" = c("species_1", "species_2", "species_3", "species_4"),
#'                                "Distr" = c("POIS", "NB", "ZIP", "ZINB"),
#'                                "mu" = c(10, 15, 20, 30),
#'                                "theta" = c(NA, 1, NA, 2),
#'                                "pi" = c(NA, NA, 0.2, 0.3))
#'
#' corr <- rcorrelation_matrix(1, 4)[ , , 1]
#' rownames(corr) <- colnames(corr) <- marginal_details[["Species"]]
#'
#' simulated_data <- generate_copula_data(10, marginal_details, corr)
#'
#' print("Simulated observed data:")
#' print(simulated_data[["observed"]])
#'
#' @export
generate_copula_data <- function(N, marginal_details, cov) {
        p <- nrow(marginal_details)
        
        if (! all(colnames(cov) == rownames(cov))) {
                stop("cov must be a square matrix with matching row and column names")
        }
        
        # Check that species listed in marginal_details matches those present in data.
        if (sum(colnames(cov) %in% marginal_details[["Species"]]) != p) {
                stop("Not all species listed in colnames(cov) are present in the marginals_details")
        }
        
        # Find the order that the species in data are listed in marginal_details.
        species_indices <- match(colnames(cov), marginal_details[["Species"]])
        
        marginals <- as.list(rep(NA, p))
        for (i in 1:p) {
                # Potentially reorder the marginals, so that the order of the species matches that in data.
                species_index <- species_indices[[i]]
                # Convert the marginal details into the internal format used by the mcem_copula function.
                marginals[[i]] <- generate_marginal(marginal_details[species_index, "Distr"],
                                                    marginal_details[species_index, c("mu", "theta", "pi")])
        }
        
        # Call the function which generates the data.
        data <- generate_data(N, marginals, cov)
        
        # Attach the species names from marginal_details to the output.
        colnames(data[["observed"]]) <- colnames(cov)
        colnames(data[["censored"]]) <- colnames(cov)
        
        return(data)
}

#' Find the MLE of the correlation matrix of a MVN model with unit variances, given a set of weighted data.
#'
#' @param data The set of observations.
#'
#' @param weights The weighting of each observations. The number of entries
#'   should match the number of rows in data.
#'
#' @param optim_params A named list, which may contain an entry "initial_cov".
#'   If it contains this entry, then it is used as a starting point for the
#'   numerical optimisation procedure.
#'
#' @return The MLE for the correlation matrix of a MVN model with with unit variances given the weighted data.
#'
#' @export
MLE_correlation_matrix <- function(data, weights = NULL, optim_params = list()) {
        N <- nrow(data)
        p <- ncol(data)
        
        if (is.null(weights)) {
                weights <- rep(1 / N, N)
        }
        else {
                weights <- weights / sum(weights)
        }
        
        if (length(weights) != N) {
                stop("Length of weights is not compatible with the dimension of data.")
        }
        
        # It is challenging to find this MLE in closed form.
        # Instead, we will using a numerical optimiser.
        
        # The MLE can be found by optimising against the sample covariance matrix,
        # so use this instead of the whole data set.
        condensed_data <- mean_covariance_matrix(data, weights)
        
        # Call a function which will perform the numerical optimisation.
        result <-  corr_ML_condensed_data(condensed_data, N, optim_params = optim_params)
        
        return(result)
}

#' Find the MLE of the covariance matrix of a MVN model given a set of weighted data.
#'
#' @param data The set of observations.
#'
#' @param weights The weighting of each observations. The number of entries
#'   should match the number of rows in data.
#'
#' @param optim_params Not used.
#'
#' @return The MLE for the covariance matrix of a MVN model given the weighted data.
#'
#' @export
MLE_covariance_matrix <- function(data, weights = NULL, optim_params = NULL) {
        N <- nrow(data)
        p <- ncol(data)
        
        if (is.null(weights)) {
                weights <- rep(1 / N, N)
        }
        else {
                weights <- weights / sum(weights)
        }
        
        if (length(weights) != N) {
                stop("Length of weights is not compatible with the dimension of data.")
        }
        
        result <- crossprod(data, sweep(data, MARGIN = 1, STATS = weights, FUN = '*'))
        
        # Ensure the result is completely symmetric
        result <- 1 / 2 * (result + t(result))
        
        return(result)
}

#' Print a brief outline of the result of an MCEM procedure
#'
#' @param x The value returned by a call to estimate_copula
#' @param ... Not used
#'
#' @return x
#'
#' @export
print.copuladiscreter_result <- function(x, ...) {
        cat("\nApproximate MLE for the gaussian copula covariance matrix ($cov_final):\n\n")
        print(x[["cov_final"]])
        cat("\n\n")
        
        cat(paste0("Log likelihood of the result ($cov_final_log_lik): ", x[["cov_final_log_lik"]], "\n\n"))
        
        cat(paste0("Convergence ($convergence): ", x[["convergence"]], "\n"))
        
        invisible(x)
}

#' Print a summary of the result of an MCEM procedure
#'
#' @param object The value returned by a call to estimate_copula
#' @param ... Not used
#'
#' @return object
#'
#' @export
summary.copuladiscreter_result <- function(object, ...) {
        cat("\nMLE for the gaussian copula covariance matrix ($cov_final):\n\n")
        print(object[["cov_final"]])
        cat("\n")
        
        cat("Variance of the MLE for the gaussian copula covariance matrix ($cov_final_var):\n\n")
        print(object[["cov_final_var"]])
        cat("\n")
        
        cat(paste0("Log likelihood of the result ($cov_final_log_lik): ", object[["cov_final_log_lik"]], "\n\n"))
        
        cat(paste0("Convergence ($convergence): ", object[["convergence"]], "\n\n"))
        
        cat(paste0("Complete list of fields:\n"))
        print(names(object))
        
        graphics::plot(unlist(object[["cov_steps_log_lik"]]), ylab = "log Likelihood", xlab = "MCEM iteration")
        
        invisible(object)
}

###############################################################
# This is the core function which executes the MCEM algorithm #
###############################################################
mcem_copula <- function(data,
                        marginals,
                        mcem_samples = NULL,
                        initial_cov = NULL,
                        max_iterations = 100,
                        unobserved_data = NULL,
                        fixed_jitters = FALSE,
                        verbosity = 1,
                        window_size = 10,
                        grad_conf_level = 0.75,
                        cov_matrix_function = MLE_correlation_matrix) {
        
        start_time <- as.numeric(as.POSIXlt(Sys.time()))
        
        N <- nrow(data)
        p <- ncol(data)
        
        # Check for some common mistakes in the input
        if (length(marginals) != p) {
                stop("Dimension of the data matrix does not match the number of marginals.")
        }
        
        if (!is.null(initial_cov) && (dim(initial_cov)[[1]] != p || dim(initial_cov)[[2]] != p)) {
                stop("Dimension of the data matrix does not match the dimension of the initial correlation matrix.")
        }
        
        if (! is.null(unobserved_data) && fixed_jitters == FALSE) {
                stop("When using predetermined unobserved_data fixed_jitters must be set to TRUE.")
        }
        
        if ((is.null(mcem_samples) + is.null(unobserved_data)) != 1){
                stop("Either mcem_samples or unobserved_data must be fixed, and not both.")
        }
        
        if (is.null(mcem_samples) && ! is.null(unobserved_data)) {
                mcem_samples <- nrow(unobserved_data) / N
        }
        
        if (max_iterations < window_size) {
                stop("max_iterations must be equal or greater than window_size")
        }
        
        if (window_size < 2) {
                stop("The window_size must be greater than 1.")
        }
        
        # Set-up the list which will be returned by this function
        result <- structure(list(), class = "copuladiscreter_result")
        result[["cov_initial"]] <- NA
        result[["cov_initial_log_lik"]] <- NA
        result[["cov_final"]] <- NA
        result[["cov_final_var"]] <- NA
        result[["cov_final_log_lik"]] <- NA
        result[["cov_final_log_lik_var"]] <- NA
        result[["convergence"]] <- 1
        result[["iterations"]] <- NA
        result[["mcem_samples"]] <- mcem_samples
        result[["unobserved_data"]] <- ifelse(is.null(unobserved_data), NA, unobserved_data)
        result[["data"]] <- data
        result[["marginals"]] <- marginals
        result[["window_size"]] <- window_size
        result[["grad_conf_level"]] <- grad_conf_level
        
        result[["cov_steps"]] <- list()
        result[["cov_steps_var"]] <- list()
        result[["cov_steps_log_lik"]] <- list()
        result[["cov_steps_log_lik_var"]] <- list()
        result[["cov_steps_start_time"]] <- list()
        
        # If no initial covariance matrix is supplied, then just use an
        # identity matrix. Don't use anything more complicated, so that
        # it is still a valid result for a wide variety of covariance
        # modeling structures.
        if (is.null(initial_cov)) {
                cov_matrix <- diag(rep(1, p))
        }
        else {
                cov_matrix <- initial_cov
        }
        
        colnames(cov_matrix) <- colnames(data)
        rownames(cov_matrix) <- colnames(data)
        
        result[["cov_initial"]] <- cov_matrix
        
        # Perform an initial MCEM expectation step using the initial covariance matrix
        ##  Uniformly sample the latent data, combine with the observed data and map into the MVN space.
        mc_sampled_data_MVN <- resample_data(data, marginals, mcem_samples, unobserved_data)
        
        ##  Calculate the weighting for each MC sample
        mcem_weights <- dgauss_ratio(mc_sampled_data_MVN, cov_matrix)
        
        ##  Calculate the normalisation constant for each set of MC samples
        ##   Store these values, so that we can use them to efficiently calculate
        ##   the log likelihood of the current cov_matrix
        mcem_norms <- calc_mcem_norms(mcem_weights, N)
        
        ## Normalise the weights
        mcem_weights <- normalise_weights(mcem_weights, mcem_norms)
        # End of initial expectation step
        
        # Use the MC sample normalisation constants to efficiently estimate the log Likelihood of cov_matrix
        result[["cov_initial_log_lik"]] <- log_likelihood(data, marginals, mcem_norms = mcem_norms)
        
        if (verbosity > 0) {
                print("Starting MCEM procedure...")
                pb <- utils::txtProgressBar(style = 3)
        }
        
        for (iteration in 1:max_iterations) {
                result[["cov_steps_start_time"]][[iteration]] <- as.numeric(as.POSIXlt(Sys.time())) - start_time
                
                if (verbosity > 0) {
                        utils::setTxtProgressBar(pb, iteration / max_iterations)
                }
                
                # The MCEM maximisation step
                #   Call a function which finds the MLE of the covariance matrix given
                #   latent data with a distribution determined by the previous covariance matrix estimate.
                cov_matrix <- cov_matrix_function(mc_sampled_data_MVN,
                                                  mcem_weights,
                                                  optim_params = list("initial_cov" = cov_matrix))
                
                result[["cov_steps"]][[iteration]] <- cov_matrix
                
                if (! fixed_jitters) {
                        mc_sampled_data_MVN <- resample_data(data, marginals, mcem_samples)
                }
                
                # The MCEM expectation step for this iteration
                mcem_weights <- dgauss_ratio(mc_sampled_data_MVN, cov_matrix)
                mcem_norms <- calc_mcem_norms(mcem_weights, N)
                mcem_weights <- normalise_weights(mcem_weights, mcem_norms)
                
                # Calculate the log Likelihood for the current cov_matrix
                result[["cov_steps_log_lik"]][[iteration]] <- log_likelihood(data, marginals, mcem_norms = mcem_norms)
                
                # Halt the algorithm if the log likelihood appears to have converged.
                if (iteration >= window_size) {
                        window <- (iteration - window_size + 1):iteration
                        result[["cov_steps_var"]][[iteration]] <- matrix_var(result[["cov_steps"]][window])
                        result[["cov_steps_log_lik_var"]][[iteration]] <- var(unlist(result[["cov_steps_log_lik"]][window]))
                        if (! is.na(grad_conf_level) && check_convergence(result[["cov_steps_log_lik"]][window], grad_conf_level)) {
                                result[["convergence"]] <- 0
                                break
                        }
                }
                else {
                        result[["cov_steps_var"]][[iteration]] <- NA
                }
        }
        
        if (verbosity > 0) {
                close(pb)
                print("Completed MCEM iterations")
        }
        
        result[["iterations"]] <- iteration
        
        
        if (result[["convergence"]] != 0) {
                result[["convergence"]] <- 1
                warning("MCEM algorithm reached the maximum number of iterations without converging.")
        }
        
        # The covariance matrix returned by the final iteration of the MCEM algorithm is not necassarily the
        #  the best one found by the algorithm. Search through the window and return the covariance matrix
        #  which has the highest Likelihood.
        best_iteration <- which.max(result[["cov_steps_log_lik"]][(iteration - window_size + 1):iteration]) + iteration - window_size
        result[["cov_final"]] <- result[["cov_steps"]][[best_iteration]]
        result[["cov_final_log_lik"]] <- result[["cov_steps_log_lik"]][[best_iteration]]
        result[["cov_final_var"]] <- result[["cov_steps_var"]][[iteration]]
        result[["cov_final_log_lik_var"]] <- result[["cov_steps_log_lik_var"]][[iteration]]
        
        result[["total_time"]] <- as.numeric(as.POSIXlt(Sys.time())) - start_time
        
        return(result)
}

##############################################################
# The rest of this file contains utility functions called by #
# the main mcem_copula() function and used for testing       #
##############################################################

#' Create random correlation matrices
#'
#' \code{rcorrelation_matrix} creates random correlation matrices
#' based on covariance matrices taken from the Wishart distribution
#' \eqn{W_p \left( I_n, p \right)}
#'
#' @param n The number of matrices to generate
#' @param p The dimension of each matrices
#'
#' @return A p by p by n multi-dimension array. The [,,1] slice of
#' the array will be the first correlation matrix generated.
#' @export
rcorrelation_matrix <- function(n, p) {
        covariance_matrices <- stats::rWishart(n, p, diag(rep(1, p)))
        correlation_values <- apply(covariance_matrices, 3, stats::cov2cor)
        return(array(correlation_values, dim = c(p, p, n)))
}

#' Generate Data from a Gaussian Copula
#'
#' \code{generate_data} generates random samples from a gaussian
#' copula with the given correlation matrix and a set of marginals.
#'
#' @param N The number of samples to generate
#' @param marginals A list of p marginals
#' @param corr A p by p matrix which represents the correlation
#'             structure of the data in the multi-variate normal space.
#'
#' @return An N by p matrix, each row of which is a sample from the gaussian copula.
#' @export
generate_data <- function(N, marginals, corr) {
        if (length(dim(corr)) != 2 || dim(corr)[[1]] != dim(corr)[[2]]) {
                stop("corr is not a square matrix.")
        }
        
        p <- dim(corr)[[1]]
        
        if (length(marginals) != p) {
                stop("Dimension of the correlation matrix does not match the number of marginals.")
        }
        
        data_MVN <- mnormt::rmnorm(n = N, varcov = corr)
        if (N == 1) {
                data_MVN <- t(data_MVN)
        }
        data_uniform <- prob_int_trans_stdnorm(data_MVN)
        data <- prob_int_trans(data_uniform, marginals, inverse = TRUE)
        
        data_observed <- floor(data)
        data_censored <- data - data_observed
        
        if (N > 1) {
                cov_actual <- stats::cov.wt(data_MVN, method = "ML", center = FALSE)[["cov"]]
        }
        else {
                cov_actual <- NA
        }
        return(list(observed = data_observed, censored = data_censored, cov_actual = cov_actual))
}

#' Resample Data with Uniform Jitters
#'
#' This function takes an integer-valued dataset and replicates
#' it (a given number of times) adding fractional components,
#' uniformly sampled from \eqn{[0, 1)}, to each value
#' and then applies the marginals probability integral transform (PIT),
#' and the inverse of the standard normal PIT returning
#' the data in the gaussian copula's multi-variate normal (MVN) space.
#'
#' @param data An N by p matrix of integer-valued samples
#' @param marginals A list of p marginals
#' @param replicates The number of times to resample each sampling unit
#' @param unobserved_data A set of predetermined latent data. If this is supplied
#'   then replicates should be set to NULL.
#'
#' @return A N * replicates by p matrix where each set of N contiguous rows
#'         represents a single replication of the original dataset with latent data added
#'         tranformed under the marginals PIT and then the inverse of the standard normal
#'         PIT into the space gaussian normal's MVN space.
resample_data <- function(data, marginals, replicates = NULL, unobserved_data = NULL) {
        N <- nrow(data)
        p <- ncol(data)
        
        if (length(marginals) != p) {
                stop("Dimension of the data matrix does not match the number of marginals.")
        }
        
        if (! is.null(unobserved_data)) {
                if (is.null(replicates)) {
                        replicates <- nrow(unobserved_data) / N
                }
                
                if (nrow(unobserved_data) != N * replicates) {
                        stop("unobserved_data does not contain the correct number of data points")
                }
                
                if (ncol(unobserved_data) != p) {
                        stop("unobserved_data does not contain the correct number of variables")
                }
        }
        
        data_uniform <- prob_int_trans(data, marginals)
        data_uniform_maxima <- prob_int_trans(data + 1, marginals)
        
        data_uniform_hypercube <- data_uniform_maxima - data_uniform
        data_uniform_hypercube <- rep(1, replicates) %x% data_uniform_hypercube
        
        data_uniform_resampled <- rep(1, replicates) %x% data_uniform
        
        if (is.null(unobserved_data)) {
                data_uniform_resampled <- data_uniform_resampled + stats::runif(N * p * replicates) * data_uniform_hypercube
        }
        else {
                data_uniform_resampled <- data_uniform_resampled + unobserved_data * data_uniform_hypercube
        }
        
        data_MVN_resampled <- prob_int_trans_stdnorm(data_uniform_resampled, inverse = TRUE)
        colnames(data_MVN_resampled) <- colnames(data)
        
        return(data_MVN_resampled)
}

corr_ML_condensed_data <- function(condensed_data, N, optim_params = list()) {
        p <- dim(condensed_data)[[1]]
        
        if (is.null(optim_params[["initial_cov"]])) {
                initial_corr <- stats::cov2cor(1 / N * condensed_data)
        }
        else {
                initial_corr <- optim_params[["initial_cov"]]
                optim_params[["initial_cov"]] <- NULL
        }
        
        initial_spherical_params <- corr_to_spherical(initial_corr)
        
        f <- function(spherical_params) {
                corr <- spherical_to_corr(spherical_params, p)
                result <- log_likelihood_MVN_condensed_data(corr, condensed_data, N)
                if (result == -Inf) {
                        # Return the negative sqrt of the largest double precision number.
                        # Don't return the actual largest double precision number so
                        # that a derivative based optimiser still has a chance of numericaly
                        # calculating a representable derivative.
                        return(-sqrt(.Machine[["double.xmax"]]))
                }
                if (result == Inf) {
                        stop("Error: found a covariance matrix with infinite likelihood!")
                }
                
                return(result)
        }
        
        # A small epsilon is used to prevent the optimiser straying into a region
        # with a correlation matrix containing a pair of completely correlated variables.
        # Such a correlation matrix will have a zero likelihood of producing most data and
        # the resulting log_likelihood = -Inf cannot be understood by the optimiser.
        epsilon <- .Machine[["double.eps"]]
        lower <- rep(0 + epsilon, length(initial_spherical_params))
        upper <- rep(pi - epsilon, length(initial_spherical_params))
        
        if (length(initial_spherical_params) > 1) {
                optim_value <- dfoptim::hjkb(initial_spherical_params,
                                             f,
                                             lower = lower,
                                             upper = upper,
                                             control = list(maximize = TRUE))
        }
        else {
                optim_params[["fnscale"]] <- -1
                optim_params[["trace"]] <- 0
                
                optim_value <- stats::optim(initial_spherical_params,
                                            f,
                                            method = "L-BFGS-B",
                                            lower = lower,
                                            upper = upper,
                                            control = optim_params)
        }
        
        if (optim_value[["convergence"]] != 0) {
                warning("convergence failure")
        }
        
        result <- spherical_to_corr(optim_value[["par"]], p)
        
        colnames(result) <- colnames(condensed_data)
        rownames(result) <- rownames(condensed_data)
        
        return(result)
}

mean_covariance_matrix <- function(data, weights = NULL) {
        N <- nrow(data)
        p <- ncol(data)
        
        if (is.null(weights)) {
                weights <- rep(1 / N, N)
        }
        else {
                weights <- weights / sum(weights)
        }
        
        if (length(weights) != N) {
                stop("Length of weights is not compatible with the dimension of data.")
        }
        
        # Simple, but slow technique
        #result <- matrix(0, nrow = p, ncol = p)
        #for (i in 1:N) {
        #   result <- result + weights[[i]] * t(data[i,, drop = FALSE]) %*% data[i,, drop = FALSE]
        #}
        
        # Fast, but hard to understand technique
        result <- crossprod(data, sweep(data, MARGIN = 1, STATS = weights, FUN = '*'))
        
        colnames(result) <- colnames(data)
        rownames(result) <- colnames(data)
        
        return(result)
}

#' The log likelihood of a MVN given data in condensed form.
#'
#' If this function is maximised with respect to the \code{cov}
#' parameter then the maximum likelihood estimator (MLE) of
#' the covariance matrix (\code{cov}) will be found. Note that the
#' data is supplied in a condensed form, not as the full data matrix.
#' All population means are assumed to be zero.
#'
#' @param cov The p by p covariance matrix
#' @param data_matrix The data condensed into a p by p matrix.
#'                    \eqn{\text{data_matrix} = \frac{1}{\sum_{i = 1}^N w_i}\sum_{i = 1}^N w_i x_i^\top x_i}
#' @param N The number of samples in the original dataset before it was condensed
#'
#' @return A value equal to the log likelihood with a constant offset.
log_likelihood_MVN_condensed_data <- function(cov, data_matrix, N) {
        if (all(dim(cov) != dim(data_matrix))) {
                stop("Dimensions of the covariance matrix do not match the condensed data matrix.")
        }
        
        # The covariance matrix is singular.
        # Therefore we will
        if (det(cov) < .Machine[["double.eps"]]) {
                return(-Inf)
        }
        
        p <- dim(cov)[1]
        
        #result <- - N / 2 * (p * log(2 * pi) + log(det(cov)) + sum(diag(data_matrix %*% solve(cov))))
        result <- try(- N / 2 * (p * log(2 * pi) + log(det(cov)) + sum(diag(data_matrix %*% solve(cov)))))
        
        if (class(result) == "try-error") {
                stop("Failed to calculate log likelihood")
        }
        
        return(result)
}

log_likelihood_MVN <- function(cov, data, weights = NULL) {
        N <- nrow(data)
        p <- ncol(data)
        
        if (is.null(weights)) {
                weights <- rep(1 / N, N)
        }
        else {
                weights <- weights / sum(weights)
        }
        
        if (length(weights) != N) {
                stop("Length of weights is not compatible with the dimension of data.")
        }
        
        cov_inv <- solve(cov)
        
        result <- -N / 2 *  p * log(2 * pi) - N / 2 * log(det(cov))
        
        for (i in 1:N) {
                result <- result - N / 2 * weights[[i]] * sum(diag(data[i,, drop = FALSE] %*% cov_inv %*% t(data[i,, drop = FALSE])))
        }
        
        return(result)
}

mag <- function(x) {
        return(sqrt(sum(x * x)))
}

spherical_to_euclidean <- function(spherical_vector) {
        N <- length(spherical_vector)
        
        if (N == 1) {
                return(spherical_vector)
        }
        
        sines   <- c(1, sin(spherical_vector[2:N]))
        cosines <- c(cos(spherical_vector[2:N]), 1)
        
        return(spherical_vector[1] * cumprod(sines) * cosines)
}

spherical_angles_to_euclidean <- function(spherical_angles) {
        sines   <- c(1, sin(spherical_angles))
        cosines <- c(cos(spherical_angles), 1)
        
        return(cumprod(sines) * cosines)
}

euclidean_to_spherical <- function(euclidean_vector) {
        N <- length(euclidean_vector)
        
        if (N == 1){
                return(euclidean_vector)
        }
        
        spherical_vector <- rep(NA, N)
        
        spherical_vector[[1]] <- mag(euclidean_vector)
        
        for (i in 2:N) {
                spherical_vector[[i]] <- acos(euclidean_vector[[i - 1]] / mag(euclidean_vector[(i - 1):N]))
        }
        
        if (euclidean_vector[[N]] < 0) {
                spherical_vector[[N]] <- 2 * pi - spherical_vector[[N]]
        }
        
        return(spherical_vector)
}

spherical_to_corr <- function(spherical_params, N) {
        if(length(spherical_params) != (N * (N - 1) / 2)) {
                stop("Incorrect number of spherical params for this dimension")
        }
        
        # Everything is shifted up 1 in effect because the first element we want must be 1
        #split_spherical <- split(spherical_params, rep(1:(N - 1), 1:(N - 1)))
        
        split_spherical <- as.list(rep(NA, N - 1))
        start <- 1
        for (i in 1:(N - 1)) {
                split_spherical[[i]] <- spherical_params[start:(start + i - 1)]
                start <- start + i
        }
        
        split_euclidean <- lapply(split_spherical, spherical_angles_to_euclidean)
        
        # Add the implicit first unit vector.
        euclidean <- c(1, unlist(split_euclidean, FALSE, FALSE))
        # Faster, but not guaranteed to work in future versions of R
        #euclidean <- c(1, .Internal(unlist(split_euclidean, FALSE, FALSE)))
        
        # Create the correlation matrix by applying the inverse Cholesky factorisation.
        U <- matrix(0, nrow = N, ncol = N)
        U[upper.tri(U, diag = TRUE)] <- euclidean
        corMatrix <- t(U) %*% U
        
        diag(corMatrix) <- rep(1, N)
        return(corMatrix)
}

corr_to_spherical <- function(corr) {
        if (length(dim(corr)) != 2 || (dim(corr)[[1]] != dim(corr)[[2]])) {
                stop("corr is not a square matrix")
        }
        
        N <- nrow(corr)
        
        if (! isTRUE(all.equal(diag(corr), rep(1, length(diag(corr))), check.names = FALSE))) {
                stop("corr does not have unit variances")
        }
        
        if (N == 1) {
                warning("datasets with only a single variable are not supported by this package")
                return(NULL)
        }
        
        U <- chol(corr)
        
        spherical_coords <- as.list(rep(NA, N))
        
        # Only start with the second dimension as the first has length zero in the spherical paramerisation
        for (i in 2:N) {
                spherical_coords[[i]] <- euclidean_to_spherical(U[,i][1:i])
                
                # Remove the implicit first entry which equals the variance (1)
                spherical_coords[[i]] <- spherical_coords[[i]][2:i]
        }
        
        return(unlist(spherical_coords[2:N]))
}

reorder_jitters <- function(jitters) {
        N <- length(jitters)
        p <- ncol(jitters[[1]])
        resamples <- nrow(jitters[[1]])
        
        for (i in 1:N) {
                if (nrow(jitters[[i]]) != resamples) {
                        stop("Different number of jitters for each data sample unit")
                }
                if (ncol(jitters[[i]]) != p) {
                        stop("Different number of variables for each data sample unit")
                }
        }
        
        result <-  do.call(rbind, jitters)[order(sequence(rep(resamples, N))), ]
        
        # Ensure the result is a matrix (and not a R-vector)
        if (is.null(dim(result))) {
                result <- t(result)
        }
        
        # The direct optimiser jitters are deducted from the maxima
        # so transform them so they can be added to the minima as in the MCEM code.
        result <- 1 - result
        
        return(result)
}

#' The Log Likelihood of a Gaussian Copula
#'
#' @param data An N by p matrix containing the observed count data
#' @param marginals A list of p marginals
#' @param cov The gaussian copula's correlation matrix
#' @param mc_samples The number of resamples to perform during MC-integration for each sample
#' @param mcem_norms The normalisation contants used in the MCEM calculation. If this parameter
#'   supplied then mc_samples should be NULL. Using precomputed mcem_norms significantly reduces
#'   the computational time required by this function.
#'
#' @return The log likelihood of the copula parameters given the data
log_likelihood <- function(data, marginals, cov = NULL, mc_samples = NULL, mcem_norms = NULL)
{
        N <- nrow(data)
        p <- ncol(data)
        
        if (! is.null(mcem_norms) && (! is.null(cov) || ! is.null(mc_samples))) {
                stop("If mcem_norms is supplied then cov and resamples should not be.")
        }
        
        if (is.null(mcem_norms) && (is.null(cov) || is.null(mc_samples))) {
                stop("cov and mc_samples must be supplied (unless mcem_norms is).")
        }
        
        if (is.null(mcem_norms)) {
                resampled_data_MVN <- resample_data(data, marginals, mc_samples)
                resampled_data_dgauss_ratio <- dgauss_ratio(resampled_data_MVN, cov)
                mcem_norms <- calc_mcem_norms(resampled_data_dgauss_ratio, N)
        }
        
        marginal_dens <- matrix(nrow = N, ncol = p)
        
        for (j in 1:p) {
                # The contribution from each marginal
                marginal_dens[ , j] <- do.call(marginals[[j]][["pdf"]], c(list(x = data[, j]), marginals[[j]][["parameters"]]))
        }
        
        marginal_dens <- apply(marginal_dens, 1, prod)
        dim(marginal_dens) <- c(N, 1)
        
        result <- sum(log(marginal_dens) + log(mcem_norms))
        
        return(result)
}

log_likelihood_from_latent <- function(data, marginals, cov, latent_data) {
        N <- nrow(data)
        data_resampled_MVN <- resample_data(data, marginals, unobserved_data = latent_data)
        sample_unit_gaussian_dens <- dgauss_ratio(data_resampled_MVN, cov)
        mcem_norms <- calc_mcem_norms(sample_unit_gaussian_dens, N)
        return(log_likelihood(data, marginals, mcem_norms = mcem_norms))
}

dgauss_ratio <- function(data, cov) {
        result <- mnormt::dmnorm(data, varcov = cov) / apply(stats::dnorm(data), 1, prod)
        
        dim(result) <- c(nrow(data), 1)
        
        return(result)
}

calc_mcem_norms <- function(dgauss_ratio_data, N) {
        resamples <- nrow(dgauss_ratio_data) / N
        
        dim(dgauss_ratio_data) <- c(N, resamples)
        result <- apply(dgauss_ratio_data, 1, sum)
        result <- 1 / resamples * result
        dim(result) <- c(N, 1)
        
        return(result)
}

normalise_weights <- function(dgauss_ratio_data, norms) {
        resamples <- nrow(dgauss_ratio_data) / nrow(norms)
        return(dgauss_ratio_data / (rep(1, resamples) %x% norms))
}

check_convergence <- function(log_likelihood, level) {
        x <- 1:length(log_likelihood)
        y <- unlist(log_likelihood)
        data_to_fit <- data.frame(x, y)
        linear_model <- stats::lm(y ~ x, data_to_fit)
        two_tailed_level <- 1 - (1 - level) * 2
        grad_LB <- stats::confint(linear_model, 2, level)[[1]]
        
        if (grad_LB <= 0) {
                return(TRUE)
        }
        else {
                return(FALSE)
        }
}

matrix_var <- function(list_of_matrices) {
        cov_array <- simplify2array(list_of_matrices)
        
        return(apply(cov_array, c(1, 2), stats::var))
}

#' Generate a marginal
#'
#' \code{generate_marginal} takes a marginal specified by a string and a
#' vector of parameters and converts it into the list representation of a
#' marginal used within this package.
#'
#' @param name A string representing the name of the marginal.
#'     Can be "Pois", "NB", "ZIP" or "ZINB"
#' @param parameters A vector of three numbers (mu, theta, pi) which are
#'     the parameters for the distribution.
#'
#' @return The list (pdf, cdf, qf, parameters) where parameters is a named list.
generate_marginal <- function(name, parameters) {
        
        if (toupper(name) == "POIS") {
                return(list(pdf = stats::dpois, cdf = stats::ppois, qf = stats::qpois, parameters = list(lambda = parameters[[1]])))
        }
        else if (toupper(name) == "NB") {
                return(list(pdf = stats::dnbinom, cdf = stats::pnbinom, qf = stats::qnbinom, parameters = list(mu = parameters[[1]],
                                                                                                               size = parameters[[2]])))
        }
        else if (toupper(name) == "ZIP") {
                if (!is.na(parameters[[2]])) {
                        stop("second parameter for a ZIP must be NA")
                }
                return(list(pdf = ZIM::dzip, cdf = ZIM::pzip, qf = ZIM::qzip, parameters = list(lambda = parameters[[1]],
                                                                                                omega = parameters[[3]])))
        }
        else if (toupper(name) == "ZINB") {
                return(list(pdf = ZIM::dzinb, cdf = ZIM::pzinb, qf = ZIM::qzinb, parameters = list(lambda = parameters[[1]],
                                                                                                   k = parameters[[2]],
                                                                                                   omega = parameters[[3]])))
        }
        else {
                stop(sprintf("%s distribution type is not yet implemented.", name))
        }
}

#######################################################
# Functions below relate to the handling of marginals #
#######################################################

#' Create a continuous CDF for a integer-valued random variable.
#'
#' \code{smoothed_cdf} creates a continuous cumulative density function from
#' an integer-valued random variable by distributing each integer's
#' probability density uniformly between it and the next integer.
#'
#' @param discrete_pdf The pdf of the integer-valued random variable.
#' @param discrete_cdf The cdf of the integer-valued random variable.
#' @param parameters The parameters for \code{discrete_pdf} and \code{discrete_cdf}.
#'
#' @return A function from the real numbers to the unit interval which
#'         is a continuous version of the discrete random variable's cumulative
#'         density function.
smoothed_cdf <- function(discrete_pdf, discrete_cdf, parameters) {
        return( function(x) {
                x_floor <- floor(x)
                
                cd_base <- do.call(discrete_cdf, c(list(q = x_floor - 1), parameters))
                cd_derivative <- do.call(discrete_pdf, c(list(x = x_floor), parameters))
                
                return(cd_base + (x - x_floor) * cd_derivative)
        })
}

#' Create a continuous quantile function for a integer-valued random variable.
#'
#' \code{smoothed_cdf} creates a continuous quantile function from
#' an integer-valued random variable by distributing each integer's
#' probability density uniformly between it and the next integer.
#'
#' @param discrete_pdf The pdf of the integer-valued random variable.
#' @param discrete_cdf The cdf of the integer-valued random variable.
#' @param discrete_qf The quantile of the integer-valued random variable.
#' @param parameters The parameters for \code{discrete_pdf},
#'                   \code{discrete_cdf} and \code{discrete_qf}.
#'
#' @return A function from the the unit interval to the real numbers which
#'         is a continuous version of the discrete random variable's quantile
#'         function.
smoothed_qf <- function(discrete_pdf, discrete_cdf, discrete_qf, parameters) {
        return( function(p) {
                x_base <- do.call(discrete_qf, c(list(p = p), parameters))
                x_derivative <- 1 / do.call(discrete_pdf, c(list(x = x_base), parameters))
                p_delta <- p - do.call(discrete_cdf, c(list(q = x_base - 1), parameters))
                
                return(x_base + p_delta * x_derivative)
        })
}

#' Standard Normal Probability Integral Transform
#'
#' \code{prob_int_trans_stdnorm} applies the standard normal distribution
#' probability integral transform (PIT), or its inverse, to a set of multivariate
#' data stored in a matrix.
#'
#' @param data the data
#' @param inverse whether to apply the PIT or its inverse
#'
#' @return the transformed data
prob_int_trans_stdnorm <- function(data, inverse = FALSE) {
        if (!inverse) {
                return(stats::pnorm(data))
        }
        else {
                return(stats::qnorm(data))
        }
}

#' Probability Integral Transform
#'
#' \code{prob_int_trans_stdnorm} applies the probability integral transform
#' of the supplied marginals to the supplied dataset
#'
#' @param data the data
#' @param marginals a list of marginal lists
#' @param inverse whether to apply the PIT or its inverse
#'
#' @return the transformed data
prob_int_trans <- function(data, marginals, inverse = FALSE) {
        N <- dim(data)[[1]]
        p <- dim(data)[[2]]
        
        if (length(marginals) != p) {
                stop("Dimension of the data matrix does not match the number of marginals.")
        }
        
        data_trans <- matrix(nrow = N, ncol = p)
        
        for (j in 1:p) {
                if (!inverse) {
                        func <- smoothed_cdf(marginals[[j]][["pdf"]], marginals[[j]][["cdf"]], marginals[[j]][["parameters"]])
                }
                else {
                        func <- smoothed_qf(marginals[[j]][["pdf"]], marginals[[j]][["cdf"]], marginals[[j]][["qf"]], marginals[[j]][["parameters"]])
                }
                data_trans[,j] <- func(data[,j])
        }
        
        return(data_trans)
}


