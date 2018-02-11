#
set.seed(12345)
X <- matrix(rnorm(10 * 3), 10, 3)

#
mdm(X, type = "asym_dcov")$stat
mdm(X, type = "sym_dcov")$stat
mdm(X, type = "comp")$stat
mdm(X, type = "comp_simp")$stat
mdm(X, type = "asym_comp")$stat
mdm(X, type = "asym_comp_simp")$stat
mdm(X, type = "sym_comp")$stat
mdm(X, type = "sym_comp_simp")$stat

# cd mmi/simulation/real_data/ff_5_year
source('dCov.R')
set.seed(12345)
X <- matrix(rnorm(10 * 3), 10, 3)
num_obs <- nrow(X)
dim_comp <- c(1, 1, 1)
num_dim <- sum(dim_comp)
num_comp <- length(dim_comp)
index_comp <- cumsum(c(1, dim_comp)) # start index and end index
X <- c(t(X))
D <- rep(0, num_comp * num_obs * num_obs)

#
dCov_asymmetric(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
dCov_symmetric(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_complete(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_complete_simple(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_asymmetric(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_asymmetric_simple(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_symmetric(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q
est_symmetric_simple(X, D, 0, num_obs, num_dim, num_comp, index_comp - 1)$Q


