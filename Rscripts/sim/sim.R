#
# set random seed for task relative to entire job
#

# set seed for job
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)

# set seed for task
taskId = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
s <- .Random.seed
for (i in 1:taskId) {
  s <- parallel::nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s

# determine simulation indices for task
nsim = 1e4
ncores = 400
tasks = parallel::splitIndices(nsim, ncores)
sim_inds = tasks[[taskId]]

#
# workspace configuration
#

library(dplyr)

# load summary statistic fn.
source('R/eda/sattag_summary.R')

# load additional components for simulation
source('R/data/segment_fn.R')
targets::tar_load(template_bins)
bin_starts = template_bins$center - template_bins$halfwidth

# load all window configurations
targets::tar_load(study_windows)

#
# simulation code
#

# number of timepoints to simulate (convert days to 5 min observations)
nbaseline = 7     * 24 * 60 / 5
nexposed = (6/24) * 24 * 60 / 5

# transition matrix for sequence of dive depths
Tm = matrix(
  data = c(.9, .1, .25, .75), nrow = 2, byrow = TRUE, 
  dimnames = list(c('shallow', 'deep'), c('shallow','deep'))
)

# transition matrix for sequence of dive depths after exposure (shallow only)
Texposed = matrix(
  data = c(1, 0, 1, 0), nrow = 2, byrow = TRUE, 
  dimnames = list(c('shallow', 'deep'), c('shallow','deep'))
)

# extract information about simulation scenario
stationary = eigen(t(Tm))$vectors[,1]
stationary = stationary / sum(stationary)
nstates = nrow(Tm)

# simulation parameters for depths
mu = c(shallow = 100, deep = 1e3)
sigma = c(shallow = 100, deep = 100)

# function to simulate a dive trajectory
sim = function(n, nexposed, response) {
  # Parameters:
  #  n - length of baseline series to simulate
  #  nexposed - length of exposed series to simulate
  #  response - TRUE to simulate data with a response
  
  # ergodic distribution for states
  p = eigen(Tm)$vectors[,1]
  
  # sample initial state and depth
  states = sample(x = nstates, size = 1, prob = stationary)
  depths = rnorm(n = 1, mean = mu[states[1]], sd = sigma[states[1]])
  
  # sample additional baseline depths
  for(i in 2:n) {
    states[i] = sample(x = nstates, size = 1, prob = Tm[states[i-1],])
    depths[i] = rnorm(n = 1, mean = mu[states[i]], sd = sigma[states[i]])
  }
  
  # set the transition matrix used for state transitions for exposed conditions
  if(response) {
    Tm2 = Texposed
  } else {
    Tm2 = Tm
  }
  
  # sample exposed depths
  if(nexposed > 0) {
    for(i in (n + 1:nexposed)) {
      states[i] = sample(x = nstates, size = 1, prob = Tm2[states[i-1],])
      depths[i] = rnorm(n = 1, mean = mu[states[i]], sd = sigma[states[i]])
    }
  }
  
  # truncate to surface
  depths[depths < 0] = 0
  
  # return depths
  depths
}

# function to analyze time series of depths
prepost = function(x, nwin, nlag, baseline_end) {
  # Parameters:
  #  x - time series of depths
  #  nwin - number of timepoints in each pre/post window
  #  nlag - number of timepoints between each pre/post window
  #  baseline_end - index of last baseline observation
  
  # observed pre/post
  pre_inds = seq(to = baseline_end, length.out = nwin)
  post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
  
  # validate that enough data is present for analysis
  if(max(c(pre_inds,post_inds)) > length(x)) {
    stop('Not enough data to run analysis')
  }
  
  #
  # enrich time series
  #
  
  # initial features
  d = data.frame(depths = x) %>% 
    mutate(
      depth.bin = findInterval(x = depths, vec = bin_starts),
      ddepths = c(0, diff(depths)),
      ddepths.sign = sign(ddepths),
      descending = ddepths.sign > 0,
      ascending = ddepths.sign < 0,
      no_change = ddepths.sign == 0
    )
  
  # segment dives
  labels = dive.segmentation(
    y = d$depth.bin, merge.ratio = .6, depth.bins = template_bins, 
    times = as.POSIXct(1:nrow(d), origin = '1970-01-01 00:00.00 UTC'), 
    timestep = 1
  )
  
  # compute dive types
  diveTypes = d %>% 
    mutate(diveId = labels) %>% 
    group_by(diveId) %>% 
    summarise(
      diveType = ifelse(max(depths) > 800, 'Deep', 'Shallow')
    )
  
  # merge dive types
  d = d %>% 
    mutate(diveId = labels) %>% 
    left_join(diveTypes, by = 'diveId')
  
  #
  # compute window summaries
  #
  
  # observed results
  observed_pre = sattag_summary(d[pre_inds, ])
  observed_post = sattag_summary(d[post_inds, ])
  
  # null distribution samples
  null_samples_pre = do.call(rbind, lapply(1:baseline_end, function(start) {
    # determine window inds
    pre_inds = seq(from = start, length.out = nwin)
    post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
    # skip processing if any indices are outside baseline window
    if(max(c(pre_inds,post_inds)) > nbaseline) {
      return(NULL)
    }
    # compute summaries
    sattag_summary(d[pre_inds,])
  }))
  null_samples_post = do.call(rbind, lapply(1:baseline_end, function(start) {
    # determine window inds
    pre_inds = seq(from = start, length.out = nwin)
    post_inds = seq(from = tail(pre_inds, 1) + nlag + 1, length.out = nwin)
    # skip processing if any indices are outside baseline window
    if(max(c(pre_inds,post_inds)) > nbaseline) {
      return(NULL)
    }
    # compute summaries
    sattag_summary(d[post_inds,])
  }))
  
  #
  # compute tail probabilities
  #
  
  probs = do.call(rbind, lapply(1:ncol(null_samples_pre), function(ind) {
    
    # bandwidths for kde, modified s.t. results are always non-zero if data
    # are not degenerate
    h = c(MASS::bandwidth.nrd(null_samples_pre[,ind]),
          MASS::bandwidth.nrd(null_samples_post[,ind]))
    if(h[1] == 0) {
      x = null_samples_pre[,ind]
      h[1] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    if(h[2] == 0) {
      x = null_samples_post[,ind]
      h[2] = 4 * 1.06 * sqrt(var(x)) * length(x)^(-1/5)
    }
    
    # bivariate kernel density estimate (KDE)
    dens = MASS::kde2d(
      x = null_samples_pre[,ind], 
      y = null_samples_post[,ind], 
      n = 500, 
      lims = c(range(c(null_samples_pre[,ind], observed_pre[[ind]])), 
               range(c(null_samples_post[,ind], observed_post[[ind]]))),
      h = h
    )
    
    # indexes to retrieve conditional distributions for KDE 
    kde_ind_pre = which.min(abs(dens$x - observed_pre[[ind]]))
    kde_ind_post = which.min(abs(dens$y - observed_post[[ind]]))
    
    # conditional distribution P(Post | Pre) and complement
    C = sum(dens$z[kde_ind_pre,])
    cdf.conditional = cumsum(dens$z[kde_ind_pre,]) / C
    ccdf.conditional = rev(cumsum(rev(dens$z[kde_ind_pre,]))) / C
    
    
    # bivariate p-values
    data.frame(
      p.left = cdf.conditional[kde_ind_post],
      p.right = ccdf.conditional[kde_ind_post],
      stat = factor(colnames(null_samples_pre)[ind])
    )
  }))
  
  probs
}

# run a number of simulations with responses
sim_res = do.call(rbind, lapply(sim_inds, function(rep_ind) {
  # local random seed, to simulate datasets from same random variates
  seed = .GlobalEnv$.Random.seed
  xr = sim(n = nbaseline, nexposed = nexposed, response = TRUE)
  .GlobalEnv$.Random.seed = seed
  xb = sim(n = nbaseline, nexposed = nexposed, response = FALSE)
  # analyze datasets for all window combinations
  do.call(rbind, lapply(1:nrow(study_windows), function(wid) {
    # extract window config
    wlen = study_windows$window_length[wid]
    wlag = study_windows$response_lag[wid]
    # convert to number of observations
    nwin = wlen * 60 / 5
    nlag = wlag * 60 / 5
    # analyze datasets
    cbind(
      # simulation details
      rep = rep_ind,
      window_len = wlen,
      lag_len = wlag,
      # paired simulation results
      rbind(
        cbind(
          sim = 'baseline',
          prepost(x = xb, nwin = nwin, nlag = nlag, baseline_end = nbaseline)
        ),
        cbind(
          sim = 'response',
          prepost(x = xr, nwin = nwin, nlag = nlag, baseline_end = nbaseline)
        )
      )
    )
  }))
}))


#
# save output
#

f = file.path('output', 'sim', 'samples')
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

fname = file.path(f, paste('samples_', taskId, '.rds', sep = ''))

saveRDS(sim_res, file = fname)


#
# plot a snippet of a simulated trajectory
#

if(taskId == 1) {
  
  library(ggplot2)
  library(ggthemes)
  
  # simulate baseline and response data
  seed = .GlobalEnv$.Random.seed
  xr = sim(n = nbaseline, nexposed = nexposed, response = TRUE)
  .GlobalEnv$.Random.seed = seed
  xb = sim(n = nbaseline, nexposed = nexposed, response = FALSE)
  
  # points to include in plot
  inds = seq(to = length(xr), length.out = 100)
  
  # merge data into plottable format
  df = rbind(
    data.frame(depth = xr[inds], t = inds, type = 'Response/Test simulation'),
    data.frame(depth = xb[inds], t = inds, type = 'Baseline/Control simulation')
  )
  
  # plot series
  pl = ggplot(df, aes(x = t, y = depth)) + 
    facet_wrap(~type, nrow = 2) + 
    geom_line() +
    geom_point() + 
    geom_hline(yintercept = 0, lty = 3) + 
    geom_vline(xintercept = nbaseline, lty = 2, alpha = .5) +
    xlab('Simulated observation num.') + 
    scale_y_reverse('Depth (m)') + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    theme_few()
  
  # save plot
  f = file.path('output', 'figures', 'simulation')
  dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
  ggsave(pl, file = file.path(f, 'example_simulation_pair.pdf'), width = 8)
  
  lines(-tail(xr,100),type='l',col=2)
  plot(-tail(xb,100),type='l')
  
}