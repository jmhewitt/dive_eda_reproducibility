pre_post_targets = list(

  # compute summary statistics for pre/post pairs at CEE and baseline times
  tar_target(
    name = sattag_eda,
    command = {
      
      # skip processing if no information
      if(is.null(pre_post_pair_defs)) {
        return(NULL)
      }
      
      # identify depth record file for tag
      tag_info = pre_post_pair_defs[[1]]$tag
      f = grep(pattern = tag_info$deployid, x = depth_files, value = TRUE)
      
      # load data
      d = read.csv(f)
      d$Date = as.POSIXct(d$Date, origin = '1970-01-01 00:00.00 UTC', 
                          tz = 'UTC')
      
      # map all depths to standardized bins
      d$depth.bin = sapply(d$Depth, function(depth) {
        which.min(abs(depth - template_bins$center))
      })
      d$depth.standardized = template_bins$center[d$depth.bin]
      
      # add dive segmentation to depths time series
      d$diveId = dive_labels[[
        which(sapply(dive_labels, function(x) x$tag == tag_info$deployid))
      ]]$labels
      
      # classify dives by max observed depth
      diveTypes = d %>%
        dplyr::group_by(diveId) %>%
        dplyr::summarise(
          diveType = ifelse(any(diveId == 0), 'Unknown', ifelse(
            max(depth.standardized) > deep_dive_depth, 'Deep', 'Shallow'
          ))
        )
      
      # merge dive type information with sattag record
      d = d %>% left_join(diveTypes, by = 'diveId')
      
      # reformat column names and add derived data to sattag record
      d = d %>% mutate(
        depths = depth.standardized,
        ddepths = c(0, diff(depths)),
        ddepths.sign = sign(ddepths),
        descending = ddepths.sign > 0,
        ascending = ddepths.sign < 0,
        no_change = ddepths.sign == 0
      )
      
      #
      # process baseline data
      #
      
      # difference-based null distribution samples
      null_samples = do.call(rbind, lapply(
        pre_post_pair_defs[[1]]$baseline_pairs, function(x) {
          sattag_summary(d[x$study_inds,]) - sattag_summary(d[x$context_inds,])
        }
      ))
      
      # bivariate null distribution samples
      null_samples_pre = do.call(rbind, lapply(
        pre_post_pair_defs[[1]]$baseline_pairs, function(x) {
          sattag_summary(d[x$context_inds,])
        }
      ))
      null_samples_post = do.call(rbind, lapply(
        pre_post_pair_defs[[1]]$baseline_pairs, function(x) {
          sattag_summary(d[x$study_inds,])
        }
      ))
      
      # observed value of difference-based test statistics
      observed = sattag_summary(d[pre_post_pair_defs[[1]]$response_inds,]) - 
        sattag_summary(d[pre_post_pair_defs[[1]]$pre_exposure_inds,])
      
      # observed value of bivariate test statistics
      observed_pre = sattag_summary(d[pre_post_pair_defs[[1]]$pre_exposure_inds,])
      observed_post = sattag_summary(d[pre_post_pair_defs[[1]]$response_inds,])
        
      # bivariate p-values
      p.bivariate = do.call(rbind, lapply(1:ncol(null_samples), function(ind) {
        
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
        
        # mahalanobis distance calculation components
        df = cbind(null_samples_pre[,ind], null_samples_post[,ind])
        md = mahalanobis(
          x = rbind(df, cbind(observed_pre[[ind]], observed_post[[ind]])), 
          center = colMeans(df), 
          cov = var(df)
        )
        
        # bivariate p-values
        data.frame(
          p.left = c(cdf.conditional[kde_ind_post], NA),
          p.right = c(ccdf.conditional[kde_ind_post],
                      mean(tail(md, 1) <= md[1:nrow(df)])),
          stat = factor(colnames(null_samples)[ind]),
          type = factor(c('bivariate_kde', 'bivariate_mdist'))
        )
      }))
      
      # right-tail probability
      p.right = rowMeans(apply(X = null_samples, MARGIN = 1, FUN = function(r) {
        r >= observed
      }))
      
      # left-tail probability
      p.left = rowMeans(apply(X = null_samples, MARGIN = 1, FUN = function(r) {
        r <= observed
      }))

      # package results
      list(list(
        tag = pre_post_pair_defs[[1]]$tag,
        window = pre_post_pair_defs[[1]]$window, 
        response_window = pre_post_pair_defs[[1]]$response_window,
        pre_exposure_window = pre_post_pair_defs[[1]]$pre_exposure_window, 
        baseline_window = pre_post_pair_defs[[1]]$baseline_window,
        conditional = pre_post_pair_defs[[1]]$baseline_pairs_are_conditional,
        null_samples = null_samples,
        null_samples_pre = null_samples_pre,
        null_samples_post = null_samples_post,
        observed = observed,
        observed_pre = observed_pre,
        observed_post = observed_post,
        p.left = p.left,
        p.right = p.right,
        p.bivariate = p.bivariate
      ))
    }, 
    pattern = map(pre_post_pair_defs), 
    deployment = 'worker'
  )
  
)
