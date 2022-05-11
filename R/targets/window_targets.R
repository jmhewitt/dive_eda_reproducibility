window_targets = list(

  # window lengths (h) to analyze
  tar_target(window_lengths, c(6,3,1)),
  
  # define nested pre/post windows to use in analysis
  tar_target(
    name = study_windows,
    command = do.call(
      rbind, lapply(window_lengths, function(window_length, max_length) {
        # response_lags to partition (0, max_length) w/"window_length" intervals
        response_lags = seq(from = 0, to = max_length, by = window_length)
        # intervals that partition (0, max_length) 
        do.call(rbind, lapply(2:length(response_lags), function(i) {
          data.frame(description = paste(response_lags[i-1], '--', 
                                         response_lags[i], 'h', sep =''),
                     response_lag = response_lags[i-1], 
                     window_length = window_length)
        }))
      }, max_length = max(window_lengths))
    )
  ),
  
  # only use baseline windows that match the context of the exposure
  tar_target(conditional, c(TRUE, FALSE)),
  
  # identify baseline and testing pre/post window pairs in sattag records
  tar_target(
    name = pre_post_pair_defs, 
    command = {
      
      # skip processing if tag is not associated with a CEE
      if(!is.finite(tag_info$cee_start)) {
        return(NULL)
      }
      
      # initialize output
      res = list(
        tag = tag_info,
        window = study_windows
      )
      
      # identify depth record file for tag
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
      
      
      #
      # data windows
      #
      
      res$response_window = tag_info$cee_start + lubridate::seconds(
        3600*(study_windows$response_lag + c(0, study_windows$window_length))
      )
      
      res$pre_exposure_window = tag_info$cee_start - lubridate::seconds(
        3600 * c(study_windows$window_length, 0)
      )
      
      res$baseline_window = c(d$Date[1], tag_info$baseline_end)

      res$response_inds = which(
        (res$response_window[1] <= d$Date) &
        (d$Date < res$response_window[2])
      )
      
      res$pre_exposure_inds = which(
        (res$pre_exposure_window[1] <= d$Date) & 
        (d$Date < res$pre_exposure_window[2])
      )
      
      res$baseline_inds = which(
        (res$baseline_window[1] <= d$Date) &
        (d$Date < res$baseline_window[2])
      )
      
      #
      # validate pre-exposure and response windows are compatible
      #
      
      # skip processing if no response observations are observed
      if(length(res$response_inds) == 0) {
        return(NULL)
      }
      
      # skip processing if pre-exposure and response windows are not balanced
      if(length(res$response_inds) != length(res$pre_exposure_inds)) {
        return(NULL)
      }
      
      #
      # determine valid baseline window pairs
      #

      if(conditional) {
        pre_exposure_deep = d$depth.standardized[
          tail(res$pre_exposure_inds,1)
        ] > deep_dive_depth
      }
        
      res$baseline_pairs = lapply(res$baseline_inds, function(ind) {
        
        # 
        # pre/post windows within baseline data
        #
        
        context_window = d$Date[ind] + lubridate::seconds(
          3600 * c(0, study_windows$window_length)
        )
        
        study_window = context_window[2] + lubridate::seconds(
          3600*(study_windows$response_lag + c(0, study_windows$window_length))
        )
        
        context_inds = which(
          # times within context window
          (context_window[1] <= d$Date) &
          (d$Date < context_window[2]) &
          # that are restricted to lie within baseline period
          (res$baseline_window[1] <= d$Date) &
          (d$Date < res$baseline_window[2])
        )
        
        study_inds = which(
          # times within study window
          (study_window[1] <= d$Date) &
          (d$Date < study_window[2]) &
          # that are restricted to lie within baseline period
          (res$baseline_window[1] <= d$Date) &
          (d$Date < res$baseline_window[2])
        )
        
        #
        # validate that baseline windows are comparable to exposure windows
        #
        
        if(length(context_inds) != length(res$pre_exposure_inds)) {
          return(NULL)
        }
        
        if(length(study_inds) != length(res$response_inds)) {
          return(NULL)
        }
        
        #
        # validate that dive context/state of baseline windows matches exposure
        #
        
        if(conditional) {
          context_deep = d$depth.standardized[
            tail(context_inds,1)
          ] > deep_dive_depth
          
          if(!identical(pre_exposure_deep, context_deep)) {
            return(NULL)
          }
        }
        
        # return indices and windows
        list(
          context_inds = context_inds,
          study_inds = study_inds,
          context_window = context_window,
          study_window = study_window
        )
        
      })
            
      # remove nulls from the output
      res$baseline_pairs = res$baseline_pairs[
        !sapply(res$baseline_pairs, is.null)
      ]
      
      # for creating a complete time series sequence of statistics
      res$all_pairs = lapply(1:length(d$Date), function(ind) {

        #
        # pre/post windows relative to data index
        #

        context_window = d$Date[ind] + lubridate::seconds(
          3600 * c(0, study_windows$window_length)
        )

        study_window = context_window[2] + lubridate::seconds(
          3600*(study_windows$response_lag + c(0, study_windows$window_length))
        )

        context_inds = which(
          # times within context window
          (context_window[1] <= d$Date) &
          (d$Date < context_window[2])
        )

        study_inds = which(
          # times within study window
          (study_window[1] <= d$Date) &
          (d$Date < study_window[2])
        )

        #
        # validate that baseline windows are comparable to exposure windows
        #

        if(length(context_inds) != length(res$pre_exposure_inds)) {
          return(NULL)
        }

        if(length(study_inds) != length(res$response_inds)) {
          return(NULL)
        }

        #
        # validate that dive context/state of baseline windows matches exposure
        #

        if(conditional) {
          context_deep = d$depth.standardized[
            tail(context_inds,1)
          ] > deep_dive_depth

          if(!identical(pre_exposure_deep, context_deep)) {
            return(NULL)
          }
        }

        # return indices and windows
        list(
          context_inds = context_inds,
          study_inds = study_inds,
          context_window = context_window,
          study_window = study_window
        )

      })
      
      # remove nulls from the output
      res$all_pairs = res$all_pairs[
        !sapply(res$all_pairs, is.null)
      ]
      
      # skip processing if no baseline data is valid
      if(length(res$baseline_pairs) == 0) {
        return(NULL)
      }
      
      res$baseline_pairs_are_conditional = conditional
      
      list(res)
    },
    pattern = cross(map(tag_info), map(study_windows), map(conditional)),
    deployment = 'worker', 
    memory = 'transient'
  )
  
)
