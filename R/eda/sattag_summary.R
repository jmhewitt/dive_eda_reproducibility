sattag_summary = function(dat) {
  # dat %>% dplyr::summarise(
  #   'prop_surface' = mean(depth.bin == 1),
  #   'prop_downward' = mean(descending),
  #   'prop_upward' = mean(ascending),
  #   'prop_no_change' = mean(no_change),
  #   'total_absolute' = sum(abs(ddepths)),
  #   'total_downward' = sum(abs(ddepths) * descending),
  #   'total_upward' = sum(abs(ddepths) * ascending),
  #   'depth_seq' = sum(1:n() * depths),
  #   'type_seq' = sum(1:n() * (diveType == 'Deep')),
  #   'iddi_obs' = mean(diveType == 'Shallow'),
  #   # sequences of differences
  #   'depth_diff_seq' = sum(1:(n()-1) * ddepths[-1]),
  #   'direction_change_seq' = sum(1:(n()-1) * (abs(diff(ddepths.sign)) > 0))
  # ) 
  dat %>% dplyr::summarise(
    'time_on_surface' = sum(depth.bin == 1),
    'avg_depth' = mean(depths),
    'time_deep_diving' = sum(diveType == 'Deep'),
    'time_ascending' = sum(ascending),
    'time_descending' = sum(descending),
    'time_level' = sum(no_change),
    'total_vertical' = sum(abs(ddepths)),
    'total_ascending' = sum(abs(ddepths) * ascending),
    'total_descending' = sum(abs(ddepths) * descending),
    'direction_changes' = sum(abs(diff(ddepths.sign)) > 0),
    
    'timing_time_on_surface' = sum((1:n()) * (depth.bin == 1)),
    'timing_avg_depth' = mean((1:n()) * depths),
    'timing_time_deep_diving' = sum((1:n()) * (diveType == 'Deep')),
    'timing_time_ascending' = sum((1:n()) * ascending),
    'timing_time_descending' = sum((1:n()) * descending),
    'timing_time_level' = sum((1:n()) * no_change),
    'timing_total_vertical' = sum((1:n()) * abs(ddepths)),
    'timing_total_ascending' = sum((1:n()) * (abs(ddepths) * ascending)),
    'timing_total_descending' = sum((1:n()) * (abs(ddepths) * descending)),
    'timing_direction_changes' = sum((1:(n()-1)) * (abs(diff(ddepths.sign)) > 0))
  ) 
}