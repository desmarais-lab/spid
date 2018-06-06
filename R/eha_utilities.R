#' Create a policy adoption event history data frame
#'
#' This the short description
#'
#' And some details
#'
#' @import dplyr
#'
#' @param  cascades an object of class cascade containing node and cascade
#'     information. See \code{\link{as_cascade_long}} and
#'     \code{\link{as_cascade_wide}} for details.
#'
#' @return Returns an object of class \code{data.frame}
#'
#' @examples
#'
#' data(cascades)
#'
#' @export
make_eha_data = function(cascades, networks, decay_parameter, min_time) {

    # Create the time x cascade_id grid (withing the time range of each cascade)
    ranges = lapply(1:length(cascades$cascade_times), function(i) {
        x = cascades$cascade_times[[i]]
        data_frame(event_time = min(x):max(x),
                   cascade_id = names(cascades$cascade_times[i]) )
    })
    ranges = do.call(rbind, ranges)

    # Create event indicator
    p <- tbl_df(as.data.frame(cascades)) %>%
        arrange(cascade_id, event_time) %>%
        mutate(event = 1)

    # Get for time t number the of adoptions until (including) t - 1
    events_by_time = group_by(p, cascade_id, event_time) %>%
        summarize(count = n()) %>%
        right_join(ranges, by = c('event_time', 'cascade_id')) %>%
        group_by(cascade_id) %>%
        mutate(count = ifelse(is.na(count), 0, count),
               events_so_far = c(0, cumsum(count)[-length(count)])) %>%
        dplyr::select(-count)

    # All node-time combinations for each cascade
    ranges_list = split(ranges, f = ranges$cascade_id)
    node_time_combos = lapply(ranges_list, function(x) {
        tbl_df(expand.grid(cascades$node_names, x$event_time,
                           stringsAsFactors = FALSE)) %>%
            rename(node_name = Var1, event_time = Var2) %>%
            mutate(cascade_id = x$cascade_id[1])
    })

    node_time_combos = do.call(rbind, node_time_combos) %>%
        # Fill in the event outcome
        left_join(p, by = c('node_name', 'event_time', 'cascade_id')) %>%
        # Fill outcome of non-matched node-times with zero
        mutate(event = ifelse(is.na(event), 0, 1)) %>%
        # For each time get how many events happened in previous times
        left_join(events_by_time, by = c('cascade_id', 'event_time')) %>%
        # Generate cumulative indicator if node had event at any preceding time
        arrange(cascade_id, node_name, event_time) %>%
        group_by(cascade_id, node_name) %>%
        mutate(node_cumsum = cumsum(event),
               node_cumsum = ifelse(event == 1, 0, node_cumsum),
               decay_weight = decay_event(event, decay_parameter)) %>%
        ungroup()

    # Calculate number of events of neighbors for each cascade-node-time
    neighbor_events = left_join(networks, node_time_combos,
                                by = c("origin_node" = "node_name",
                                       "time" = "event_time")) %>%
        filter(!is.na(node_cumsum)) %>%
        dplyr::select(-event, -events_so_far) %>%
        group_by(cascade_id, destination_node, time) %>%
        summarize(n_neighbor_events = sum(node_cumsum),
                  n_neighbor_events_decay = sum(decay_weight))

    eha_data = filter(node_time_combos, event_time >= min_time, node_cumsum == 0) %>%
        left_join(neighbor_events, by = c('cascade_id' = 'cascade_id',
                                          'node_name' = 'destination_node',
                                          'event_time' = 'time')) %>%
        dplyr::select(-node_cumsum, -decay_weight) %>%
        mutate(n_neighbor_events = ifelse(is.na(n_neighbor_events), 0,
                                          n_neighbor_events),
               n_neighbor_events_decay = ifelse(is.na(n_neighbor_events_decay),
                                                0, n_neighbor_events_decay),
               cascade_id = as.factor(cascade_id))
    return(eha_data)
}

decay_event = function(x, lambda) {
    if(sum(x) == 0) return(x)
    idx = which(x == 1)
    y = seq(0,(length(x)-idx))
    dec_seq = lambda * exp(-y/lambda)
    # Make the weight 0 in the actual adoption year and shift the sequence
    dec_seq = c(0, dec_seq[-length(dec_seq)])
    return(c(rep(0, idx-1), dec_seq))
}


#' Evaluate a parameter combination
#'
#' This the short description
#'
#' And some details
#'
#' @import speedglm
#'
#' @param  cascades an object of class cascade containing node and cascade
#'     information. See \code{\link{as_cascade_long}} and
#'     \code{\link{as_cascade_wide}} for details.
#'
#' @return Returns an object of class \code{data.frame}
#'
#' @examples
#'
#' data(cascades)
#'
#' @export
evaluate_grid_point <- function(n_edges, params, time_window, cascades,
                                min_time = NULL, decay_parameter) {

    max_time = max(unlist(cascades$cascade_times))
    if(is.null(min_time)) min_time = time_window + 1
    times = seq(min_time, max_time, 1)

    # Infer the network for all years
    networks = do.call(rbind, lapply(times, infer_network,
                                     time_window = time_window,
                                     cascades = cascades, params = params,
                                     n_edges = n_edges))
    if(nrow(networks) == 0) return(NA)

    # Build dataset for event history model
    event_data = make_eha_data(cascades, networks, decay_parameter, min_time)

    # Fit model
    mod <- event ~ events_so_far + n_neighbor_events_decay + cascade_id
    res <- speedglm::speedglm(mod, data = droplevels(event_data),
                              family = binomial(link = "logit"))
    return(BIC(res))
}



#' Evaluate a parameter grid
#'
#' This the short description
#'
#' And some details
#'
#' @import doParallel
#'
#' @param  cascades an object of class cascade containing node and cascade
#'     information. See \code{\link{as_cascade_long}} and
#'     \code{\link{as_cascade_wide}} for details.
#'
#' @return Returns an object of class \code{data.frame}
#'
#' @examples
#'
#' data(cascades)
#'
#' @export
grid_search_eha = function(cascades, n_jobs, n_edges, params, time_windows) {

    grid_points <- expand.grid(params, n_edges, time_windows)

    # Run
    cl <- makeCluster(n_jobs)
    registerDoParallel(cl)

    results <- foreach(i = 1:nrow(grid_points),
                       .packages = c("dplyr", "NetworkInference", "spid"),
                       .combine = c) %dopar% {
                           evaluate_grid_point(n_edges = grid_points[i, 2],
                                               params = grid_points[i, 1],
                                               time_window = grid_points[i, 3],
                                               cascades = cascades,
                                               decay_parameter = grid_points[i, 1])
                       }
    stopCluster(cl)

    out = as.data.frame(cbind(grid_points, results))
    colnames(out) = c('params', 'n_edges', 'time_window', 'bic')
    return(out)
}

#' Infer diffusion network on time slice of cascades
#'
#' This the short description
#'
#' And some details
#'
#@import checkmate
#'
#' @param  cascades an object of class cascade containing node and cascade
#'     information. See \code{\link{as_cascade_long}} and
#'     \code{\link{as_cascade_wide}} for details.
#'
#' @return Returns an object of class \code{data.frame}
#'
#' @examples
#'
#' data(cascades)
#'
#' @export
infer_network <- function(time, time_window, cascades, params, n_edges) {

    casc = subset_cascade_time(cascades, (time - time_window), time)

    network = try(netinf(cascades = casc, trans_mod = "exponential",
                         n_edges = n_edges, params = params, quiet = TRUE),
                  silent = TRUE)

    if(inherits(network, 'try-error')){
        warning(paste0('Could not fit netinf for parameter combination:',
                       ' n_edges: ', n_edges, ', params: ', params,
                       ', time: ', time, ', time_window: ', time_window,
                       '. Reason: ', attr(network, 'condition'),
                       '. Returning empty network'))
        return(data.frame(origin_node = c(), destination_node = c(),
                          year = c()))
    }

    network = network[, -c(3, 4)]
    network$time = time
    return(network)
}
