
library(NetworkInference)
data("cascades")

for(c in 1:length(cascades$cascade_times)){
    cascades$cascade_times[[c]] <- round(10*cascades$cascade_times[[c]])
}

params = c(0.5, 1, 2)
time_windows = c(10, 50, 100)
n_edges = c(5, 10, 20)
library(spid)

out = grid_search_eha(cascades, 4, n_edges, params, time_windows)
