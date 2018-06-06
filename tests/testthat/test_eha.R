library(NetworkInference)
library(spid)

context("Event history tuning")

test_that("evaluate grid point produces correct BIC", {
    data(cascades)
    for(c in 1:length(cascades$cascade_times)){
        cascades$cascade_times[[c]] <- round(10*cascades$cascade_times[[c]])
    }
    out = evaluate_grid_point(n_edges = 20, params = 1, time_window = 100,
                              cascades = cascades, decay_parameter = 1)
    expect_equal(round(out, 2), 1056.97)
})

