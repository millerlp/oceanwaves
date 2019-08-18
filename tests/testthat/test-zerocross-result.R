context("Testing zero-cross results")

test_that("Returns list of 6 numeric values", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(length(res), 6)
  for (i in 1:length(res)){
    expect_type(res[[i]],'double')
  }
})

test_that("Returns 12 second mean period", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(res$Tmean, 12.01042, tolerance = 0.001)
})

test_that("Returns 12 second sig. period", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(res$Tsig, 12)
})

test_that("Returns ~3.002xx Hsig.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(res$Hsig, 3.002, tolerance = 0.0001)
})


test_that("Returns 10.75 Tsig.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 5 second period and 12 second cycle 
  y =  sin(2*pi*(1/5)*x) + 1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(res$Tsig, 10.75, tolerance = 0.00001)
})

test_that("Returns 4.74 Hsig.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 5 second period and 12 second cycle 
  y =  sin(2*pi*(1/5)*x) + 1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsZC(y,Fs = 4)
  expect_equal(res$Hsig, 4.74, tolerance = 0.001)
})

