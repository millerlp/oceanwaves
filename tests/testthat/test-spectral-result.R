context("Testing spectral analysis results")

test_that("Returns list of 8 numeric values", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(length(res), 8)
  for (i in 1:length(res)){
    expect_type(res[[i]],'double')
  }
})

test_that("Returns 12 second peak period", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(res$Tp, 12)
})

test_that("Returns expected Hm0", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(res$Hm0, 4.23, tolerance = 0.0001)
})

test_that("Returns approximate expected T_0_1.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 12 second period cycle only
  y =  1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(res$T_0_1, 11.9387, tolerance = 0.0001)
})


test_that("Returns 12 sec peak period on mixed signal.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 5 second period and 12 second cycle 
  y =  sin(2*pi*(1/5)*x) + 1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(res$Tp, 12)
})

test_that("Returns approximate 5.076 Hm0.", {
  # Generate set of x values, at 0.25sec interval
  x = seq(0,(5*60)-0.25,by = 0.25)
  # 5 second period and 12 second cycle 
  y =  sin(2*pi*(1/5)*x) + 1.5*sin(2*pi*(1/12)*(x))
  res = waveStatsSP(y,Fs = 4)
  expect_equal(res$Hm0, 5.076, tolerance = 0.001)
})
