
production_model <- function(data) {
  # Convert input to a standard data frame
  df <- as.data.frame(data)

  # -------------------------------------------------------------------
  # 1) Create all derived columns exactly as in final model
  # -------------------------------------------------------------------
  df$x1_squared                  <- df$x1^2
  df$x2_squared                  <- df$x2^2
  df$x1_x_x2_squared             <- df$x1 * (df$x2^2)
  df$x1                          <- df$x1
  df$x2_x_x3                     <- df$x2 * df$x3
  df$x2_x_x3_squared             <- df$x2 * (df$x3^2)
  df$x3_squared                  <- df$x3^2
  df$x3_x_x3_squared             <- df$x3 * (df$x3^2)    # x3^3
  df$x3                          <- df$x3
  df$x4_x_x2_squared             <- df$x4 * (df$x2^2)
  df$x1_x_x5                     <- df$x1 * df$x5
  df$x1_x_x4                     <- df$x1 * df$x4
  df$x1_squared_x_x2_squared     <- (df$x1^2) * (df$x2^2)
  df$x1_squared_x_x3_squared     <- (df$x1^2) * (df$x3^2)
  df$x3_x_x1_squared             <- df$x3 * (df$x1^2)
  df$x5_x_x2_squared             <- df$x5 * (df$x2^2)
  df$x5                          <- df$x5
  df$x2_x_x4_squared             <- df$x2 * (df$x4^2)
  df$x4                          <- df$x4
  df$x4_x_x5                     <- df$x4 * df$x5
  df$x5_squared                  <- df$x5^2
  df$x1_x_x4_squared             <- df$x1 * (df$x4^2)
  df$x1_squared_x_x5_squared     <- (df$x1^2) * (df$x5^2)
  df$x4_x_x4_squared             <- df$x4 * (df$x4^2)    # x4^3
  df$x5_x_x4_squared             <- df$x5 * (df$x4^2)
  df$x3_squared_x_x4_squared     <- (df$x3^2) * (df$x4^2)
  df$x3_x_x4_squared             <- df$x3 * (df$x4^2)
  df$x2_x_x5_squared             <- df$x2 * (df$x5^2)
  df$x4_x_x5_squared             <- df$x4 * (df$x5^2)
  df$x4_squared_x_x5_squared     <- (df$x4^2) * (df$x5^2)
  df$x3_squared_x_x5_squared     <- (df$x3^2) * (df$x5^2)
  df$x3_x_x5_squared             <- df$x3 * (df$x5^2)
  df$x2_squared_x_x5_squared     <- (df$x2^2) * (df$x5^2)
  df$x1_squared_x_x4_squared     <- (df$x1^2) * (df$x4^2)
  df$x4_squared                  <- df$x4^2
  df$x5_x_x5_squared             <- df$x5 * (df$x5^2)    # x5^3
  df$x2_x_x1_squared             <- df$x2 * (df$x1^2)
  df$x4_x_x1_squared             <- df$x4 * (df$x1^2)
  df$x3_x_x5                     <- df$x3 * df$x5
  df$x5_x_x3_squared             <- df$x5 * (df$x3^2)
  df$x2_squared_x_x4_squared     <- (df$x2^2) * (df$x4^2)
  df$x4_x_x3_squared             <- df$x4 * (df$x3^2)
  df$x3_x_x4                     <- df$x3 * df$x4
  df$x2_x_x5                     <- df$x2 * df$x5
  df$x2_x_x4                     <- df$x2 * df$x4
  df$x1_x_x5_squared             <- df$x1 * (df$x5^2)
  df$x1_x_x3                     <- df$x1 * df$x3
  df$x1_x_x3_squared             <- df$x1 * (df$x3^2)
  df$x3_x_x2_squared             <- df$x3 * (df$x2^2)
  df$x2_squared_x_x3_squared     <- (df$x2^2) * (df$x3^2)
  df$x2_x_x2_squared             <- df$x2 * (df$x2^2)    # x2^3
  df$x1_x_x2                     <- df$x1 * df$x2
  df$x2                          <- df$x2
  df$x5_x_x1_squared             <- df$x5 * (df$x1^2)
  df$x1_x_x1_squared             <- df$x1 * (df$x1^2)    # x1^3

  # -------------------------------------------------------------------
  # 2) Build the design matrix in the EXACT order matching coefficient list
  # -------------------------------------------------------------------
  X_mat <- cbind(
    df$x1_squared,
    df$x2_squared,
    df$x1_x_x2_squared,
    df$x1,
    df$x2_x_x3,
    df$x2_x_x3_squared,
    df$x3_squared,
    df$x3_x_x3_squared,
    df$x3,
    df$x4_x_x2_squared,
    df$x1_x_x5,
    df$x1_x_x4,
    df$x1_squared_x_x2_squared,
    df$x1_squared_x_x3_squared,
    df$x3_x_x1_squared,
    df$x5_x_x2_squared,
    df$x5,
    df$x2_x_x4_squared,
    df$x4,
    df$x4_x_x5,
    df$x5_squared,
    df$x1_x_x4_squared,
    df$x1_squared_x_x5_squared,
    df$x4_x_x4_squared,
    df$x5_x_x4_squared,
    df$x3_squared_x_x4_squared,
    df$x3_x_x4_squared,
    df$x2_x_x5_squared,
    df$x4_x_x5_squared,
    df$x4_squared_x_x5_squared,
    df$x3_squared_x_x5_squared,
    df$x3_x_x5_squared,
    df$x2_squared_x_x5_squared,
    df$x1_squared_x_x4_squared,
    df$x4_squared,
    df$x5_x_x5_squared,
    df$x2_x_x1_squared,
    df$x4_x_x1_squared,
    df$x3_x_x5,
    df$x5_x_x3_squared,
    df$x2_squared_x_x4_squared,
    df$x4_x_x3_squared,
    df$x3_x_x4,
    df$x2_x_x5,
    df$x2_x_x4,
    df$x1_x_x5_squared,
    df$x1_x_x3,
    df$x1_x_x3_squared,
    df$x3_x_x2_squared,
    df$x2_squared_x_x3_squared,
    df$x2_x_x2_squared,
    df$x1_x_x2,
    df$x2,
    df$x5_x_x1_squared,
    df$x1_x_x1_squared
  )

  # -------------------------------------------------------------------
  # 3) Define the coefficients (in the same order) and intercept
  # -------------------------------------------------------------------
  coefs <- c(
    8.111226,     # x1_squared
    8.044411,     # x2_squared
    4.133054,     # x1_x_x2_squared
    3.317359,     # x1
    2.245588,     # x2_x_x3
    2.245588,     # x2_x_x3_squared
    1.795218,     # x3_squared
    1.795218,     # x3_x_x3_squared
    1.795218,     # x3
    1.298972,     # x4_x_x2_squared
    1.077665,     # x1_x_x5
    0.713558,     # x1_x_x4
    0.679454,     # x1_squared_x_x2_squared
    0.507564,     # x1_squared_x_x3_squared
    0.507564,     # x3_x_x1_squared
    0.312359,     # x5_x_x2_squared
    0.277009,     # x5
    0.236092,     # x2_x_x4_squared
    0.171768,     # x4
    0.057962,     # x4_x_x5
    0.034288,     # x5_squared
    0.032483,     # x1_x_x4_squared
    0.006601,     # x1_squared_x_x5_squared
    0.005040,     # x4_x_x4_squared
    0.002527,     # x5_x_x4_squared
    0.002464,     # x3_squared_x_x4_squared
    0.002464,     # x3_x_x4_squared
    0.001668,     # x2_x_x5_squared
    0.000880,     # x4_x_x5_squared
    0.000052,     # x4_squared_x_x5_squared
   -0.004075,     # x3_squared_x_x5_squared
   -0.004075,     # x3_x_x5_squared
   -0.004558,     # x2_squared_x_x5_squared
   -0.010061,     # x1_squared_x_x4_squared
   -0.010395,     # x4_squared
   -0.041476,     # x5_x_x5_squared
   -0.048625,     # x2_x_x1_squared
   -0.082624,     # x4_x_x1_squared
   -0.146048,     # x3_x_x5
   -0.146048,     # x5_x_x3_squared
   -0.181626,     # x2_squared_x_x4_squared
   -0.207875,     # x4_x_x3_squared
   -0.207875,     # x3_x_x4
   -0.395717,     # x2_x_x5
   -1.135356,     # x2_x_x4
   -1.462084,     # x1_x_x5_squared
   -1.832480,     # x1_x_x3
   -1.832480,     # x1_x_x3_squared
   -2.884968,     # x3_x_x2_squared
   -2.884968,     # x2_squared_x_x3_squared
   -3.221832,     # x2_x_x2_squared
   -5.140215,     # x1_x_x2
   -5.614728,     # x2
   -17.104797,    # x5_x_x1_squared
   -66.345228     # x1_x_x1_squared
  )

  # Intercept
  intercept <- -0.06228481654693496

  # -------------------------------------------------------------------
  # 4) Compute predictions = intercept + X_mat * coefs
  # -------------------------------------------------------------------
  y_pred <- intercept + X_mat %*% coefs

  # Return predictions as a numeric vector
  return(as.vector(y_pred))
}
