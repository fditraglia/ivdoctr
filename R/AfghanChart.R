# Generates outer envelope for Afghan data
# buildschool is Instrument
# enrolled07 is treatment
# enrolled is outcome
library(plotly)

p <- mean(afghan$buildschool)
r_Ty <- cor(afghan$enrolled07, afghan$enrolled)
r_Tz <- cor(afghan$enrolled, afghan$buildschool)
r_zy <- cor(afghan$buildschool, afghan$enrolled)
s2_T <- var(afghan$enrolled07)
draws <- data.frame(r_Ty = r_Ty, r_Tz = r_Tz, r_zy = r_zy)
L <- get_L(draws)
max_a0 <- s2_T * (1 - L) / (1 - p)
a0 <- seq(0, max_a0, by = 0.0001)
a1 <- (s2_T * (1 - L) - (1 - p) * a0) / (p - a0)
plot_ly(x = ~a0) %>%
  add_lines(y = ~a1) %>%
  layout(xaxis = list(range = c(0, p)),
         yaxis = list(range = c(0, 1 - p)))
