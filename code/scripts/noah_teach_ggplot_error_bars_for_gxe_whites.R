library(data.table)
library(ggplot2)


r1 <- c(1,2,3,1,2,1)
r2 <- c(1,2,1,2,1,2)
r3 <- c(1,1,1,2,2,3)

r1 model
r2 mean
r3 VC
rdata <- data.table(r1, r2, r3)

r1 <- c(1,2,3,1,2,1)
r2 <- c(1,2,1,2,1,2)
r3 <- c(1,1,1,2,2,3)
rdata <- data.table(r1, r2, r3)
ggplot(rdata, aes(color = r3, x = r1, y= r2)) +
  geom_point(stat = "identity",
           position = position_dodge2(width = 1, preserve = "single"))

ggplot(df_long_ci, aes(fill = VC, x = Model, y= mean)) +
  geom_point(colour = "black", stat = "identity",
           position = position_dodge2(width=1,preserve = "single"))+
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 1,
                position = position_dodge2(width = 1, preserve = "single"), linewidth = 0.2)
