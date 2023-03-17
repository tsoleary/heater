# Test script to try the scheduler out

require(tidyverse)

p <- ggplot(iris) +
  geom_point(aes(x = Sepal.Length,
                 y = Sepal.Width))

ggsave("testplot.pdf")
