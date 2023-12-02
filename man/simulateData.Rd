\name{simulateData}
\alias{simulateData}
\title{simulateData}
\usage{
simulateData()
}
\description{
Generate simulation data for test
}
\examples{
# ## Examples
#
# inputdata_example <- data.frame(
#   time = c(0, 8, 8, 17, 24),
#   n = c(0, 15, 2, 0, 0),
#   age = c(0, 8, 50, 0, 0),
#   gender = rep(NA, 5),
#   MET = c(NA, 3, 1.8, NA, NA),
#   CO2rate = rep(NA, 5)
# )
#
# exdata <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01)
# exdata2 <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01)
#
# library(ggplot2)
# library(patchwork)
# plot1 <- ggplot(exdata, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot2 <- ggplot(exdata2, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot1 + plot2 + plot_annotation(title='Euler method')
#
# # compare to using Exponential method
# exdata_exp <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01, method='Exponential')
# exdata2_exp <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01, method='Exponential')
# plot1_exp <- ggplot(exdata_exp, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot2_exp <- ggplot(exdata2_exp, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot1_exp + plot2_exp + plot_annotation(title='Exponential method')
#
#
# inputdata_example_staggered_entrance <- data.frame(
#   time = c(0, 8, 8, 10,10, 12,12, 17, 24),
#   n = c(0, 5, 2, 10, 2, 15, 2, 0, 0),
#   age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
#   gender = rep(NA, 9),
#   MET = c(NA, 3, 1.8, 3, 1.8, 3, 1.8, NA, NA),
#   CO2rate = rep(NA, 9)
# )
# exdata3 <- simulateData(persondata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .1, CO2var = 0.01)
# exdata4 <- simulateData(persondata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .5, CO2var = 0.01)
#
# plot3 <- ggplot(exdata3, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot4 <- ggplot(exdata4, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot3 + plot4
#
#
# inputdata_example_moving <- data.frame(
#   time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
#   n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
#   age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
#   gender = rep(NA, 9),
#   MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
#   CO2rate = rep(NA, 9)
# )
# exdata5 <- simulateData(persondata = inputdata_example_moving, volume=500, ventilation_rate = .1, CO2var = 0.01)
# ggplot(exdata5, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; go in and out') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
}
