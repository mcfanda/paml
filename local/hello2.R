library(gamlj)
library(lmerTest)
data("beers_bars")
data<-beers_bars
data$bar<-factor(data$bar)
model<-lmer(smile~beer+(1+beer|bar),data=data)
summary(model)

amod<-simr::extend(model,along = "bar", n=65)
modelpower(amod,effects = "random", nsim=1000)

q<-required_n(model,expand = "bar",direction = "between",effects = "random", nsim=1000)
q
q1<-required_n(model,expand = "bar",direction = "between",effects = c("fixed"))

q1


