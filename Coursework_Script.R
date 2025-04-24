library(ggplot2)
library(tidyverse)
library(tune)
#PART 1

experiment = read_csv("experiment.csv")
head(experiment)

# a)
ggplot(data=experiment, aes(x=HRR_delta180, y=MAS)) +
  geom_point() +
  labs(x="HRR_delta180 (bpm)", y="MAS (kph)")
  
ggplot(data=experiment) +
  geom_histogram(aes(x=HRR_delta180, Y=..density..)) +
  geom_vline(xintercept=mean(experiment$HRR_delta180), colour="red") +
  labs(x="HRR_delta180 (bpm)", y="Density")


ggplot(data=experiment) +
  geom_histogram(aes(x=MAS, y=..density..)) +
  geom_vline(xintercept=mean(experiment$MAS), colour="red") +
  labs(x="MAS (kph)", y="Density")


# b)
simpleModel = lm(MAS ~ HRR_delta180, data=experiment)
summary(simpleModel)

ggplot(data=experiment, aes(x=HRR_delta180, y=MAS)) +
  geom_point() +
  stat_smooth(method="lm") +
  stat_regline_equation(label.x.npc = 0, label.y.npc = 1) +
  labs(x="HRR_delta180 (bpm)", y="MAS (kph)")


summary(aov(simpleModel))

extendedModel = lm(MAS ~ BMI + Age + Sex + HRR_delta180, data=experiment)
summary(extendedModel)

summary(aov(extendedModel))


# c)
extendedModel2 = lm(MAS ~ (BMI + Age + HRR_delta180)*Sex, data=experiment)
summary(extendedModel2)

summary(aov(extendedModel2))


# d)
# look at assumptions of true residuals in linear model

# standardized residuals vs fitted values
ggplot (extendedModel2, mapping=aes(x=.fitted, y=.stdresid)) +
  geom_point() +
  geom_smooth(method="loess", se=FALSE) +
  geom_hline(yintercept=0, linetyp="dotted") +
  labs(x="Fitted values", y="Standardised residuals")

# justifying suggesting of outliers:

# q-q plot to compare standardised residuals to the expected standard normal model
ggplot(extendedModel2, mapping= aes(sample=.stdresid)) +
  geom_qq_line(colour="red", linewidth=1) +
  geom_qq() +
  labs(x="Theoretical quantiles", y="Standardised residuals") +
  tune::coord_obs_pred()

# histogram to compare with standard normal curve to show this
ggplot(extendedModel2) +
  geom_histogram(aes(x = .stdresid, y = ..density..)) +
  stat_function(geom="line", fun=dnorm, colour="red", linewidth=1) +
  labs(x="Standardised Residuals", y="Density")

# finding those outliers

# standardised residuals vs leverage
cd_cont_pos <- function (leverage, level, model)
{sqrt (level*length (coef (model))* (1- leverage)/ leverage)}

cd_cont_neg <- function (leverage, level, model)
{-cd_cont_pos (leverage, level, model)}

level <- 4/nrow(extendedModel2$model)

ggplot (data=extendedModel2, mapping=aes(x=.hat, y=.stdresid)) +
  geom_point() +
  geom_smooth(method="loess", se=FALSE) +
  stat_function(fun=cd_cont_pos, args=list(level=level, model=extendedModel2), colour="red") +
  stat_function(fun=cd_cont_neg, args=list(level=level, model=extendedModel2), colour="red") +
  geom_hline(yintercept=0, linetype="dotted") +
  labs(x="Leverage", y="Standardised residuals")


# e)
female_data = experiment[experiment$Sex == 1,]
head(female_data)

# in order of increasing p-value for HRR index:

HRR_delta180Model = lm(MAS ~ Age + BMI + HRR_delta180, data=female_data)
summary(HRR_delta180Model)

HRR_delta180_holm = min((4-1+1)*2.7e-7, 1)

HRR180perModel = lm(MAS ~ Age + BMI + HRR180per, data=female_data)
summary(HRR180perModel)

HRR180per_holm = min((4-2+1)*3.48e-7, 1)

HRR10perModel = lm(MAS ~ Age + BMI + HRR10per, data=female_data)
summary(HRR10perModel)

HRR10per_holm = min((4-3+1)*0.0135, 1)

HRR_delta10Model = lm(MAS ~ Age + BMI + HRR_delta10, data=female_data)
summary(HRR_delta10Model)

HRR_delta10_holm = min((4-4+1)*0.0157, 1)

# bonferroni correction to 0.05 sig lvl
bonferroni_sig_level = 0.05/4



# PART 2

drake = read_csv("drake_music.csv")
head(drake)

kendrick = read_csv("kendrick_music.csv")
head(kendrick)

# a)
loglik = function(p, data) {
  a = p[1]
  b = p[2]
  n = nrow(data)
  
  #summation term:
  s = sum((a-1)*log(data$danceability) + (b-1)*log(1-data$danceability))
  
  #remaining terms:
  x = (- n*lgamma(a)-n*lgamma(b) 
       + n*lgamma(a+b))
  
  return(x+s)
}

results = optim(c(10, 10), loglik, control=list(fnscale=-1), data=kendrick)

alpha_hat = results$par[1]
beta_hat = results$par[2]

# b)
ggplot(kendrick, mapping = aes(sample=danceability)) +
  geom_abline(color="red", linewidth=1) +
  stat_qq(distribution=qbeta, dparams=list(alpha_hat, beta_hat)) +
  labs(x="Theoretical quantiles", y="Sample Quantiles") +
  xlim(c(0, 1)) +
  ylim(c(0, 1))

ggplot(kendrick) +
  geom_histogram(aes(x = danceability, y = ..density..)) +
  stat_function(geom="line", fun=dbeta,
                args =list(alpha_hat, beta_hat),
                colour="red", linewidth=1) +
  geom_vline(xintercept=mean(kendrick$danceability), colour="blue") +
  labs(x="Danceability", y="Density")

# c)
bins = c(0, 0.15, 0.3, 0.45, 0.6, 0.75, 1)
observed = table(cut(kendrick$danceability, bins))
cdf = pbeta(bins, alpha_hat, beta_hat)
expected = nrow(kendrick)*diff(cdf)

test_stat = sum(((observed-expected)^2)/expected)

dof = length(bins)-1 -2 -1
p_value = 1-pchisq(test_stat, dof)
critical = qchisq(0.95, 3)


# d)
area = function(x1, x2, x3) {
  0.5*abs(x1[1]*(x2[2]-x3[2]) +
            x2[1]*(x3[2]-x1[2]) +
            x3[1]*(x1[2]-x2[2]))
}

iter = 1000
kendrick_areas = numeric(iter)
drake_areas = numeric(iter)

for(i in 1:iter) {
  kendrick_sample = kendrick[sample(nrow(kendrick), 3), ]
  drake_sample = drake[sample(nrow(drake), 3), ]
  
  kendrick_areas[i] = area(kendrick_sample[1, ],
                           kendrick_sample[2, ],
                           kendrick_sample[3, ])
  drake_areas[i] = area(drake_sample[1, ],
                        drake_sample[2, ],
                        drake_sample[3, ])
}

kendrick_areas = unlist(kendrick_areas)
drake_areas = unlist(drake_areas)

areas = data.frame(area = c(kendrick_areas, drake_areas),
                    artist = c(rep("kendrick", 1000),
                              rep("drake", 1000)))

ggplot(areas, aes(x=area, color=artist)) +
  stat_ecdf(linewidth=1) +
  labs(x="Triangle Area Measure", y="ECDF")

ks.test(drake_areas, kendrick_areas, alternative="greater")

