# Load Data, convert peerID, get escapees, etc.

library(igraph)
library(lubridate)
library(Rfast)

data <- read.csv("~/Desktop/OO-paper/histories_BYU.csv")
leger <- read.csv("~/Desktop/OO-paper/participants_BYU.csv")

leger[,"p2p_id"] <- paste(leger[,"p2p_id"])
data[,"peer_id"] <- paste(data[,"peer_id"])

contacts <- subset(data, type == "contact")
end_time <- contacts$time + contacts$contact_length/1000
contacts <- contacts[, c("user_id", "peer_id", "time")]
contacts <- cbind(contacts, end_time)

convert_peer_ID <- function(x, leger){
  return(leger[which(leger[,"p2p_id"] == x), "id"][1])
}

contacts[,"peer_id"] <- sapply(contacts[,"peer_id"], convert_peer_ID, leger = leger)

contacts$time <- as.POSIXct(contacts$time, origin = "1970-01-01", tz = "America/Denver")
contacts$end_time <- as.POSIXct(contacts$end_time, origin = "1970-01-01", tz = "America/Denver")

lower <- pmin(contacts$user_id, contacts$peer_id)
upper <- pmax(contacts$user_id, contacts$peer_id)

contacts$user_id <- lower
contacts$peer_id <- upper

contacts <- contacts[with(contacts, order(user_id, peer_id, time)), ]

gap <- 30 #seconds

i=2
while (i < nrow(contacts)) {
  if(
    contacts[i,"user_id"] == contacts[i-1,"user_id"] & 
    contacts[i,"peer_id"] == contacts[i-1,"peer_id"] & 
    contacts[i,"time"] < (contacts[i-1,"end_time"] + gap)
  ){
    contacts[i-1, "end_time"] <- max(contacts[i-1, "end_time"], contacts[i, "end_time"])
    contacts <- contacts[-i, ]
  }else{
    i <- i+1
  }
}

contacts <- contacts[with(contacts, order(time)), ]
row.names(contacts) <- 1:nrow(contacts)
duration <- as.numeric(contacts$end_time - contacts$time)
contacts <- cbind(contacts, duration)

min_duration <- 60 #seconds
contacts <- subset(contacts, duration >= min_duration)

#Users
users <- leger$id

# Escapees
escape <- subset(data, out == "ESCAPED")
escape <- escape[, c("user_id", "time")]
escape$time <- as.POSIXct(escape$time, origin = "1970-01-01", tz = "America/Denver")

# Time info for escape

start = as.POSIXct(min(contacts$time), origin = "1970-01-01", tz = "America/Denver")
end = as.POSIXct(max(contacts$time), origin = "1970-01-01", tz = "America/Denver")

user_end <- function(i){
  if(users[i] %in% escape$user_id){
    return(max(escape$time[which(escape$user_id == users[i])]))
  }else{
    return(end)
  }
}

### STEP 1: model # contacts after D days on # of contacts after Y days, where the latter is known. 
#(This may not even be possible, we will try)




### try to solve model: total contacts(t) = A(1 - exp(-S/A t))

crv <- function(t, params){
  params[2]*(1 - exp(-params[1]*t/params[2]))
}

obj <- function(params){
  sum((sapply(xs, crv, params = params) - ys)^2)
}

xs <- c()
ys <- c()

for (i in 1:length(users)) {
  
  print(i)

  user = users[i]
  dates <- seq.POSIXt(start, user_end(i), "hour")
  
  y = c()
  for (d in 1:length(dates)) {
    subsample <- subset(contacts, time < dates[d] & (user_id == user | peer_id == user))
    subsample <- unique(subsample[,1:2])
    y[d] <- nrow(subsample)
  }
  
  xs <- c(xs, 0:(length(dates) - 1))
  ys <- c(ys, y)
  
}

nlm(obj, p = c(1,1))

# better idea
# simply fit curve to aggregate data, i.e. curve of best fit for every user at once. This will regress toward the mean by OLS
## this will make for a nice boxplot :)

# what about variance? could simply fit variance as function of modeled mean (or real mean)
means <- c()
vars <- c()

for (i in 1:length(seq.POSIXt(start, end, "hour"))) {
  means[i] <- mean(ys[which(xs == (i-1))])
  vars[i] <- var(ys[which(xs == (i-1))])
}

summary(lm(vars~means + 0))

ts <- 0:(length(means) - 1)

summary(lm(means~ts + 0))

# Step 2: test proportionate mixing assumption a few different ways
#way 1: node degree product as predictor of quantile of interaction length

# no escapees (for now)
full_users <- setdiff(users, escape$user_id)

real_adj <- matrix(0, nrow = length(full_users), ncol = length(full_users))

for (i in 1:length(full_users)) {
  for (j in i:length(full_users)) {
    sub <- subset(
      contacts,
      (user_id == full_users[i] & peer_id == full_users[j]) | (user_id == full_users[j] & peer_id == full_users[i])
    )
    real_adj[i,j] <- sum(sub$duration)
    real_adj[j,i] <- sum(sub$duration)
  }
}

bin_adj <- real_adj > 1

sort(upper_tri(real_adj))

deg <- colsums(bin_adj)

pred_adj <- outer(deg,deg)

preds <- sqrt(as.numeric(pred_adj))
deps <- as.numeric(bin_adj)

#mod <- glm(as.numeric(bin_adj)~preds, family = binomial(link='logit'))
#summary(mod)

#sum((as.numeric(bin_adj) - as.numeric(predict(lmod)))^2)

lmod <- lm(deps~preds+0)
summary(lmod)


## what doesnt work

preds2 <- preds[deps>0]
deps2 <- as.numeric(real_adj)[deps>0]

lmod2 <- lm(deps2~preds2+0)
summary(lmod2)

cor(deps2,preds2)

# this correlation is low :) so we proceed by 
## bootstrapping degrees, 
## bootstrapping edges based on lmod, 
## then bootstrapping interaction patterns based on 1 week data (repeated)



g <- graph_from_adjacency_matrix(bin_adj, mode = "undirected")
plot(g, vertex.size = 2, vertex.label = "")

### testing new package

library(distr)

dumb <- function(x){
  if(x<0){
    0
  }else if(x<1){
    x/2
  }else if(x<2){
    1/2
  }else if(x<3){
    x/2 - 1/2
  }else{
    1
  }
}

library(GoFKernel)

f1 <- inverse(dumb)

dumber <- function(x){
  sapply(x,dumb)
}

dist <-AbscontDistribution(d=dumber, withStand = TRUE)  # signature for a dist with pdf ~ p
rdist <- r(dist)
hist(rdist(1000))




### prior on vaccinated proportion

betapar <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha, beta))
}

bdens <- function(x){dbeta(x, 5.5, 1.833)}

plot(seq(0,1,0.01), dbeta(seq(0,1,0.01), 13.3125, 4.4375))

rbeta(1000, 13.3125, 4.4375)


##### OD
# target variance
v = 8

mom2 <- (p_contact_dist^2*v + (-1 + p_contact_dist)*p_contact_dist*(1 + r_contact_dist)*psi + (-1 + p_contact_dist)*(1 + r_contact_dist)*(-2 + p_contact_dist + (-1 + p_contact_dist)*r_contact_dist)*psi^2)/((-1 + p_contact_dist)*(1 + r_contact_dist)*(-3 - r_contact_dist + p_contact_dist*(2 + r_contact_dist)))

alpha <- ((-mom2)*psi + psi^2)/(mom2 - psi^2)
beta <- ((mom2 - psi)*(-1 + psi))/(mom2 - psi^2)


tries <- rlnorm(100000, -.5, sqrt(.5*2))
mean(tries)
var(tries)








