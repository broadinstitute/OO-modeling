library(igraph)
library(lubridate)

## fucking around

g <- erdos.renyi.game(358, 1/358, "gnp")

el <- as_edgelist(g)

secondaries <- ego(g, order = 2, mindist = 2)
secondaries <- lapply(secondaries, as.numeric)


get_edges <- function(v){
  
}

adj <- as.matrix(as_adjacency_matrix(g))

for (i in 1:358) {
  new_edges <- rbinom(length(secondaries[[i]]), 1, 0.9)
  adj[i, secondaries[[i]]] <- new_edges
  adj[secondaries[[i]], i] <- new_edges
}

g2 <- graph_from_adjacency_matrix(adj, mode = "undirected")

hist(degree(g2))
mean(degree(g2))
var(degree(g2))
transitivity(g2)


###### CMU

data <- read.csv("~/Desktop/OO-paper/histories_CMU.csv")
leger <- read.csv("~/Desktop/OO-paper/participants_CMU.csv")

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

gap <- 2 #seconds

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


####### get weighted average contact length
escape <- subset(data, out == "ESCAPED")
escape <- escape[, c("user_id", "time")]
escape$time <- as.POSIXct(escape$time, origin = "1970-01-01", tz = "America/Denver")

start = as.POSIXct(min(data$time), origin = "1970-01-01", tz = "America/Denver")
end = as.POSIXct(max(data$time), origin = "1970-01-01", tz = "America/Denver")

contacts_per_day <- function(v){
  
  if(v %in% escape$user_id){
    escape_time <- min(escape[which(escape$user_id == v), "time"])
  }else{
    escape_time <- end
  }
  
  nrow(subset(contacts, (user_id == v | peer_id == v) & time <= escape_time)) / 
    as.numeric(difftime(escape_time, start, units = "days"))
}

get_weight <- function(v){
  if(v %in% escape$user_id){
    escape_time <- min(escape[which(escape$user_id == v), "time"])
  }else{
    escape_time <- end
  }
  
  as.numeric(difftime(escape_time, start, units = "days"))
}


daily_contacts <- sapply(leger$id, contacts_per_day)
weight <- sapply(leger$id, get_weight)
weight <- weight / max(weight)
hist(daily_contacts, breaks = seq(0,300,1))

### test proportionate mixing hypothesis - VERSION 1, number of contacts only

people <- leger$id
people <- people[weight == 1]

unique_contacts <- contacts[, 1:2]
unique_contacts <- unique(unique_contacts)
unique_contacts <- subset(unique_contacts, user_id %in% people & peer_id %in% people)
row.names(unique_contacts) <- 1:nrow(unique_contacts)

promiscuity <- c()
for (i in 1:length(people)) {
  promiscuity[i] <- nrow(subset(unique_contacts, user_id == people[i] | peer_id== people[i]))
}
promiscuity <- promiscuity/max(promiscuity)

predicted_adj <- outer(promiscuity, promiscuity) * 1.5
diag(predicted_adj) <- 0
new <- runif(271*135) < predicted_adj[upper.tri(predicted_adj)]
predicted_adj[upper.tri(predicted_adj)] <- new
predicted_adj <- t(predicted_adj)
predicted_adj[upper.tri(predicted_adj)] <- new
sum(predicted_adj)

gp <- graph_from_adjacency_matrix(predicted_adj, mode = "undirected")
largest_cliques(gp)
transitivity(gp)

plot(gp, vertex.label = "", vertex.size = 3)



real_adj <- matrix(0, nrow = length(people), ncol = length(people))
for (r in 1:nrow(unique_contacts)) {
  real_adj[which(people == unique_contacts[r,1]), which(people == unique_contacts[r,2])] <- 1
  real_adj[which(people == unique_contacts[r,2]), which(people == unique_contacts[r,1])] <- 1
}

gg <- graph_from_adjacency_matrix(real_adj, mode = "undirected")
plot(gg, vertex.label = "", vertex.size = 3)

rg <- erdos.renyi.game(271, 1125/(271*135), "gnp")
E(rg)

transitivity(rg)
transitivity(gg)

clique_data <- c()
for (i in 1:11) {
  clique_data[i] <- length(cliques(gg, min = i, max = i))
}
clique_data

mc <- max_cliques(gg)
nmc <- c()
for (i in 1:length(mc)) {
  nmc[i] <- length(as.numeric(mc[[i]]))
}


plot(as.numeric(predicted_adj), as.numeric(real_adj))
summary(lm(as.numeric(real_adj) ~ as.numeric(predicted_adj)))

### VERSION 2 - total interaction time

# new idea
# if you have one heavy edge, how good is that at predicting dist of other edges?
# could help contradict iid assumption
# could inform iterative edge weight development

complete_contacts <- subset(contacts, user_id %in% people & peer_id %in% people)

promiscuity <- c()

for (i in 1:length(people)) {
  promiscuity[i] <- sum(subset(complete_contacts, user_id == people[i] | peer_id == people[i])$duration)
}

predicted_adj <- outer(promiscuity, promiscuity)

real_adj <- matrix(0, nrow = length(people), ncol = length(people))
for (r in 1:nrow(unique_contacts)) {
  sub <- subset(complete_contacts, (user_id == unique_contacts[r,1] & peer_id == unique_contacts[r,2]) | (user_id == unique_contacts[r,2] & peer_id == unique_contacts[r,1]))
  real_adj[which(people == unique_contacts[r,1]), which(people == unique_contacts[r,2])] <- sum(sub$duration)
  real_adj[which(people == unique_contacts[r,2]), which(people == unique_contacts[r,1])] <- sum(sub$duration)
}

# a/l = 2770

library(ggplot2)

crv <- function(x){dlnorm(x, 6.1, 1.74)}

dat <- data.frame(edw = edw)
ggp <- ggplot(data = dat, aes(x = edw)) + geom_histogram(aes(y = ..density..), binwidth = 10) + xlim(0,2000) + stat_function(fun = crv)
ggp

plot(as.numeric(predicted_adj), as.numeric(real_adj))
summary(lm(as.numeric(real_adj) ~ as.numeric(predicted_adj)))

### Version 3 - iid? testing


df <- data.frame()
for (i in 1:nrow(real_adj)) {
  
    n_contacts <- sum(real_adj[i,] > 0)
    lens <- real_adj[i,][which(real_adj[i,] > 0)]
    
    newdata <- data.frame(n = rep(n_contacts, length(lens)), length = lens)
    
    
    df <- rbind(df, newdata)
    
}
colnames(df) <- c("n", "length")

summary(lm(length ~ n, df))
plot(df$n, df$length)

# Hmmmm, this isn't very promising. Let's analyze 2 edges as predictor of 3rd edge

prd <- data.frame()
for (i in 1:nrow(real_adj)) {
  for (j in i:nrow(real_adj)) {
    mutual <- intersect(which(real_adj[i,] > 0), which(real_adj[j,] > 0))
    
    pp <- sqrt(real_adj[i, mutual] * real_adj[j, mutual])
    
    newdata <- data.frame(predictor = pp, actual = rep(real_adj[i,j], length(pp)))
    prd <- rbind(prd, newdata)
  }
}

plot(prd$actual, prd$predictor)



##############

weighted_count <- function(n){
  indices <- which(n <= daily_contacts & daily_contacts < n + 1)
  sum(weight[indices])
}

plot(0:300, sapply(0:300, weighted_count))

weighted.mean(daily_contacts, weight)






# 
unique_contacts <- unique(contacts[,1:2])

graph <- graph_from_data_frame(unique_contacts, directed = FALSE)
graph <- add_vertices(graph, 358 - length(degree(graph)))
subgraph <- induced_subgraph(graph, which(rbinom(length(degree(graph)), 1, 0.5)==1))

var(degree(graph))
var(degree(subgraph))

sum(degree(graph) == 0)/length(degree(graph))
sum(degree(subgraph) == 0)/length(degree(subgraph))

hist(degree(graph), breaks = seq(0,100,1))
hist(degree(subgraph))


hist(
  c(
    as.numeric(degree(graph)), 
    rep(0, 358 - length(degree(graph)))
  ),
  breaks = seq(0,100,1)
)



############








# write.csv(contacts, "~/Desktop/CMU_processed.csv")

hist(as.numeric(duration)[as.numeric(duration) < 1000], breaks = seq(0,1000,10))











