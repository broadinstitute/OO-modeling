# Load and clean data

library(igraph)
library(lubridate)
library(Rfast)
library(mixdist)
library(parallel)
library(ggplot2)

# Set directory to the location of "OO-modeling" folder
setwd("~/Desktop/OO-modeling")

# Read in OO simulation backend
data <- read.csv("./histories_CMU.csv")

# Read in key to convert between simulation IDs and P2P IDs
leger <- read.csv("./participants_CMU.csv")
leger[,"p2p_id"] <- paste(leger[,"p2p_id"])
data[,"peer_id"] <- paste(data[,"peer_id"])

# Limit data to contacts only
contacts <- subset(data, type == "contact")

# Reformat contact duration
contacts$contact_length <- round(milliseconds(contacts$contact_length))

# Compute end time per contact
end_time <- contacts$time + contacts$contact_length

# Subset dataset to relevant categories
contacts <- contacts[, c("user_id", "peer_id", "time")]
contacts <- cbind(contacts, end_time)

# Convert P2P ID to sim ID
convert_peer_ID <- function(x, leger){
  return(leger[which(leger[,"p2p_id"] == x), "id"][1])
}

# Convert old ID to new sim ID
convert_user_ID <- function(x, leger){
  return(
    leger$id[which(leger$p2p_id == leger$p2p_id[which(leger$id == x)])[1]]
  )
}

contacts[,"user_id"] <- sapply(contacts[,"user_id"], convert_user_ID, leger = leger)
contacts[,"peer_id"] <- sapply(contacts[,"peer_id"], convert_peer_ID, leger = leger)

# Reformat contact time
contacts$time <- as.POSIXct(contacts$time, origin = "1970-01-01", tz = "America/Denver")
contacts$end_time <- as.POSIXct(contacts$end_time, origin = "1970-01-01", tz = "America/Denver")

# Make user ID consistently less than peer ID
lower <- pmin(contacts$user_id, contacts$peer_id)
upper <- pmax(contacts$user_id, contacts$peer_id)
contacts$user_id <- lower
contacts$peer_id <- upper

contacts <- contacts[with(contacts, order(user_id, peer_id, time)), ]

# Length between contacts to merge into 1 contact
gap <- 30 #seconds

# Merge same contacts
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

# Clean up datasheet and merge in relevant metrics
contacts <- contacts[with(contacts, order(time)), ]
row.names(contacts) <- 1:nrow(contacts)
duration <- contacts$end_time - contacts$time
contacts <- cbind(contacts, duration)

###

# Who escaped?
escape <- subset(data, out == "ESCAPED")
escape <- escape[, c("user_id", "time")]
escape$time <- as.POSIXct(escape$time, origin = "1970-01-01", tz = "America/Denver")
escape <- escape[order(escape$time),]
escape[,"user_id"] <- sapply(escape[,"user_id"], convert_user_ID, leger = leger)
row.names(escape) <- 1:nrow(escape)

# When did the simulation start?
start = as.POSIXct(min(contacts$time), origin = "1970-01-01", tz = "America/Denver")

# When did the simulation end?
end <- as.POSIXct(max(contacts$time), origin = "1970-01-01", tz = "America/Denver")

# Who escaped BEFORE the end?
escapees <- escape$user_id[escape$time < end]

# All participants
leger[,"id"] <- sapply(leger[,"id"], convert_user_ID, leger = leger)
users <- setdiff(leger$id, escapees)

# Subset contacts
contacts <- subset(contacts, (user_id %in% users & peer_id %in% users))
contacts <- subset(contacts, time < end)

# Modify contacts to break by hour

i=1
while (i <= nrow(contacts)) {
  if(
    (floor_date(contacts[i,"time"], unit = "hours") != floor_date(contacts[i,"end_time"], unit = "hours")) & (contacts[i,"end_time"] != floor_date(contacts[i,"end_time"], unit = "hours"))
  ){
    oldrow <- contacts[i, ]
    oldrow$end_time <- ceiling_date(oldrow$time, unit = "hours", change_on_boundary = TRUE)
    oldrow$duration <- difftime(oldrow$end_time, oldrow$time, units = "secs")
    
    newrow <- contacts[i, ]
    newrow$time <- ceiling_date(newrow$time, unit = "hours", change_on_boundary = TRUE)
    newrow$duration <- difftime(newrow$end_time, newrow$time, units = "secs")
    
    contacts[i, ] <- oldrow
    contacts <- rbind(contacts, newrow)
    
  }
  
  i <- i+1
  
}

contacts <- contacts[with(contacts, order(time)), ]

# Merge contacts within an hour

breaks <- seq(floor_date(min(contacts$time), unit = "hours"), ceiling_date(max(contacts$end_time), unit = "hours"), "hours")

df <- data.frame()

for (i in 1:(length(breaks) - 1)) {
  subcontacts <- subset(contacts, breaks[i] <= time & end_time <= breaks[i+1])
  subcontacts <- subcontacts[with(subcontacts, order(user_id, peer_id, time)), ]
  j=2
  while (j <= nrow(subcontacts)) {
    if(
      subcontacts[j,"user_id"] == subcontacts[j-1,"user_id"] & 
      subcontacts[j,"peer_id"] == subcontacts[j-1,"peer_id"]
    ){
      subcontacts[j-1, "end_time"] <- max(subcontacts[j-1, "end_time"], subcontacts[j, "end_time"])
      subcontacts[j-1, "duration"] <- round(subcontacts[j-1, "end_time"] - subcontacts[j-1, "time"])
      subcontacts <- subcontacts[-j, ]
    }else{
      j <- j+1
    }
  }
  
  hour_start <- rep(breaks[i], nrow(subcontacts))
  subcontacts <- subcontacts[, c("user_id", "peer_id", "duration")]
  subcontacts <- cbind(subcontacts, hour_start)
  
  df <- rbind(df, subcontacts)
}

row.names(df) <- 1:nrow(df)

# Rename the people

#users <- sort(unique(c(df$user_id, df$peer_id)))

rename <- function(user){
  which(users == user)
}

df$user_id <- sapply(df$user_id, rename)
df$peer_id <- sapply(df$peer_id, rename)


df$duration <- as.numeric(df$duration)

# Contacts and durations. 1st index is timeblock, 2nd index is duration of contact with each

adj <- list()
dur <- list()

N <- length(users)

for (i in 1:(length(breaks) - 1)) {
  print(i)
  adj[[i]] <- list()
  dur[[i]] <- list()
  sub <- subset(df, hour_start == breaks[i])
  mtx <- matrix(0, N, N)
  mtx[as.matrix(sub[, 1:2])] <- sub$duration
  mtx[as.matrix(sub[, 2:1])] <- sub$duration
  for (j in 1:N) {
    adj[[i]][[j]] <- which(mtx[,j] > 0)
    dur[[i]][[j]] <- mtx[,j][which(mtx[,j] > 0)]
  }
}



# how many total contacts per person?
contacts_per_person <- c()
unique_contacts_per_person <- c()
dur_contacts_per_person <- c()

for (i in 1:N) {
  sdf <- subset(df, user_id==i | peer_id==i)
  contacts_per_person[i] <- nrow(sdf)
  unique_contacts_per_person[i] <- nrow(unique(sdf[,1:2]))
  dur_contacts_per_person[i] <- sum(sdf$duration)
}

##### now for the model

avg_daily_contacts <- mean(contacts_per_person / as.numeric(difftime(breaks[length(breaks)], breaks[1], units = "days")))

get_lambda <- function(lambda, R0, avg_infection){
  (R0 - mean(1 - exp(-lambda*df$duration)) * avg_daily_contacts * avg_infection)^2
}

pickup <- function(user, time, I, lambda){
  neighbors <- adj[[time-1]][[user]]
  durations <- dur[[time-1]][[user]]
  prob <- 1 - prod(1 - I[time-1, neighbors]*(1 - exp(-lambda*durations))) # given susceptible
}


###

delta <- 1/(5*24) #hours
gamma <- 1/(10*24) #hours
avg_infection <- 1/(gamma*24) # days duration

R0 <- 2.5


lambda <- suppressWarnings(
  nlm(get_lambda, R0 / (mean(df$duration) * avg_daily_contacts * avg_infection), R0 = R0, avg_infection = avg_infection)$estimate
)

#lambda <- lambda * (1 - (input$vax * input$p_vax)/10000) * (1 - (input$mask * input$p_mask)/10000)^2


S <- matrix(ncol = N, nrow = length(breaks))
E <- matrix(ncol = N, nrow = length(breaks))
I <- matrix(ncol = N, nrow = length(breaks))
R <- matrix(ncol = N, nrow = length(breaks))

I0 <- 0.02

S[1, ]  <- rep(1 - I0, N) # probability susceptible
E[1, ] <- rep(0, N) # probability exposed
I[1, ] <- rep(I0, N) # probability infectious
R[1, ] <- rep(0, N) # probability recovered

for (i in 2:length(breaks)) {
  
  subdf <- subset(df, hour_start == breaks[i-1])
  subusers <- unique(c(subdf$user_id, subdf$peer_id))
  
  p_transmit <- rep(0, N)
  if(length(subusers)>0){
    p_transmit[subusers] <- sapply(subusers, pickup, time = i, I = I, lambda = lambda)
  }
  
  
  S[i, ] <- S[i-1, ] - S[i-1, ]*p_transmit
  E[i, ] <- E[i-1, ] + S[i-1, ]*p_transmit - delta*E[i-1, ]
  I[i, ] <- I[i-1, ] + delta*E[i-1, ] - gamma*I[i-1, ]
  R[i, ] <- R[i-1, ] + gamma*I[i-1, ]
  
}

final <- 1 - S[nrow(S), ]

save(final, file = "final_CMU.RData")

