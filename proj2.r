# Chao You s2785107; Chenyi Jia s2792046; Kangxin Wei s2817655
# Our team members completed this assignment together, completing each step 
# through discussion.In the early stage, we read the notes together to 
# understand the requirements and learn functions. Later, we program together, 
# find bugs, and solve them together.

# https://github.com/sibei-git/Assignment2_group10.git
#-------------------------------------------------------------------------------

# SEIR Network Simulation Model
# This script simulates infectious disease spread in a structured population
# using an SEIR (Susceptible–Exposed–Infectious–Recovered) model. 

# Individuals can become infected through three transmission routes:
#   (1) household contacts, (2) network connections, and (3) random mixing.

# The model compares four scenarios to show how structure and heterogeneity
# affect epidemic dynamics.
# At the end, the code compares four different epidemic scenarios.
#-------------------------------------------------------------------------------

##### 1. households_maker()
# Creates a vector 'h' that assigns each individual to a household.
# Each household has a random size (between 1 and hmax).
households_maker <- function(n, hmax = 5) {
  h <- sample(rep(1:n, sample(1:hmax, n, replace = TRUE))[1:n])
  return(h)
}

##### 2. get.net() 
# Constructs a contact network ('alink') among individuals based on:
# Their transmission rate (beta)，Their household (h)
# Desired average number of connections (nc)

# The resulting 'alink' is a list where each element contains the indices
# of individuals connected to that person (undirected links).

# Each individual's household links are pre-defined to capture intra-family
# transmission. When generating additional contact edges, the algorithm ensures
# that no household members are included as network contacts. 
# This separation enforces mutual exclusivity between household and network
# transmission routes.
get.net <- function(beta, h, nc = 15) {
  n <- length(beta)
  beta_bar <- mean(beta)
  alink <- vector("list", n)
  # Loop summary:
  # This loop goes through each individual (i) and probabilistically creates
  # links to others (j) who are NOT in the same household. The probability of
  # forming a connection depends on both individuals' betas.
  for (i in 1:(n - 1)) {
    j <- (i + 1):n                    # Potential contacts for individual i
    valid <- (h[j] != h[i])           # Exclude same-household members
    if (!any(valid)) next             # Skip if no valid contacts exist
    
    j_valid <- j[valid]               # Valid contact indices
    prij <- nc * beta[i] * beta[j_valid] / (beta_bar^2 * (n - 1))  # Link prob.
    draws <- runif(length(prij)) < prij
    linked <- j_valid[draws]          # Realized links
    
    # Store symmetric links (both directions)
    if (length(linked) > 0) {
      alink[[i]] <- c(alink[[i]], linked)
      for (m in linked) {
        alink[[m]] <- c(alink[[m]], i)
      }
    }
  }
  
  return(alink)
}

##### 3. nseir() 
# Simulates SEIR disease dynamics in a structured population.
# Transmission occurs via three mechanisms:
# (1) Household (alpha[1]), (2) Network contacts (alpha[2])
# (3) Random mixing (alpha[3])

# Returns a list of time series (S, E, I, R).
nseir <- function(beta, h, alink,
                  alpha = c(.1, .01, .01),
                  delta = .2, gamma = .4,
                  nc = 15, nt = 100, pinf = .005) {
  
  beta_bar <- mean(beta)
  n <- length(beta)
  x <- rep(0, n)  # 0=S, 1=E, 2=I, 3=R
  
  # Initialize a small proportion of individuals as infectious
  ni <- ceiling(pinf * n)
  inf_init <- sample(1:n, ni)
  x[inf_init] <- 2
  
  # Storage for population in each state
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- n - ni
  I[1] <- ni
  
  # Main simulation loop over time
  for (i in 2:nt) {
    u <- runif(n)
    
    # Transition updates
    x[x == 2 & u < delta] <- 3   # I -> R with prob delta
    x[x == 1 & u < gamma] <- 2   # E -> I with prob gamma
    
    I_now <- which(x == 2)       # Current infectious individuals
    
    # Reset "not infected" probabilities
    qh_notinfect <- qc_notinfect <- qr_notinfect <- rep(1, n)
    
    # Loop summary:
    # For each infectious individual:
    #   1. Infects susceptible household members (alpha[1])
    #   2. Infects network contacts (alpha[2])
    #   3. Infects random population members (alpha[3])
    
    # The probabilities accumulate multiplicatively for each susceptible person.
    
    # At each step, calculate the exposure risk for the 3 scenarios separately:
    # (1) household contact, (2) network contact, (3) random mixing.
    # Household or network infection routes are mutually exclusive.
    # Because in the function get.net, it is ensured that the infected person's 
    # contact network does not include their family members.
    for (k in I_now) {
      # Household transmission
      hmates <- which(h == h[k] & x == 0)
      if (length(hmates) > 0) {
        qh_notinfect[hmates] <- qh_notinfect[hmates] * (1 - alpha[1])
      }
      
      # Network transmission
      contacts <- alink[[k]]
      if (length(contacts) > 0) {
        sus_contact <- contacts[x[contacts] == 0]
        if (length(sus_contact) > 0) {
          qc_notinfect[sus_contact] <- qc_notinfect[sus_contact] * (1 - alpha[2])
        }
      }
      
      # Random mixing transmission
      sus_random <- which(x == 0)
      if (length(sus_random) > 0) {
        qr_notinfect[sus_random] <- qr_notinfect[sus_random] *
          (1 - (alpha[3] * nc * beta[k] * beta[sus_random]) /
             (beta_bar^2 * (n - 1)))
      }
    }
    
    # Calculate total infection probability
    q_infect <- 1 - qh_notinfect * qc_notinfect * qr_notinfect
    
    # Infect new susceptibles
    sus <- which(x == 0)
    if (length(sus) > 0) {
      infected_new <- sus[runif(length(sus)) < q_infect[sus]]
      x[infected_new] <- 1
    }
    
    # Record population state counts
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  
  return(list(S = S, E = E, I = I, R = R, t = 1:nt))
}

##### 4. plot_nseir() 
# Generates a time series plot of the SEIR simulation results.
# Displays the counts of S, E, I, R over time.
plot_nseir <- function(epi, title = "SEIR Simulation Dynamics") {
  t <- epi$t
  nmax <- max(c(epi$S, epi$E, epi$I, epi$R))
  cols <- c(S = "green", E = "blue", I = "red", R = "orange")
  
  plot(t, epi$S, type = "n", ylim = c(0, nmax),
       xlab = "Time (days)", ylab = "Population count", main = title)
  grid(col = "gray85", lty = "dotted")
  
  lines(t, epi$S, col = cols["S"], lwd = 2)
  lines(t, epi$E, col = cols["E"], lwd = 2)
  lines(t, epi$I, col = cols["I"], lwd = 2)
  lines(t, epi$R, col = cols["R"], lwd = 2)
  
  legend("right", legend = c("S", "E", "I", "R"), col = cols, lwd = 2, bty = "n")
}

##### 5. Scenario Comparison 
# Compare epidemic dynamics under four different scenarios
n <- 10000
h <- households_maker(n)
beta <- runif(n, 0, 1)
alink <- get.net(beta, h, nc = 15)
beta_mean <- rep(mean(beta), n)

# Scenario 1: Full model (household + network + random mixing)
epi1 <- nseir(beta, h, alink, alpha = c(0.1, 0.01, 0.01), nt = 100)

# Scenario 2: Random mixing only
epi2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04), nt = 100)

# Scenario 3: Constant beta (structured contacts)
epi3 <- nseir(beta_mean, h, alink, alpha = c(0.1, 0.01, 0.01), nt = 100)

# Scenario 4: Constant beta + random mixing
epi4 <- nseir(beta_mean, h, alink, alpha = c(0, 0, 0.04), nt = 100)


##### 6. Plot Results 
# Display all four scenarios side-by-side
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_nseir(epi1, "1. Full model (household + network)")
plot_nseir(epi2, "2. Random mixing only")
plot_nseir(epi3, "3. Constant beta (structured)")
plot_nseir(epi4, "4. Constant beta + random mixing")

#Analysis:

# Figure 1: Full model (household + network)
# With both household and network structure, infection spreads more slowly 
# and peaks later with a lower amplitude.
# The prolonged infectious period reflects localize transmission within 
# connected groups.
# Structural clustering limits transmission and flattens the epidemic curve.

# Figure 2: Random mixing only
# Only assuming random contacts accelerates the outbreak dramatically.
# The number of susceptible individuals drops rapidly, 
# and the infection peak occurs earlier and higher.
# Random mixing eliminates local constraints, 
# leading to faster and broader epidemics.

# Figure 3: Constant β (structured)
# Keeping the structure but fixing β makes the curves smoother 
# and slightly lowers the peak.
# Without variability in infectiousness, super-spreader effects disappear 
# and transmission becomes more uniform.
# Removing individual effects reduces stochastic fluctuations 
# and produces more predictable dynamics.

# Figure 4: Constant β + Random mixing
# With no structure and constant β, the system behaves as a basic SEIR model.
# The curves are symmetric and smooth, with the highest and earliest 
# infection peak.
# Represents the theoretical upper limit of epidemic speed in a well-mixed
# population.

# Comparative Conclusion
# Social and household structures slower transmission and lower epidemic peaks.
# Individual variability in β increases epidemic intensity and accelerates
# early spread.
# Social structure governs spread speed; 
# Individual variability in the transmission rate parameter β shapes peak height.
# Real epidemics lie between fully structured and fully mixed extremes.

# Conclusion
# Structural contact patterns mitigate transmission speed, 
# while heterogeneity amplifies it.
# Incorporating both mechanisms yields a more realistic representation of 
# epidemic dynamics.