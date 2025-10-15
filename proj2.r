##### Generate the households
households_maker <- function(n,hmax=5){
  h <- sample(rep(1:n, sample(1:hmax, n, replace = TRUE))[1:n])
  return(h)
}
h<-households_maker(1000)
##### get.net()
get.net <- function(beta,h,nc = 15){
  n <- length(beta)
  beta_bar <- mean(beta)
  alink <- vector("list",n)
  for(i in 1:(n-1)){
    j <- (i+1):n
    valid <- (h[j] != h[i])
    if(!any(valid)) next
    j_valid <- j[valid]
    prij <- nc*beta[i]*beta[j_valid]/(beta_bar^2*(n-1))
    draws <- runif(length(praij))<prij
    linked <- j_valid[draws]
    if(length(linked)>0){
      alink[[i]] <- c(alink[[i]],linked)
      for (m in linked){
        alink[[m]] <- c(alink[[m]],i)  
      }
    }
  }
  return(alink)
}
######
nseir <- function(beta, h, alink,
                  alpha = c(.1, .01, .01),
                  delta = .2, gamma = .4,
                  nc = 15, nt = 100, pinf = .005) {
  
  beta_bar <- mean(beta)
  n <- length(beta)
  x <- rep(0, n)  ## initialize to S
  
  ## 初始感染者
  ni <- ceiling(pinf * n)
  inf_init <- sample(1:n, ni)
  x[inf_init] <- 2
  
  ## 统计量
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- n - ni
  I[1] <- ni
  
  ## 主时间循环
  for (i in 2:nt) {
    u <- runif(n)
    x[x == 2 & u < delta] <- 3   ## I -> R
    x[x == 1 & u < gamma] <- 2   ## E -> I
    
    ## 找出当前感染者
    I_now <- which(x == 2)
    ## 初始化未感染概率（每轮都重置）
    qh_notinfect <- qc_notinfect <- qr_notinfect <- rep(1, n)
    
    ## 对所有感染者累积传播作用
    for (k in I_now) {
      ## 家庭传播
      hmates <- which(h == h[k] & x == 0)
      if (length(hmates) > 0) {
        qh_notinfect[hmates] <- qh_notinfect[hmates] * (1 - alpha[1])
      }
      
      ## 网络传播
      contacts <- alink[[k]]
      if (length(contacts) > 0) {
        sus_contact <- contacts[x[contacts] == 0]
        if (length(sus_contact) > 0) {
          qc_notinfect[sus_contact] <- qc_notinfect[sus_contact] * (1 - alpha[2])
        }
      }
      
      ## 随机群体传播
      sus_random <- which(x == 0)
      if (length(sus_random) > 0) {
        qr_notinfect[sus_random] <- qr_notinfect[sus_random] *
          (1 - (alpha[3] * nc * beta[k] * beta[sus_random]) /
             (beta_bar^2 * (n - 1)))
      }
    }
    
    ## 综合感染概率
    q_infect <- 1 - qh_notinfect * qc_notinfect * qr_notinfect
    
    ## 对易感者进行随机感染判断
    sus <- which(x == 0)
    if (length(sus) > 0) {
      infected_new <- sus[runif(length(sus)) < q_infect[sus]]
      x[infected_new] <- 1
    }
    
    ## 记录群体状态
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  
  return(list(S = S, E = E, I = I, R = R, t = 1:nt))
}
# n <- 1000
# h <- households_maker(n)
# beta <- runif(n, 0, 1)
# alink <- get.net(beta, h, nc=15)

# out <- nseir(beta, h, alink, nt=100)
# plot(out$t, out$I, type="l", col="red", lwd=2, ylab="Count", xlab="Day")
# lines(out$t, out$S, col="blue")
# lines(out$t, out$R, col="green")
# legend("topright", legend=c("S","I","R"), col=c("blue","red","green"), lwd=2)
##### 4. Plotting function
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

##### 5. Compare 4 scenarios
set.seed(123)
n <- 1000
h <- households_maker(n)
beta <- runif(n, 0, 1)
alink <- get.net(beta, h, nc = 15)
beta_mean <- rep(mean(beta), n)

# Scenario 1: Full model
epi1 <- nseir(beta, h, alink, alpha = c(0.1, 0.01, 0.01), nt = 100)
# Scenario 2: Random mixing only
epi2 <- nseir(beta, h, alink, alpha = c(0, 0, 0.04), nt = 100)
# Scenario 3: Constant beta (structured)
epi3 <- nseir(beta_mean, h, alink, alpha = c(0.1, 0.01, 0.01), nt = 100)
# Scenario 4: Constant beta + random mixing
epi4 <- nseir(beta_mean, h, alink, alpha = c(0, 0, 0.04), nt = 100)

##### 6. Plot all 4 side by side #####
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_nseir(epi1, "1.Full model (household + network)")
plot_nseir(epi2, "2.Random mixing only")
plot_nseir(epi3, "3.Constant beta (structured)")
plot_nseir(epi4, "4.Constant beta + random mixing")
