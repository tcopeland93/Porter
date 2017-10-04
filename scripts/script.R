'                                      Tim Copeland                                      '
'                            Empirical Industrial Organization                           '
'                                      Problem Set 2                                     '

# To run code: source("./scripts/script.R")

suppressMessages(library(systemfit))
suppressMessages(library(stats4))

# Question 2: Produce Summary Statistics

railway <- read.csv('./data/railway.csv',header = TRUE)
gr <- railway$gr
tqg <- railway$tqg			# Abbreviate variable names for convenience
lakes <- railway$lakes
po <- railway$po
dm1 <- railway$dm1
dm2 <- railway$dm2
dm3 <- railway$dm3
dm4 <- railway$dm4
s1 <- railway$seas1
s2 <- railway$seas2
s3 <- railway$seas3
s4 <- railway$seas4
s5 <- railway$seas5
s6 <- railway$seas6
s7 <- railway$seas7
s8 <- railway$seas8
s9 <- railway$seas9
s10 <- railway$seas10
s11 <- railway$seas11
s12 <- railway$seas12



specify_decimal <- function(x,k) trimws(format(round(x,k), nsmall = k))

stats <-  function(a){
	c(
	specify_decimal(mean(a), 4), 
	specify_decimal(sd(a), 4), 
	min(a), 
	max(a)
	)
	}

# Summary statistics table

sumstats <- matrix(c(stats(gr), stats(tqg), stats(lakes), stats(po)), ncol = 4, byrow = TRUE)
colnames(sumstats) <- c("mean", "standard deviation", "min value", "max value")
rownames(sumstats) <- c("gr", "tqg", "lakes", "po")



# Question 3: Produce IV Estimates


demand.equation <- log(tqg)~log(gr)+lakes+s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12
supply.equation <- log(gr)~log(tqg)+dm1+dm2+dm3+dm4+po+s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12
inst <- ~lakes+s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+dm1+dm2+dm3+dm4+po

twosls <- systemfit(list(demand.equation,supply.equation), inst = inst, method = '2SLS')

IV.table.estimates <- matrix (specify_decimal(
	c(
	twosls$coefficients[1], twosls$coefficients[3], twosls$coefficients[2], NA, NA, NA, 
	NA, NA, NA, twosls$coefficients[16], NA, NA, twosls$coefficients[18], 
	twosls$coefficients[19], twosls$coefficients[20], twosls$coefficients[21],
	twosls$coefficients[22], twosls$coefficients[17]), 
	3),
	ncol = 2)

rownames(IV.table.estimates) <- c("Intercept", "Lakes", "GR", "DM1", "DM2", "DM3", "DM4", "PO", "TQG")
colnames(IV.table.estimates) <- c("Demand", "Supply")




# Question 4: Maximum Likelihood Estimation

y <- matrix(c(log(tqg),log(gr)),ncol=2)
X <- matrix(c(rep(1,328),lakes,dm1,dm2,dm3,dm4),ncol=6)
I <- matrix(po,ncol=1)



LogLL.givenI <- function(a1,a2,a3,b1,b2,b3,b4,b5,b6,delta){
	B <- matrix(c(1,-a2,-b2,1),ncol=2)
	gma <- matrix(c(a1,a3,0,0,0,0,b1,0,b3,b4,b5,b6),ncol=2)
	u <- y %*% B - X %*% gma - I %*% matrix(c(0,delta),ncol=2)
	Sigma <<- t(u) %*% u
	-sum(log((1/(2*pi))*(det(Sigma)**-.5)*norm(B, type = "f")*exp(-0.5*diag(u %*% solve(Sigma) %*% t(u)))))
}

starting <-  list(a1 = 9.177, a2 = -.735, a3 = -.437, b1 = -3.975, b2 = 0.253, b3 = -0.202, 
b4 = -.173, b5 = -.319, b6 = -.208, delta = 0.368)




#LLest <- mle(LogLL.givenI, start = starting, method = "Nelder-Mead")
#estimates <- coef(LLest)

#Code to test LL function: LogLL.estllambda(9.1,-.735,-.437,-3.975,.253,-.202,-.173,-.319,-.208,.368)

llambda = .5           #sum(I)/length(I)
I <- rep(c(.25,.75), 164)


LogLL.estllambda <- function(a1,a2,a3,b1,b2,b3,b4,b5,b6,delta){
	B <- matrix(c(1,-a2,-b2,1),ncol=2)
	gma <- matrix(c(a1,a3,0,0,0,0,b1,0,b3,b4,b5,b6),ncol=2)
	u <- y %*% B - X %*% gma - I %*% matrix(c(0,delta),ncol=2)
	Sigma <<- t(u) %*% u
	-sum(log((1/(2*pi))*(det(Sigma)**-.5)*norm(B, type = 'f')*
	(llambda*exp(-0.5*diag(u %*% solve(Sigma) %*% t(u)))+(1-llambda)*exp(-0.5*diag(u %*% solve(Sigma) %*% t(u))))))
}

LLest2 <- mle(LogLL.estllambda, start=starting, method = 'Nelder-Mead')
print(LLest2)
estimates <- coef(LLest2)

B <- matrix(c(1,-estimates[2],-estimates[5],1),ncol=2)
gma <- matrix(c(estimates[1],estimates[3],0,0,0,0,estimates[4],0,estimates[6],estimates[7],estimates[8],estimates[9]),ncol=2)
delta <- estimates[10]
omega <- list(a1 = as.numeric(estimates[1]), a2 = as.numeric(estimates[2]), a3 = as.numeric(estimates[3]), 
		b1 = as.numeric(estimates[4]), b2 = as.numeric(estimates[5]), b3 = as.numeric(estimates[6]), 
		b4 = as.numeric(estimates[7]), b5 = as.numeric(estimates[8]), b6 = as.numeric(estimates[9]), 
		delta = as.numeric(estimates[10]))

w <- function(){
	u1 <<- y %*% B - X %*% gma - matrix(rep(1,328), ncol=1) %*% matrix(c(0,delta),ncol=2)
	u0 <<- y %*% B - X %*% gma
	Sigma.1 <<- t(u1) %*% u1
	Sigma.0 <<- t(u0) %*% u0
	num <<- llambda*(1/(2*pi))*(det(Sigma.1)**-.5)*norm(B, type = "f")*(exp(-0.5*diag(u1 %*% solve(Sigma.1) %*% t(u1))))
	den <<- num + (1-llambda)*(1/(2*pi))*(det(Sigma.0)**-.5)*norm(B, type = "f")*(exp(-0.5*diag(u0 %*% solve(Sigma.0) %*% t(u0))))
	num/den
}

I <- w()



print(w())

estimates.log <- matrix(c(rep(c(0,1),5),estimates), ncol=10, byrow = TRUE)
I.log <- matrix(c(rep(c(0,1), 164),I),ncol = 328, byrow = TRUE)

while(cor(I.log[dim(I.log)[1],],I.log[(dim(I.log)[1]-1),])<0.999){
	LLest <- mle(LogLL.estllambda, start = omega, method = "Nelder-Mead")
	estimates <- coef(LLest)
	B <- matrix(c(1,-estimates[2],-estimates[5],1),ncol=2)
	gma <- matrix(c(estimates[1],estimates[3],0,0,0,0,estimates[4],0,estimates[6],estimates[7],estimates[8],estimates[9]),ncol=2)
	delta <- estimates[10]
	omega <- list(a1 = as.numeric(estimates[1]), a2 = as.numeric(estimates[2]), a3 = as.numeric(estimates[3]), 
		b1 = as.numeric(estimates[4]), b2 = as.numeric(estimates[5]), b3 = as.numeric(estimates[6]), 
		b4 = as.numeric(estimates[7]), b5 = as.numeric(estimates[8]), b6 = as.numeric(estimates[9]), 
		delta = as.numeric(estimates[10]))
	
	print(estimates)
	
	estimates.log <- rbind(estimates.log,estimates)
	
	
	
	w()
	I <- w()
	I.log <- rbind(I.log, I)
	print(w()[1:20])
	print(llambda)
	llambda = sum(I)/length(I)
}




"cat('\n
------------------------------------  Tim Copeland  ------------------------------------
\n
------------------------------------  Empirical IO  ------------------------------------
\n
------------------------------------  Problem Set 2  -----------------------------------
\n\n

Problem 1: Model:
-----------------

The authors use a 2SLS technique, using a reported presence of collusion, to estimate a 
supply and demand model with the given data. 

Next, they construct a Maximum Likelihood Estimation model, allowing the model to detect
or predict collusive behavior that fits the data.  


Problem 2: Summary Statistics:
------------------------------

')


print(sumstats)

cat('\n

Problem 3: IV Estimation:
------------------------- 

'

)

print(IV.table.estimates)"



