
X = matrix(
	c(2.771244718,1.784783929,0,
		1.728571309,1.169761413,0,
		3.678319846,2.81281357,0,
		3.961043357,2.61995032,0,
		2.999208922,2.209014212,0,
		7.497545867,3.162953546,1,
		9.00220326,3.339047188,1,
		7.444542326,0.476683375,1,
		10.12493903,3.234550982,1,
		6.642287351,3.319983761,1),
	nrow=10,
	ncol=3,
	byrow=TRUE)


sgn<-function(W,X,row){
	# W is the weight matrix
	# x is one row of the the innput matrix
	# sgn is the m-bit m-bit vector of potential split decisions (h)
	x = X[row,c(1:ncol(X)-1)]
	out_mat = W %*% x -1
	sign  = sign(out_mat)
	return(sign)}

f <- function(h){
	# function takes as arguments:
	# an m-bit vector of potential split decisions (h)
	# function returns: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	h[h==-1]<-0
	h<-list(h[1],h[2:3])
	
	# TODO: Rewrite this^ so it works for a longer m decision tree

	# find number of rows from input
	N <- length(h)

	# go through and pick out the values in each tree that are valid based
	# on previous route
	out <- c(h[[1]], rep(0, N-1))
	for(i in 2:N){
	out[i] <- h[[i]][sum(out[i:(i-1)] * 2^(i-1)/(2^((i-1):1))) + 1]
	}

	# now find the final position in the bottom row and return as a vector
	out_pos <- sum(out * 2^N/(2^(1:N))) + 1
	full_vec <- rep(0, 2^N)
	full_vec[out_pos] <- 1

	return(full_vec)
}

update_probs<-function(theta,ID_vec){
	# function takes as arguments: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	# an 1 by k length matrix of probabilities
	# function returns: 
	# an m+1 by k matrix of softmax processed probabilities of p(y= l|j)	
	theta_vec<-theta[grep(1,ID_vec),]
	denominator<-do.call(sum,lapply(theta_vec,exp))
	new_probs<-matrix(unlist(lapply(theta_vec, function(x) exp(x)/denominator)))
	new_probs = (ID_vec) %*% t(new_probs)
	return(new_probs)
}


loss<-function(X,samp_row,theta,g){
	# function takes as arguments: 
	# the dataset and row being evaluated
	# predicted probabilities from the update_probs softmax
	# function returns: 
	# log loss value
	true_base = as.numeric(X[samp_row,ncol(X)])
	probs<-t(theta)%*%f(g)
	prob<-probs[2]
	# log_loss = -true_base + log(sum(exp(probs)))
	log_loss = - (sum(true_base * log(prob) + (1 - true_base) * log(1 - prob))) / length(true_base)
	return(log_loss)
}

loss_prime<-function(theta){
	a = 1/theta
	a = -1/(1-theta)
	return(a)
}

objective<-function(W,X,theta,samp_row){
	# this is the surrogate objective function
	# function takes as arguments: 
	# Weights (W)
	# row being tested (row)
	# dataset (X)
	# theta vector
	# function returns the argument (g) that maximize the expression: 
	library(gtools)
	gs <-permutations(2,3,c(-1,1),repeats.allowed=T)
	for(row in c(1:nrow(gs))) {
		g<-gs[row,]
		x<-X[samp_row,c(1:ncol(X)-1)]
		first_term = t(g) %*% W %*% x
		# u_p <-update_theta(theta)
		second_term = loss(X,samp_row,theta,g)
		# func<-first_term+second_term
		func<-second_term
		# cat(paste(first_term,second_term,"\n"))
		# cat(g)
		# cat("\n")
		if(exists("argmax")){
			if(func<argmax){
				argmax<-func
				argmax_row<-row
			}

		}	else{
			argmax<-func
			argmax_row<-row
			}
	}

	return(gs[argmax_row,])
}

update_theta<-function(theta){
	# function takes as arguments: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	# an 1 by k length matrix of probabilities
	# function returns: 
	# an m+1 by k matrix of softmax processed probabilities of p(y= l|j)	
	for (i in seq(1,nrow(theta))){
		theta_vec<-theta[i,]
		denominator<-do.call(sum,lapply(theta_vec,exp))
		new_probs<-matrix(unlist(lapply(theta_vec, function(x) exp(x)/denominator)))
		theta[i,]<-new_probs
	}
	return(theta)
}



theta<-matrix(c(
				0.8,0.2,
				0.1,0.9,
				0.7,0.3,
				0.4,0.6
				),nrow=4,ncol=2,byrow=TRUE)
w = c(-3,0.5,.5)
# w = c(-0.6,0.15,.5) OPTIMAL
col = ncol(X)-1
W = rep(w,col)
W = matrix(W,nrow=length(w),ncol=col)

tau = 1000
batch = 3
alpha = 0.1 # learning rate
v = 02 # regularization parameter
for(t in seq(0,tau)){
	samp_row <-sample(1:nrow(X),1)
	# current path
	h = sgn(W,X,samp_row) 
	
	g = objective(W,X,theta,samp_row)
	W = W+ (alpha*(g-h)%*%t(X[samp_row,0:col]))
	# if g gives the worst path then use:
	# W = W+ (alpha*-(g+h)%*%t(X[samp_row,0:col]))


	# subtract the worst and add the current
	# W_temp = W+(alpha*g%*%t(X[samp_row,0:col]))
	# +(alpha*h%*%t(X[samp_row,0:col]))

# for(i in seq(1,nrow(W_temp))) {
# 	a = min(1, v**(1/2) / (sum(W_temp[i,]**2)**(1/2))) %*% W_temp[i,]
# 	W_temp[i,]<-a
# }
	
	# W<- (W[,1]+W[,2])/2
	# W<-rep(W,col)
	# W = matrix(W,nrow=length(w),ncol=col)

# delta_3 <- (-(Y - Y_hat) * sigmoidprime(Z_3))
# djdw2 <- t(A_2) %*% delta_3
# W_2 <- W_2 - scalar * djdw2

# true_base = as.numeric(X[samp_row,ncol(X)])
# probs<-theta[,true_base+1]
# if (true_base==1){
# 	true_probs<-c(0,1)
# } else {
# 	true_probs<-c(1,0)
# }

# error = true_probs - f(g) %*% loss_prime(theta)
# # gradient = error %*% true_probs
# # theta = theta -alpha * as.numeric(gradient)

# r<-grep(1,f(sgn(W,X,samp_row)))
# theta[r,] = theta[r,]+ alpha * error

}
# theta<-update_theta(theta)
# W_temp = W+ (alpha*(-g) %*%t(X[samp_row,0:col]))
# W=W_temp
for(i in seq(1,10)){
	cat(paste(f(sgn(W,X,i))),"\n")
	# cat(W %*% X[i,0:col])
	cat(paste(f(sgn(W,X,i))%*%theta),"\n")
}

W
# W<- (W[,1]+W[,2])/2
# W<-rep(W,col)
# W = matrix(W,nrow=length(w),ncol=col)

# delta<- (-(true_probs - f(g) %*% theta ) %*% t(loss_prime(theta)))
# # djdw2<- t(probs) %*% delta
# theta<- theta + alpha %*% t(delta)
# theta<-update_theta(theta)
# cat(theta)
	

# update_theta<-function(theta){
# 	# function takes as arguments: 
# 	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
# 	# an 1 by k length matrix of probabilities
# 	# function returns: 
# 	# an m+1 by k matrix of softmax processed probabilities of p(y= l|j)	
# 	for (i in seq(1,ncol(theta))){
# 		theta_vec<-theta[,i]
# 		denominator<-do.call(sum,lapply(theta_vec,exp))
# 		new_probs<-matrix(unlist(lapply(theta_vec, function(x) exp(x)/denominator)))
# 		theta[,i]<-new_probs
# 	}
# 	return(theta)
# }