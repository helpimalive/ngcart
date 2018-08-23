
sgn<-function(W,X,row){
	# W is the weight matrix
	# x is one row of the the innput matrix
	# sgn is the m-bit m-bit vector of potential split decisions (h)
	x = X[row,c(1:ncol(X)-1)]
	# x<-as.matrix(x)
	out_mat = W %*% x -1
	sign  = sign(out_mat)
	return(sign)}

f <- function(h){
	# function takes as arguments:
	# an m-bit vector of potential split decisions (h)
	# function returns: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	
	h[h==-1]<-0
	i=1
	l<-c(h[1])

	while(i<length(h)){
		from<-i
		upto<- i*2
		cat(from," ",upto,"\n")
		l<-append(l,list(h[from:upto]))
		i=i*2 +1
	}
	h<- l

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


loss<-function(X,samp_row,theta,g){
	# function takes as arguments: 
	# the dataset and row being evaluated
	# predicted probabilities from the update_theta softmax
	# function returns: 
	# log loss value
	true_base = as.numeric(X[samp_row,ncol(X)])
	probs<-t(theta)%*%f(g)
	prob<-probs[2]
	log_loss = - (sum(true_base * log(prob) + (1 - true_base) * log(1 - prob))) / length(true_base)
	return(log_loss)
}

loss_prime<-function(theta){
	a = 1/theta
	# a = -1/(1-theta)
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
		first_term = g %*% W %*% x
		u_p <-update_theta(theta)
		second_term = loss(X,samp_row,u_p,g)
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

objective_verbatum<-function(W,X,theta,samp_row){
	# this is the surrogate objective function
	# function takes as arguments: 
	# Weights (W)
	# row being tested (row)
	# dataset (X)
	# theta vector
	# function returns the argument (g) that maximize the expression: 
	library(gtools)
	gs <-permutations(2,3,c(-1,1),repeats.allowed=T)
	second_max<-0
	first_max<-0
	for(row in c(1:nrow(gs))) {
		g<-gs[row,]
		x<-X[samp_row,c(1:ncol(X)-1)]

		first_term = sum((sgn(W,X,samp_row) - g)^2)
		u_p <-update_theta(theta)
		second_term = loss(X,samp_row,u_p,g)
		# cat(row," ",first_term," ",second_term," ",g,"\n")
		val <- first_term+second_term
			if(
				second_term>=second_max 
				& 
				first_term>first_max
				){
				argmax_row<-row
				argmax<-val
				second_max<-second_term
				first_max<-first_term
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
	if(!is.null(nrow(theta))){
		for (i in seq(1,nrow(theta))){
			theta_vec<-theta[i,]
			denominator<-do.call(sum,lapply(theta_vec,exp))
			new_probs<-matrix(unlist(lapply(theta_vec, function(x) exp(x)/denominator)))
			theta[i,]<-new_probs}
			} else {
			denominator<-do.call(sum,lapply(theta,exp))
			theta<-matrix(unlist(lapply(theta, function(x) exp(x)/denominator)))	
			}
	return(theta)
}

total_loss<-function(theta,W){
	test_theta<-update_theta(theta)
	total_loss<-0
	for(i in seq(1,10)){
		true_base = as.numeric(X[i,ncol(X)])
		total_loss = total_loss + 1- (f(sgn(W,X,i))%*%test_theta)[true_base+1]
	}
		return(total_loss)
}


greedy<-function(theta,W,tau,alpha,v){
	for(t in seq(0,tau)){
		samp_row <-sample(1:nrow(X),1)
		# current path
		h = sgn(W,X,samp_row) 

		# optimal path based on cost function (not on the other term)
		g = objective(W,X,theta,samp_row)
		W = W+ (alpha*(g-h)%*%t(X[samp_row,0:col]))

	for(i in seq(1,nrow(W))) {
		a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
		W[i,]<-a
		}

	true_base = as.numeric(X[samp_row,ncol(X)])
	r<-grep(1,f(sgn(W,X,samp_row)))
	probs<-theta[r,]

	if (true_base==1){
		true_probs<-c(0,1)
		} else {
		true_probs<-c(1,0)
		}
	# Gradient Version with Partial Derivative
	gradient = loss_prime(probs)
	theta[r,] = theta[r,] - alpha* t(update_theta(gradient))
	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
	}


non_greedy<-function(theta,W,tau,alpha,v){
	for(t in seq(0,tau)){
		samp_row <-sample(1:nrow(X),1)
		# current path
		h = sgn(W,X,samp_row) 

		# using the paper's loss-augmented inference
		g_v = objective_verbatum(W,X,theta,samp_row)
		g_o = objective(W,X,theta,samp_row)
		W = W 
			# + (alpha*(g_o-h)%*%t(X[samp_row,0:col]))
			- (alpha*(g_v-h)%*%t(X[samp_row,0:col]))
		# g-h is worst
		# h-g is best
		# g = objective(W,X,theta,samp_row)
		# W = W+ (alpha*(g-h)%*%t(X[samp_row,0:col]))

	for(i in seq(1,nrow(W))) {
		a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
		W[i,]<-a
		}

	true_base = as.numeric(X[samp_row,ncol(X)])
	r<-grep(1,f(sgn(W,X,samp_row)))
	probs<-theta[r,]

	if (true_base==1){
		true_probs<-c(0,1)
		} else {
		true_probs<-c(1,0)
		}

	# Gradient Version with Partial Derivative
	gradient = loss_prime(probs)
	theta[r,] = theta[r,] - alpha* t(update_theta(gradient))

	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
}


X<-read.csv2('C:\\users\\mlarriva\\desktop\\banknotes.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
X<-as.matrix(X)

# X = matrix(
# 	c(2.771244718,1.784783929,0,
# 		1.728571309,1.169761413,0,
# 		3.678319846,2.81281357,0,
# 		3.961043357,2.61995032,0,
# 		2.999208922,2.209014212,0,
# 		7.497545867,3.162953546,1,
# 		9.00220326,3.339047188,1,
# 		7.444542326,0.476683375,1,
# 		10.12493903,3.234550982,1,
# 		6.642287351,3.319983761,1),
# 	nrow=10,
# 	ncol=3,
# 	byrow=TRUE)

initialize_theta<-function(X,depth){
	if(missing(depth)){
		depth=ncol(X)-1
	}
	theta<- matrix(c(runif(2*depth**2)),nrow=depth**2,ncol=2,byrow=T)
	theta<-update_theta(theta)
	return(theta)
}

# theta<-matrix(c(
# 				0.65,0.35,
# 				0.55,0.45,
# 				0.45,0.55,
# 				0.35,0.65
# 				),nrow=4,ncol=2,byrow=TRUE)

initialize_weights<-function(X,depth){
# 	initialize a matrix of m (nodes) by X-1 cols
	if(missing(depth)){
		depth<-0
		i<-1
		while(i<=ncol(X)-1){
			depth=depth+2**i
			i=i+1
		}
	}
	W<-matrix(c(runif(depth*(ncol(X)-1))),nrow=depth,ncol=(ncol(X)-1),byrow=T)
	return(W)
}

# w = c(3,-2,3)
# col = ncol(X)-1
# W = rep(w,col)
# W = matrix(W,nrow=length(w),ncol=col)
W<-initialize_weights(X)
theta<-initialize_theta(X)
tau = 80
batch = 3
alpha = 0.01 # learning rate
v = 02 # regularization parameter


results<-data.frame(g_tau=integer(),greedy_loss=double(),ng_tau=double(),ng_loss=double())
for(tau in seq(0,100,2)){
	out<-greedy(theta,W,tau,alpha,v)
	g_loss<- total_loss(out$theta,out$W)
	g_tau<- tau
	out<-non_greedy(out$theta,out$W,tau,alpha,v)
	ng_loss<- total_loss(out$theta,out$W)
	ng_tau<- 2*tau
	results<-rbind(results,c(g_tau,g_loss,ng_tau,ng_loss))
	# cat(g_tau," ",g_loss," ",ng_tau," ",ng_loss,"\n")
}
names(results)<-c('g_tau','greedy_loss','ng_tau','ng_loss')
res_a<-results[c('g_tau','greedy_loss')]
names(res_a)<-c('tau','greedy_loss')
res_b<-results[c('ng_tau','ng_loss')]
names(res_b)<-c('tau','non_greedy_loss')
res_c<- merge(res_a,res_b)
sum(res_c$non_greedy_loss<res_c$greedy_loss)/nrow(res_c)
# test_theta<-update_theta(theta)
# for(i in seq(1,10)){
# 	cat("\n","predicted leaf =",f(sgn(W,X,i)),"\n")
# 	# cat(W %*% X[i,0:col])
# 	cat("probabilities at leaf= ",f(sgn(W,X,i))%*%test_theta,"\n")
# 	true_base = as.numeric(X[i,ncol(X)])
# 	cat(" probability of correct assignment = ",(f(sgn(W,X,i))%*%test_theta)[true_base+1],"\n")
# }