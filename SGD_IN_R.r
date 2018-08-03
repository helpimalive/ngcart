
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

f<-function(h){
	# function takes as arguments:
	# an m-bit vector of potential split decisions (h)
	# function returns: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	theta_vec = c(rep(0,length(h)+1))
	position = length(h)+1
	for(bit in seq(1,length(h),2)){
		if(h[bit]>0){
			position=position
		}
		else{
			position=position/2
		}
	}
	theta_vec[position]=1
	return(theta_vec)
}


update_probs<-function(theta_vec,probs){
	# function takes as arguments: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	# an 1 by k length matrix of probabilities
	# function returns: 
	# an m+1 by 1 matrix of softmax processed probabilities of p(y= l|j)	
	denominator<-do.call(sum,lapply(probs,exp))
	new_probs<-matrix(unlist(lapply(probs, function(x) exp(x)/denominator)))
	new_probs = new_probs %*% t(matrix(unlist(theta_vec)))
	new_probs = new_probs[,grep(1,theta_vec)]
	return(new_probs)
}


loss<-function(X,row,update_probs){
	# function takes as arguments: 
	# the dataset and row being evaluated
	# predicted probabilities from the update_probs softmax
	# function returns: 
	# log loss value
	true_base = as.numeric(X[row,ncol(X)])
	log_loss = -true_base + log(sum(exp(update_probs)))
	return(log_loss)
}

loss_prime<-function(probs){
	a = exp(sum(probs))/do.call(sum,lapply(probs,exp))**2
	return(a)
}

objective<-function(W,X,theta_vec){
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
		x<-X[row,c(1:ncol(X)-1)]
		first_term = t(g) %*% W %*% x
		u_p <-update_probs(theta_vec,probs)
		second_term = loss(X,row,u_p)
		func<-first_term+second_term
		if(exists("argmax")){
			if(func>argmax){
				argmax<-func
				row<-row
			}}
		else{
			argmax<-func
			}
	}
	return(gs[row,])
}

probs<-matrix(c(0.00,1),nrow=1,ncol=2)
theta_vec<-c(0,0,0,1)
w = c(0.5,0.3,0.2)
col = ncol(X)-1
W = rep(w,col)
W = matrix(W,nrow=length(w),ncol=col)

tau = 10
batch = 3
alpha = 0.1 # learning rate
v = 2 # regularization parameter

# for(t in seq(0,tau)){}
samp_row <-sample(1:nrow(X),1)
	h = sgn(W,X,samp_row)
	g = objective(W,X,theta_vec)
	W_temp = W-(alpha*g%*%t(X[samp_row,0:col])+alpha*h%*%t(X[samp_row,0:col]))
	# i'm concerned that I don't have to transpose W before multiplying here. 

for(i in seq(1,nrow(W_temp))) {
	a = rep(min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))),col)
	W_temp[1,]<-a
}

y-y_hat

	

