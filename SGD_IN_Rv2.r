
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
	ID_vec = c(rep(0,length(h)+1))
	position = length(h)+1
	for(bit in seq(1,length(h),2)){
		if(h[bit]>0){
			position=position
		}
		else{
			position=position/2
		}
	}
	ID_vec[position]=1
	return(ID_vec)
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


loss<-function(X,row,theta){
	# function takes as arguments: 
	# the dataset and row being evaluated
	# predicted probabilities from the update_probs softmax
	# function returns: 
	# log loss value
	true_base = as.numeric(X[row,ncol(X)])
	log_loss = -true_base + log(sum(exp(theta)))
	return(log_loss)
}

loss_prime<-function(theta,h){
	ID_vec<-f(h)
	probs<-theta[grep(1,ID_vec),]
	a = exp(sum(probs))/do.call(sum,lapply(probs,exp))**2
	return(a)
}

objective<-function(W,X,theta,h){
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
		u_p <-update_probs(theta,f(h))
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

probs<-matrix(c(0.5,0.5),nrow=1,ncol=2)
theta<-matrix(c(0,0,0,0,0,0,0.25,0.75),nrow=4,ncol=2,byrow=TRUE)
ID_vec<-c(0,0,1,0)
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
	g = objective(W,X,theta,h)
	W_temp = W-(alpha*g%*%t(X[samp_row,0:col])+alpha*h%*%t(X[samp_row,0:col]))
	# i'm concerned that I don't have to transpose W before multiplying here. 

for(i in seq(1,nrow(W_temp))) {
	a = rep(min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))),col)
	W_temp[1,]<-a
}

# delta_3 <- (-(Y - Y_hat) * sigmoidprime(Z_3))
# djdw2 <- t(A_2) %*% delta_3
theta<- theta+alpha*loss_prime(theta,h)

	

