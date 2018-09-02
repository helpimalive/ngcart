library(gtools)

sgn<-function(W,data,row){
	# W is the weight matrix
	# x is one row of the the innput matrix
	# sgn is the m-bit m-bit vector of potential split decisions (h)
	x = data[row,c(1:ncol(data)-1)]
	out_mat = t(W) %*% x -1
	sign  = sign(out_mat)
	return(sign)}

f <- function(h){
	# function takes as arguments:
	# an m-bit vector of potential split decisions (h)
	# function returns: 
	# an m+1-length one-hot indicator vector, which is only non-zero at the index of the selected leaf
	h[h==-1]<-0
	l<-c(h[1])
	i=1
	while(i<=length(h)){
		n_i<- 2*i+h[i]
		i<-n_i
		l<-append(l,list(h[i]))
		}
	
	j<-0
	while(2**j<= i){
		j=j+1
		}
	out_pos<-  i - 2**(j-1) +1 

	full_vec <- rep(0, length(h)+1)
	full_vec[out_pos] <- 1

	return(full_vec)
}


loss<-function(data,samp_row,theta,g){
	# function takes as arguments: 
	# the dataset and row being evaluated
	# predicted probabilities from the update_theta softmax
	# function returns: 
	# log loss value
	true_base = as.numeric(data[samp_row,ncol(data)])
	probs<- t(theta) %*%f(g)
	prob<-probs[2]
	# cat(f(g)," ",prob,"\n")
	log_loss = - (sum(true_base * log(prob) + (1 - true_base) * log(1 - prob))) / length(true_base)
	return(log_loss)
}

loss_prime<-function(theta){
	a = 1/theta
	# a = -1/(1-theta)
	return(a)
}

objective<-function(W,data,theta,samp_row){
	# this is the surrogate objective function
	# function takes as arguments: 
	# Weights (W)
	# row being tested (row)
	# dataset (data)
	# theta vector
	# function returns the argument (g) that maximize the expression: 
	
	width<-dim(W)[2]
	argmax<- 5e10
	memo<-list()
	gs <-permutations(2,width,c(-1,1),repeats.allowed=T)
	for(row in c(1:nrow(gs))) {
		g<-gs[row,]
		if(!(list(f(g)) %in% memo)){ 		
			x<-data[samp_row,c(1:ncol(data)-1)]
			first_term = g %*% t(W) %*% x
			u_p <-update_theta(theta)
			second_term = loss(data,samp_row,u_p,g)
			# func<-first_term+second_term
			func<-second_term
			# cat(paste(row,first_term,second_term,"\n"))
			# cat(g)
			# cat("\n")
			if(func<argmax){
				argmax<-func
				argmax_row<-row
				}
			memo<-append(memo,list(f(g)))
		}
	}

	return(gs[argmax_row,])
}

objective_verbatum<-function(W,data,theta,samp_row){
	# this is the surrogate objective function
	# function takes as arguments: 
	# Weights (W)
	# row being tested (row)
	# dataset (data)
	# theta vector
	# function returns the argument (g) that maximize the expression: 
	width<-dim(W)[2]
	memo<-list()
	gs <-permutations(2,width,c(-1,1),repeats.allowed=T)
	second_max<-0
	first_max<-0
	for(row in c(1:nrow(gs))) {
		g<-gs[row,]
		if(!(list(f(g)) %in% memo)){ 		
			x<-data[samp_row,c(1:ncol(data)-1)]
			first_term = sum((sgn(W,data,samp_row) - g)^2)
			u_p <-update_theta(theta)
			second_term = loss(data,samp_row,u_p,g)
			val <- first_term+second_term
				if(second_term>=second_max & first_term>first_max){
					argmax_row<-row
					argmax<-val
					second_max<-second_term
					first_max<-first_term
				} 
			memo<-append(memo,list(f(g)))
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
	# theta<-t(theta)
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

total_loss<-function(theta,W,test_data){
	test_theta<-update_theta(theta)
	total_loss<-0
	for(i in seq(1,nrow(test_data))){
		true_base = as.numeric(test_data[i,ncol(test_data)])
		total_loss = total_loss + 1 - (f(sgn(W,test_data,i))%*%(test_theta))[true_base+1]
		# cat(true_base," ",f(sgn(W,data,i))," ",(f(sgn(W,data,i))%*%(test_theta))[true_base+1],"\n")
	}
		return(total_loss)
}

accuracy<-function(theta,W,test_data){
	test_theta<-update_theta(theta)
	acc<-0
	for(i in seq(1,nrow(test_data))){
		true_base = as.numeric(test_data[i,ncol(test_data)])
		predicted = grep(1,round((f(sgn(W,test_data,i))%*%(test_theta))))-1
		if(true_base==predicted){
			acc = acc+1
		}
	}
	return(acc/nrow(test_data))
}

greedy<-function(theta,W,tau,alpha,v,train_data){
	cols<-dim(train_data)[2]-1
	for(t in seq(0,tau)){
		samp_row <-sample(1:nrow(train_data),1)
		# current path
		h = sgn(W,train_data,samp_row) 
		# cat("current path ",h,"\n")

		# optimal path based on cost function (not on the other term)
		g = objective(W,train_data,theta,samp_row)
		# cat("optimal path ",g,"\n")

		W = W + t((alpha*(g-h)%*%t(train_data[samp_row,0:cols])))

		for(i in seq(1,nrow(W))) {
			a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
			W[i,]<-a
			}

		true_base = as.numeric(train_data[samp_row,ncol(train_data)])
		r<-grep(1,f(sgn(W,train_data,samp_row)))
		probs<-theta[r,]

		if (true_base==1){
			true_probs<-c(0,1)
			} else {
			true_probs<-c(1,0)
			}
		
		# cat("current prob of correct ID ",probs[true_base+1],"\n")

		# Gradient Version with Partial Derivative
		gradient = loss_prime(probs)
		theta[r,] = theta[r,] - alpha * (update_theta(t(gradient))-true_probs)
		# cat("udpated prob of correct ID ",theta[r,true_base+1],"\n")
	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
	}


non_greedy<-function(theta,W,tau,alpha,v,train_data){
	cols<-dim(train_data)[2]-1
	for(t in seq(0,tau)){
		samp_row <-sample(1:nrow(train_data),1)
		# current path
		h = sgn(W,train_data,samp_row) 

		# using the paper's loss-augmented inference
		g_v = objective_verbatum(W,train_data,theta,samp_row)
		g_o = objective(W,train_data,theta,samp_row)
		W = W 
			+ t(alpha*(g_o-h)%*%t(train_data[samp_row,0:cols]))
			- t(alpha*(g_v-h)%*%t(train_data[samp_row,0:cols]))
		# g-h is worst
		# h-g is best
		# g = objective(W,train_data,theta,samp_row)
		# W = W+ (alpha*(g-h)%*%t(train_data[samp_row,0:col]))

	for(i in seq(1,nrow(W))) {
		a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
		W[i,]<-a
		}

	true_base = as.numeric(train_data[samp_row,ncol(train_data)])
	r<-grep(1,f(sgn(W,train_data,samp_row)))
	probs<-theta[r,]

	if (true_base==1){
		true_probs<-c(0,1)
		} else {
		true_probs<-c(1,0)
		}

	# Gradient Version with Partial Derivative
	gradient = loss_prime(probs)
	theta[r,] = theta[r,] - alpha * (update_theta(t(gradient))-true_probs)

	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
}

initialize_theta<-function(data,depth){
#	function takes: 
# 		data the train data
# 		an optional depth which is the number of leaves (m+1) or a specified depth if 
# 			not using the full depth of tree. Depth must be of 2**n value (1,2,4,8...)

	if(missing(depth)){
		depth<-0
		i<-0
		while(i<ncol(data)-1){
			depth=depth+2**i
			i=i+1
		}
		depth<-depth+1
	}
	theta<- matrix(c(runif(2*depth**2)),nrow=depth,ncol=2,byrow=T)
	theta<-update_theta(theta)
	return(theta)
}

# theta<-matrix(c(
# 				0.65,0.35,
# 				0.55,0.45,
# 				0.45,0.55,
# 				0.35,0.65
# 				),nrow=4,ncol=2,byrow=TRUE)

initialize_weights<-function(data,depth){
# 	initialize a matrix of m (nodes) by data-1 cols
	if(missing(depth)){
		depth<-1
		i<-1
		while(i<ncol(data)-1){
			depth =depth+2**i
			i=i+1
		}
	}
	W<-matrix(c(runif(depth*(ncol(data)-1))),nrow=depth,ncol=(ncol(data)-1),byrow=T)
	return(t(W))
}


X<-read.csv2('C:\\users\\matt\\desktop\\banknotes.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
X<-X[,c(1,2,5)]
X<-as.matrix(X)

# data = matrix(
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


# w = c(3,-2,3)
# col = ncol(data)-1
# W = rep(w,col)
# W = matrix(W,nrow=length(w),ncol=col)


train_index<-sample(nrow(X),nrow(X)*.8)
train_data<-X[train_index,]
test_data<-X[!(seq(1,nrow(X)) %in% train_index),]

W<-initialize_weights(train_data)
theta<-initialize_theta(train_data)
alpha = 0.05 # learning rate
v = 02 # regularization parameter

results<-data.frame(g_tau=integer(),greedy_acc=double(),ng_tau=double(),ng_acc=double())
for(tau in seq(0,1000,100)){
	W<-initialize_weights(train_data)
	theta<-initialize_theta(train_data)
	out<-greedy(theta,W,tau,alpha,v,train_data)
	g_acc<- accuracy(out$theta,out$W,test_data)
	g_loss<- total_loss(out$theta,out$W,test_data)
	g_tau<- tau
	out<-non_greedy(out$theta,out$W,tau,alpha,v,train_data)
	ng_acc<- accuracy(out$theta,out$W,test_data)
	ng_tau<- 2*tau
	results<-rbind(results,c(g_tau,g_acc,ng_tau,ng_acc))
	cat(
		g_tau
		," "
		,g_acc
		," "
		,ng_tau
		," "
		,ng_acc
		,"\n")
}
names(results)<-c('g_tau','greedy_acc','ng_tau','ng_acc')
res_a<-results[c('g_tau','greedy_acc')]
names(res_a)<-c('tau','greedy_acc')
res_b<-results[c('ng_tau','ng_acc')]
names(res_b)<-c('tau','non_greedy_acc')
res_c<- merge(res_a,res_b)
sum(res_c$non_greedy_acc>=res_c$greedy_acc)/nrow(res_c)

# test_theta<-update_theta(theta)
# test_theta<-update_theta(theta)
# for(i in seq(1,10)){
# 	cat("\n","predicted leaf =",f(sgn(W,data,i)),"\n")
# 	# cat(W %*% data[i,0:col])
# 	cat("probabilities at leaf= ",f(sgn(W,data,i))%*%test_theta,"\n")
# 	true_base = as.numeric(data[i,ncol(data)])
# 	cat(" probability of correct assignment = ",(f(sgn(W,data,i))%*%test_theta)[true_base+1],"\n")
# }