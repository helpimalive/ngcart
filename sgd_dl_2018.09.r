library(gtools)
library(rpart)
library(tree)

########################
### HELPER FUNCTIONS ###
########################


initialize_theta<-function(data,depth){
	#	function takes: 
	# 		the train data
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
	W<-matrix(c(runif(depth*(ncol(data)-1),-10,10)),nrow=depth,ncol=(ncol(data)-1),byrow=T)
	return(t(W))
}


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

#########################
### TESTING FUNCTIONS ###
#########################

rpart_pred<-function(train_data,test_data,which_cols){
	train_data<-as.data.frame(train_data)
	test_data<-as.data.frame(test_data)
	factors<-which_cols
	formula<- as.formula(paste("base~",paste(factors,collapse="+")))
	fit<- rpart(formula,data = train_data, method="class")
	preds<-as.data.frame(predict(fit,newdata=test_data,type=c("class")))
	names(preds)<-c('pred_val')
	acc<-test_data$base==preds$pred_val
	return(sum(acc)/length(acc))
}

tree_test<-function(train_data,test_data,which_cols){
	train_data<-as.data.frame(train_data)
	test_data<-as.data.frame(test_data)
	factors<-which_cols
	formula<- as.formula(paste("as.factor(base)~",paste(factors,collapse="+")))
	fit<- tree(formula,data = train_data)
	preds<-predict(fit, newdata = test_data,type = "class")
	acc<- sum(preds==test_data[,'base'])/length(preds)
	return(acc)
}

total_loss<-function(theta,W,test_data){
	test_theta<-update_theta(theta)
	total_loss<-0
	for(i in seq(1,nrow(test_data))){
		true_base = as.numeric(test_data[i,ncol(test_data)])
		total_loss = total_loss + 1 - (f(sgn(W,test_data,i))%*%(test_theta))[true_base+1]
	}
		return(total_loss)
}

accuracy<-function(theta,W,test_data){
	test_theta<-update_theta(theta)
	acc<-0
	for(i in seq(1,nrow(test_data))){
		true_base = as.numeric(test_data[i,ncol(test_data)])
		predicted = grep(1,round((f(sgn(W,test_data,i))%*%(test_theta))))-1
		# cat(true_base," ",predicted,"\n")
		if(true_base==predicted){
			acc = acc+1
		}
	}
	return(acc/nrow(test_data))
}

est<-function(theta,W,data){
	test_theta<-update_theta(theta)
	pred<-data.frame()
	for(i in seq(1,nrow(test_data))){
		predicted = grep(1,round((f(sgn(W,test_data,i))%*%(test_theta))))-1
		pred<-rbind(pred,predicted)
	}
	return(pred)
}


#####################
### TREE BUILDERS ###
#####################

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

objective<-function(W,data,theta,samp_row,exhaustive){
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
	# gs <-permutations(2,width,c(-1,1),repeats.allowed=T)
	gs<-matrix(rep(0,width**2),nrow=width,ncol=width)
	for(row in seq(0,width)){
		gs[row,row]<-1
	}

	for(row in c(1:nrow(gs))) {
		
		# if(exhaustive==FALSE){	
		# 	samp_row <-sample(1:nrow(gs),1)
		# 	g<-gs[samp_row,]
		# 	if(!(list(f(g)) %in% memo)){ 		
		# 		x<-data[samp_row,c(1:ncol(data)-1)]
		# 		first_term = g %*% t(W) %*% x
		# 		u_p <-update_theta(theta)
		# 		second_term = loss(data,samp_row,u_p,g)
		# 		if(second_term<argmax){
		# 			argmax<-second_term
		# 			argmax_row<-row
		# 			}
		# 		memo<-append(memo,list(f(g)))
		# 	}
		# } else {
			g<-gs[row,]	
			# if(!(list(f(g)) %in% memo)){ 		
			x<-data[samp_row,c(1:ncol(data)-1)]
			first_term = g %*% t(W) %*% x
			u_p <-update_theta(theta)
			second_term = loss(data,samp_row,u_p,g)
			if(second_term<argmax){
				argmax<-second_term
				argmax_row<-row
				}
			}
		# }
	# }

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
	# gs <-permutations(2,width,c(-1,1),repeats.allowed=T)
	gs<-matrix(rep(0,width**2),nrow=width,ncol=width)
	for(row in seq(0,width)){
		gs[row,row]<-1
	}
	second_max<-0
	first_max<-0
	argmax<- -1e10
	for(row in c(1:nrow(gs))) {
		g<-gs[row,]
		# if(!(list(f(g)) %in% memo)){ 		
			x<-data[samp_row,c(1:ncol(data)-1)]
			first_term = sum((sgn(W,data,samp_row) - g)^2)
			first_term = g%*%t(W)%*%x
			u_p <-update_theta(theta)
			second_term = loss(data,samp_row,u_p,g)
			# cat("\n",first_term,second_term)
			val <- second_term
				if(
					# second_term>=second_max 
					# & first_term=>first_max
					val>argmax
					){
				# if(val>argmax){
					argmax_row<-row
					argmax<-val
					second_max<-second_term
					first_max<-first_term
				} 
			# memo<-append(memo,list(f(g)))
		# }
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

greedy<-function(theta,W,tau,alpha,v,train_data,exhaustive){
	cols<-dim(train_data)[2]-1
	for(t in seq(1,tau)){
		# samp_row <-sample(1:nrow(train_data),1)
		samp_row<-t
		
		# current path
		h = sgn(W,train_data,samp_row) 
		# optimal path based on cost function (not on the other term)
		g = objective(W,train_data,theta,samp_row,exhaustive)
		W = W + t((alpha*(g-h)%*%t(train_data[samp_row,0:cols])))

		for(i in seq(1,nrow(W))) {
			a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
			W[i,]<-a
			}

		h = sgn(W,train_data,samp_row) 
		# if(all(g==h)){
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
		theta[r,] = theta[r,] - alpha * (update_theta(t(gradient))*(theta[r,]-true_probs))
		theta[is.nan(theta)]=0.01
		# }
	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
}


non_greedy<-function(theta,W,tau,alpha,v,train_data){
	cols<-dim(train_data)[2]-1
	for(t in seq(1,tau)){
		samp_row <-sample(1:nrow(train_data),1)
		# samp_row<-t

		# current path
		h = sgn(W,train_data,samp_row) 

		# using the paper's loss-augmented inference
		g_v = objective_verbatum(W,train_data,theta,samp_row)
		g_o = objective(W,train_data,theta,samp_row)
		W = (W 
			# + t(alpha*(h)%*%t(train_data[samp_row,0:cols]))
			+ t(alpha*(g_o-h)%*%t(train_data[samp_row,0:cols]))
			- t(alpha*(g_v)%*%t(train_data[samp_row,0:cols]))
			# - t(alpha*10*(g_v-g_o)%*%t(train_data[samp_row,0:cols]))
			)
		h = sgn(W,train_data,samp_row)

	
	for(i in seq(1,nrow(W))) {
		a = min(1, v**(1/2) / (sum(W[i,]**2)**(1/2))) %*% W[i,]
		W[i,]<-a
		}
	
	# if(all(h==g_o)){
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
		theta[r,] = theta[r,] - alpha * (update_theta(t(gradient))*(theta[r,]-true_probs))
		theta[is.nan(theta)]=0.01
		# }
	}
	ret <-list("theta" = theta, "W"=W)
	return(ret)
}

#############
##BANKNOTES##
#############
# X<-read.csv('C:\\users\\matth\\Documents\\banknotes.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-read.csv('C:\\users\\mlarriva\\desktop\\banknotes.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-read.csv('C:\\users\\Matt\\desktop\\banknotes.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-as.matrix(X)
# names(X)<-c('a','b','c','d','base')
# results<-data.frame(method=character(),accuracy=double())
# which_cols<-c("a","b","c","d")
# X<-X[,c("a","b","c","d","base")]

####
## CODE BREAKS DOWN IN CASE OF SPARSE OR NON_CONTINUOUS DATA
####

###############
## CANCER    ##
###############
X<-read.csv('C:\\users\\mlarriva\\desktop\\cancer_data.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-read.csv('C:\\users\\matth\\desktop\\cancer_data.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")

X<-as.matrix(X)
names(X)<-c('a','b','c','d','base')
results<-data.frame(greedy_acc=double(),ng_acc=double(),rpart=double())
results<-rbind(results,c('greedy','non_greedy','rpart'))
results<-results[-1,]
which_cols<-c("a","b","c","d","e")
X<-X[,c(which_cols,"base")]

#################
## IRIS		   ##
#################
# X<-read.csv('C:\\users\\matt\\desktop\\iris.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-read.csv('C:\\users\\matth\\desktop\\iris.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-read.csv('C:\\users\\mlarriva\\desktop\\iris.csv',header=TRUE,sep=",",stringsAsFactors=F, dec=".")
# X<-as.data.frame(X)
# names(X)<-c('a','b','c','d','base')

# results<-data.frame(method=character(),accuracy=double())
# X[X=='Iris-setosa']<-1
# X[X=='Iris-versicolor']<-0

# X<-subset(X,X$base!='Iris-virginica')
# # X[X=='Iris-setosa']<-0
# X<-data.matrix(X)

# which_cols<-c("a","b","c","d")
# # X[,which_cols]<-scale(X[,which_cols])
# X<-X[,c(which_cols,"base")]



depth<-0
i<-0
while(i<ncol(X)-1){
	depth =depth+2**i
	i=i+1
}

for(case in seq(1,20)){
	
	train_index<-sample(nrow(X),round(nrow(X)*0.80))
	test_data<-X[-train_index,]
	train_data<-X[train_index,]

	i_weights<-apply(train_data[,1:(ncol(train_data)-1)],2,sd)
	W<-matrix(rep(c(i_weights),depth),nrow=dim(X)[2]-1,ncol=depth,byrow=F)

	theta<-initialize_theta(train_data)

	alpha<-0.1
	tau = dim(train_data)[1]
	v<-mean(W)

	it<-0
	g_max_acc<-0
	g_max_theta<-0
	g_max_w<-0

	out<-greedy(theta,W,tau,alpha,v,train_data)
	while(it<10){
		alpha<-alpha/2
		out<-greedy(out$theta,out$W,tau,alpha,v,train_data)
		g_acc<- accuracy(out$theta,out$W,train_data)
		if(g_acc>g_max_acc){
				cat("\n","alpha=",1.0000*round(alpha,4),"accuracy=",1.0000*round(g_acc,2))
				g_max_theta<-out$theta
				g_max_w<-out$W
				g_max_acc<-g_acc
				g_max_alpha<-alpha
				}
		it=it+1
	}

	out<-greedy(g_max_theta,g_max_w,tau,g_max_alpha,v,train_data)
	g_acc_1<- accuracy(g_max_theta,g_max_w,train_data)
	cat("\n","alpha=",1.0000*round(g_max_alpha,4),"accuracy=",1.0000*round(g_acc_1,2))
	
	a_cycle<-0
	ng_acc<-0
	ng_acc_max<-g_acc_1
	# ng_acc_max_theta<-out$theta
	# ng_acc_max_W<-out$W
	ng_acc_max_theta<-g_max_theta
	ng_acc_max_W<-g_max_w
	ng_min_alpha<-alpha

	v2<-abs(mean(out$W))
	if(v<0.01){
		v<-0.1
	}
	vs<-c(v,v2)
	# vs<-c(v)
	# THIS LINE BEATS BANKNOTES
	# vs<-c(v)
	# THIS LINE BEATS CANCER
	# vs<-c(v2)
	original_out<-out
	for(try_v in vs){
		out<-original_out
		## I FORGOT THIS LINE BUT IT NEEDS TO BE IN AND RERUN FOR BANKNOTES AND CANCER
		## USUALLY 0.1 but changed for iris
		alpha<-0.1
		while(a_cycle<=10){

			out<-non_greedy(out$theta,out$W,tau,alpha,try_v,train_data)
			ng_acc<- accuracy(out$theta,out$W,train_data)
			cat("\n",
				"alpha=",1.0000*round(alpha,4),
				"v=",1.0000*round(try_v,2),
				"g_acc=",1.0000*round(g_acc_1,2),
				"ng_acc=",1.0000*round(ng_acc,2)
				# "rpart=",1.0000*round(rpart,2),
				# ng_acc>=rpart
				)

			 if(ng_acc>=ng_acc_max ){
				cat("\n","saving this^ as best")
				# out<-non_greedy(out$theta,out$W,tau,alpha,try_v,rbind(validate_data,train_data))
				ng_acc_max_theta<-out$theta
				ng_acc_max_W<-out$W
				ng_acc_max<-ng_acc
				ng_min_alpha<-alpha
				}
			
			alpha<-alpha/2
			a_cycle=a_cycle+1
			}
			a_cycle=0
	}

	rpart<-rpart_pred(train_data,test_data,which_cols)
	cat("\n ",rpart)
	rpart<-data.frame(method="rpart", accuracy=rpart)
	results<-rbind(results,rpart)
	tree_acc<-tree_test(train_data,test_data,which_cols)
	cat("\n ",tree_acc)
	tree_acc<-data.frame(method="tree_acc", accuracy=as.numeric(tree_acc))
	results<-rbind(results,tree_acc)
	ng_acc<-accuracy(ng_acc_max_theta,ng_acc_max_W,test_data)
	cat("\n ","ng_acc = ",ng_acc)
	ng_acc<-data.frame(method="ng_acc", accuracy=ng_acc)
	results<-rbind(results,ng_acc)
	
}

aggregate(results,list(results$method),mean)
summary(aov(accuracy~method,data=results))

summary(aov(accuracy~method,data=results))
TukeyHSD(aov(accuracy~method,data=results))
# Iris-versicolor vs Iris-virginica
#      method accuracy
# 1     rpart     0.90
# 2  tree_acc     0.95
# 3    ng_acc     0.95
# 4     rpart     0.95
# 5  tree_acc     1.00
# 6    ng_acc     0.95
# 7     rpart     0.95
# 8  tree_acc     1.00
# 9    ng_acc     0.90
# 10    rpart     0.95
# 11 tree_acc     1.00
# 12   ng_acc     0.95
# 13    rpart     1.00
# 14 tree_acc     1.00
# 15   ng_acc     0.95
# 16    rpart     0.95
# 17 tree_acc     0.85
# 18   ng_acc     0.90
# 19    rpart     0.90
# 20 tree_acc     0.90
# 21   ng_acc     0.95
# 22    rpart     0.95
# 23 tree_acc     1.00
# 24   ng_acc     0.95
# 25    rpart     0.95
# 26 tree_acc     0.95
# 27   ng_acc     0.90
# 28    rpart     0.90
# 29 tree_acc     0.90
# 30   ng_acc     0.95
# 31    rpart     0.90
# 32 tree_acc     0.90
# 33   ng_acc     0.90
# 34    rpart     0.80
# 35 tree_acc     0.80
# 36   ng_acc     1.00
# 37    rpart     1.00
# 38 tree_acc     1.00
# 39   ng_acc     1.00
# 40    rpart     1.00
# 41 tree_acc     1.00
# 42   ng_acc     1.00
# 43    rpart     1.00
# 44 tree_acc     1.00
# 45   ng_acc     1.00
# 46    rpart     0.90
# 47 tree_acc     0.90
# 48   ng_acc     0.90
# 49    rpart     0.90
# 50 tree_acc     0.90
# 51   ng_acc     0.95
# 52    rpart     0.85
# 53 tree_acc     0.85
# 54   ng_acc     0.80
# 55    rpart     0.90
# 56 tree_acc     0.95
# 57   ng_acc     0.95
# 58    rpart     0.95
# 59 tree_acc     0.95
# 60   ng_acc     0.95

# Iris-setosa Iris-virginica
#      method accuracy
# 1     rpart     1.00
# 2  tree_acc     1.00
# 3    ng_acc     1.00
# 4     rpart     1.00
# 5  tree_acc     1.00
# 6    ng_acc     1.00
# 7     rpart     1.00
# 8  tree_acc     1.00
# 9    ng_acc     1.00
# 10    rpart     1.00
# 11 tree_acc     1.00
# 12   ng_acc     1.00
# 13    rpart     1.00
# 14 tree_acc     1.00
# 15   ng_acc     1.00
# 16    rpart     1.00
# 17 tree_acc     1.00
# 18   ng_acc     1.00
# 19    rpart     1.00
# 20 tree_acc     1.00
# 21   ng_acc     1.00
# 22    rpart     1.00
# 23 tree_acc     1.00
# 24   ng_acc     1.00
# 25    rpart     1.00
# 26 tree_acc     1.00
# 27   ng_acc     1.00
# 28    rpart     1.00
# 29 tree_acc     1.00
# 30   ng_acc     1.00
# 31    rpart     1.00
# 32 tree_acc     1.00
# 33   ng_acc     0.95
# 34    rpart     1.00
# 35 tree_acc     1.00
# 36   ng_acc     1.00
# 37    rpart     1.00
# 38 tree_acc     1.00
# 39   ng_acc     1.00
# 40    rpart     1.00
# 41 tree_acc     1.00
# 42   ng_acc     1.00
# 43    rpart     1.00
# 44 tree_acc     1.00
# 45   ng_acc     1.00
# 46    rpart     1.00
# 47 tree_acc     1.00
# 48   ng_acc     0.95
# 49    rpart     1.00
# 50 tree_acc     1.00
# 51   ng_acc     1.00
# 52    rpart     1.00
# 53 tree_acc     1.00
# 54   ng_acc     1.00
# 55    rpart     1.00
# 56 tree_acc     1.00
# 57   ng_acc     0.95
# 58    rpart     1.00
# 59 tree_acc     1.00
# 60   ng_acc     1.00

# Iris-setosa 	Iris-versicolor
#      method accuracy
# 1     rpart     1.00
# 2  tree_acc     1.00
# 3    ng_acc     1.00
# 4     rpart     1.00
# 5  tree_acc     1.00
# 6    ng_acc     1.00
# 7     rpart     1.00
# 8  tree_acc     1.00
# 9    ng_acc     0.95
# 10    rpart     1.00
# 11 tree_acc     1.00
# 12   ng_acc     1.00
# 13    rpart     1.00
# 14 tree_acc     1.00
# 15   ng_acc     0.95
# 16    rpart     1.00
# 17 tree_acc     1.00
# 18   ng_acc     1.00
# 19    rpart     1.00
# 20 tree_acc     1.00
# 21   ng_acc     0.95
# 22    rpart     1.00
# 23 tree_acc     1.00
# 24   ng_acc     1.00
# 25    rpart     1.00
# 26 tree_acc     1.00
# 27   ng_acc     1.00
# 28    rpart     1.00
# 29 tree_acc     1.00
# 30   ng_acc     1.00
# 31    rpart     1.00
# 32 tree_acc     1.00
# 33   ng_acc     1.00
# 34    rpart     1.00
# 35 tree_acc     1.00
# 36   ng_acc     1.00
# 37    rpart     1.00
# 38 tree_acc     1.00
# 39   ng_acc     1.00
# 40    rpart     1.00
# 41 tree_acc     1.00
# 42   ng_acc     1.00
# 43    rpart     1.00
# 44 tree_acc     1.00
# 45   ng_acc     1.00
# 46    rpart     1.00
# 47 tree_acc     1.00
# 48   ng_acc     1.00
# 49    rpart     1.00
# 50 tree_acc     1.00
# 51   ng_acc     1.00
# 52    rpart     1.00
# 53 tree_acc     1.00
# 54   ng_acc     1.00
# 55    rpart     1.00
# 56 tree_acc     1.00
# 57   ng_acc     1.00
# 58    rpart     1.00
# 59 tree_acc     1.00
# 60   ng_acc     1.00

# banknotes
# rpart	0.9781022
# tree	0.9708029
# ng_cart	0.9562044
# rpart	0.9708029
# tree	0.9854015
# ng_cart	0.9890511
# rpart	0.9890511
# tree	0.9817518
# ng_cart	0.9927007
# rpart	0.9562044
# tree	0.9890511
# ng_cart	0.9927007
# rpart	0.9379562
# tree	0.9890511
# ng_cart	1
# rpart	0.9306569
# tree	0.9635036
# ng_cart	0.9963504
# rpart	0.9489051
# tree	0.9744526
# ng_cart	0.9927007
# rpart	0.9635036
# tree	0.9781022
# ng_cart	0.9890511
# rpart	0.9160584
# tree	0.9708029
# ng_cart	0.9817518
# rpart	0.9781022
# tree	0.9635036
# ng_cart	0.9744526
# rpart	0.9817518
# tree	0.9927007
# ng_cart	0.9927007
# rpart	0.9671533
# tree	0.9890511
# ng_cart	1
# rpart	0.959854
# tree	0.9744526
# ng_cart	0.9854015
# rpart	0.9525547
# tree	0.9781022
# ng_cart	0.9744526
# rpart	0.9671533
# tree	0.9781022
# ng_cart	0.9927007
