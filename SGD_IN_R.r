
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

W = matrix(c(0.5,0.3,0.2, 0.1,0.2,0.3),nrow=3,ncol=2)
b_i = matrix(c(-1,-1,-1),nrow=3,ncol=1)

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

tau = 10
batch = 3
# for t in seq(0,tau)
h = sgn(W,X,sample(1:nrow(X),1))


def f(base_class,h):
# tree navigation function f :# Hm → Im+1 that maps an m-bit sequence of split decisions (Hm ≡ {−1, +1}
# m) to an indicator vector that specifies a 1-of-(m + 1) encoding. 
# Such an indicator vector is only non-zero at the index of the selected leaf.
	final_direc = list(h).pop()
	if final_direc==-1:
		return(base_class0))
	else:
		return(base_class1))
