#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
int read_in(NumericMatrix x) {
  int nr = x.nrow();
  int nc = x.ncol();
  std::vector<std::vector<double> >vec(nc);
  for( int i =0; i <nc; i++){
    NumericMatrix::Column col = x(_,i);
    vec[i].assign( col.begin(), col.end());
  }
  NumericMatrix y(nr,nc);
  for( int c=0; c<nc; c++){
    for( int r=0; r<nr; r++){
      y(r,c)=vec[r][c];
    }
  }
  int xx = y(0,1);
  return xx;
}

// [[Rcpp::export]]
float gini_index(NumericVector y, NumericVector g){
  float gini = 0.0;
  NumericVector classes = unique(y);
  int class_count = classes.size();
  int n_instances = y.size();

  float p;
  for(int group=0; group<class_count; group++){
    float size = 0;
    for(int i=0; i<g.size();i++){
      if(g[i]==classes[group]){
        size++;
      }
    }
    if(size==0){
      p=0;
    }
    float score =0.0;
    if(size!=0){
      for(int class_val=0; class_val<classes.size(); class_val++){
        float correctly_assigned = 0;
        for(int i=0; i<g.size();i++){
          if(g[i]==classes[group] && y[i]==classes[class_val])
            correctly_assigned++;
        }
        p = correctly_assigned/size;
        score = score +p*p;
      }
      gini = gini+ (1.0-score)*(size / n_instances);
    }
  }
  return gini;
}

// [[Rcpp::export]]
NumericVector test_split(double value, NumericVector feature){
  int feature_length = feature.size();
  NumericVector grouping;

  for (int i=0; i<feature_length; i++){
    if(feature[i]<value){
      grouping.push_back(1);
    }
    else{
      grouping.push_back(0);
    }
  }

  return(grouping);
}

// [[Rcpp::export]]
List get_split(NumericMatrix m){
  NumericVector  y_vals = m(_,m.ncol()-1);
  NumericVector class_values = unique(y_vals);
  int b_index = 999;
  double b_value=999.0;
  double b_score=999.0;
  NumericVector b_groups;

  for(int col=0; col<m.ncol()-1; col++){
    NumericVector  column = m(_,col);
    for(int row =0; row<m.nrow(); row++){
      NumericVector grouping = test_split(m(row,col),column);
      float gini = gini_index(y_vals,grouping);
        std::cout<<gini<<" "<<row<<" "<<m(row,col)<<std::endl;
        std::cout<<grouping<<std::endl;
        if(gini<b_score){
          b_index = row;
          b_value = m(row,col);
          b_score = gini;
          b_groups = grouping;
        }
    }
}
List ret;
ret["index"] = b_index;
ret["value"] = b_value;
ret["score"] = b_score;
ret["groups"] = b_groups;
ret["left"] = 0;
ret["right"]= 0;
return ret;
}

// [[Rcpp::export]]
NumericVector to_terminal(NumericMatrix m){
  NumericVector  y_vals = m(_,m.ncol()-1);
  NumericVector outcome = unique(y_vals);
  NumericMatrix counts(outcome.size(),2);
  counts.fill(0);
  for(int i=0; i<outcome.size();i++){
    counts(i,0)=outcome[i];
    for(int n=0; n<y_vals.size(); n++){
      if(y_vals[n]==outcome[i]){
        counts(i,1) = counts(i,1)+1;
      }
    }
  }

  float max_val=max(counts(_,1));
  NumericVector y;


  for(int z=0; z<outcome.size();z++){
    if(counts(z,1)==max_val){
      y.push_back(counts(z,0));
    }
  }

  return y;}

// [[Rcpp::export]]
void split(List node, int max_depth,int min_size, int depth,NumericMatrix data){
  List ret;
  int c = data.ncol();
  NumericVector no = as<NumericVector>(node["groups"]);

  int left_len = sum(no);
  int right_len = no.size()-left_len;
  NumericMatrix left(left_len,c);
  NumericMatrix right(right_len,c);

  //left_matrix //left node
  int j=0;
  for(int i =0; i<no.size(); i++){
    if(no[i]==1){
      left(j,_) = data(i,_);
      j=j+1;
    }
  }
  //right matrix //right node
  j=0;
  for(int i =0; i<no.size(); i++){
    if(no[i]==0){
      right(j,_) = data(i,_);
      j=j+1;
    }
  }
  //SKIPPING delete node[groups]

  //check for a no split
  if ((left_len==0) | (right_len==0)){
    node["index"]  = as<int>(node["index"]);
    node["value"]  = as<float>(node["value"]);
    node["score"]  = as<float>(node["score"]);
    node["groups"] = as<NumericVector>(node["groups"]);
    node["left"]   = to_terminal(data);
    node["right"]  = to_terminal(data);
    std::cout<<"no split; node now "<<std::endl;
    std::cout<< as<int>(node["index"]) <<std::endl;
    return;
  }

  //check for max depth
  if(depth> max_depth){
    node["left"]=to_terminal(left);
    node["right"]=to_terminal(right);
    return;
  }

  //process left child
  if(left(_,0).size()<=min_size){
    node["left"] = to_terminal(left);
  }
  else{
    node["left"]=get_split(left);
    split(node["left"],max_depth,min_size,depth+1,data);
  }
  //process right child
  if(right(_,0).size()<=min_size){
    node["right"] = to_terminal(right);
  }
  else{
    node["right"]=get_split(right);
    split(node["right"],max_depth,min_size,depth+1,data);
  }

}

// [[Rcpp::export]]
List build_tree(NumericMatrix train, int max_depth, int min_size){
  List root = get_split(train);
  split(root, max_depth, min_size, 1, train);
  return root;
}

// [[Rcpp::export]]
int predict(List node, NumericVector row){
  int index  = as<int>(node["index"]);
  float value  = as<float>(node["value"]);
  if(row[index] < value){
    if(node.containsElementNamed("right")){
      switch(TYPEOF(node["right"])){
      case REALSXP:
      case INTSXP:
        return node["right"];
      default:
        return predict(node["right"],row);
      }
    }
  }
  else{
    if(node.containsElementNamed("left")){
      switch(TYPEOF(node["left"])){
      case REALSXP:
      case INTSXP:
        return node["left"];
      default:
        return predict(node["left"],row);
      }
    }
  }
}

