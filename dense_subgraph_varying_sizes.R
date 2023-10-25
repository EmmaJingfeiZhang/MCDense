library(igraph)

temperature<-function(x,xmax,T0,Tn){
	T<-T0*(Tn/T0)^(x/xmax)
	T
}

energy_average_degree<-function(s,adj){
	if(length(s)==1){
		e<-0
	}
	else{
		n<-length(s)
		link_new<-adj[s,s]
		e<-sum(link_new)/n
	}
	-e
}	

global_move<-function(k,nei_matrix){
	cluster<-clusters(G)
	s<-sample(c(1:nrow(nei_matrix))[cluster$membership==which.max(cluster$csize)],1)
    for(i in 2:k){
		vector_id<-as.vector(nei_matrix[s,])
		vector_id<-vector_id[vector_id>0]
		vector_id<-unique(setdiff(vector_id,s))
		if(length(vector_id)==1){
			s<-c(s,vector_id)
		}
		else{
		    s<-c(s,sample(vector_id,1))
		}
	}
	s
}



#minsize/maxsize is the lower/upper bound on the sizes, beta1 is the prob for decrease sizes, beta2 is prob for increase sizes, alpha is prob for global moves
dense_subgraph_varying_sizes<-function(G,minsize,maxsize,M=100,beta1=0.25,beta2=0.25,alpha=0.1,T0=10,Tn=0.001,xmax=10){	
	adj<-get.adjacency(G,sparse=F)
	v<-c(1:length(V(G)))
	neighbor<-neighborhood(G,order=1)
	nei_matrix<-matrix(-1,length(v),max(degree(G)+1))
	for(i in 1:length(v)){
		nei_matrix[i,1:length(neighbor[[i]])]<-neighbor[[i]]
		}
	s0<-global_move(minsize,nei_matrix) ##start subgraph 
	e<-energy_average_degree(s0,adj)
	trace<-c()
	pointer<-c()	
	x=1
	while(x<=xmax){
		t1<-temperature(x,xmax,T0,Tn)
		for(r in 1:M){
			k<-length(s0)
			ind<-sample(c(1,2,3),1,prob=c(beta1,beta2,1-beta1-beta2))
			if(ind==1){  #delete a node
				non<-articulation.points(graph.adjacency(adj[s0,s0]))
				if(length(non)==0){
					delete_node<-sample(s0,1)
				}
				else{
					delete_node<-sample(s0[-non],1)
				}
				s_new<-c(s0[s0!=delete_node])
				e_new<-energy_average_degree(s_new,adj)
				p<-min(c(1,exp((e-e_new)/t1))*(k>minsize))
				if(p>runif(1)){
					s0<-s_new
					e<-e_new
				}
				pointer<-c(pointer,length(trace))
				trace<-c(trace,e,s0)
			}
			else if(ind==2){#add a node
				vector_id<-as.vector(nei_matrix[s0,])
				vector_id<-vector_id[vector_id>0]
				vector_id<-unique(setdiff(vector_id,s0))
				replace_node<-sample(vector_id,1)
				s_new<-c(s0,replace_node)
				e_new<-energy_average_degree(s_new,adj)
				p<-min(c(1,exp((e-e_new)/t1))*(k<maxsize))
				if(p>runif(1)){
					s0<-s_new
				    e<-e_new
				}
				pointer<-c(pointer,length(trace))
				trace<-c(trace,e,s0)
			}
			else if(ind==3&runif(1)>alpha){#localmove
				vector_id<-as.vector(nei_matrix[s0,])
				vector_id<-vector_id[vector_id>0]
				vector_id<-unique(setdiff(vector_id,s0))
				new_node<-sample(vector_id,1)
			    s<-c(s0,new_node)
				nono<-articulation.points(graph.adjacency(adj[s,s]))
				if(length(nono)==0){
					int<-sample(c(1:(k+1)),1)
				}
				else{
					int<-sample(c(1:(k+1))[-nono],1)
				}
				s_new<-s[-int]
				e_new<-energy_average_degree(s_new,adj)
				p<-min(c(1,exp((e-e_new)/t1)))			
				if(p>runif(1)){
					s0<-s_new
					e<-e_new
				}
				pointer<-c(pointer,length(trace))
				trace<-c(trace,e,s0)
			}
			else{ #global move
				s_new<-global_move(k,nei_matrix)
				e_new<-energy_average_degree(s_new,adj)
				p<-min(c(1,exp((e-e_new)/t1)))			
				if(p>runif(1)){
					s0<-s_new
					e<-e_new
				}
				pointer<-c(pointer,length(trace))
				trace<-c(trace,e,s0)
			}
		}
		x<-x+1
	}
	density<-c()
	for(i in 1:length(pointer)){
		density<-c(density,trace[pointer[i]+1])
	}
	opt<-which.max(-density)
	member<-trace[(pointer[opt]+2):(pointer[opt+1])]
	Results=list(MaximumDensity=max(-density),Membership=member,Size=length(member))
	return(Results)
}


##example
G<-erdos.renyi.game(1000,0.05)
dense_subgraph_varying_sizes(G,400,500)
