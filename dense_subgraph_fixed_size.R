#the temperature function in the simulated annealing algorithm. 
#Here x is number of the current cooling steps, xmax is the total number of cooling steps, T0 is the initial temperature and Tn is the final temperature.

library(igraph)

temperature<-function(x,xmax,T0,Tn){
	T<-T0*(Tn/T0)^(x/xmax)
	T
}

energy_edge_density<-function(s,adj){
	if(length(s)==1){
		e<-0
	}
	else{
		n<-length(s)
		link_new<-adj[s,s]
		e<-sum(link_new)/(n*(n-1))
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


#G is the input graph (an igraph object), M is the number of iteration at each temperature, alpha is the ratio of global move
dense_subgraph_fixed_size<-function(G,k,M=100,alpha=0.1,T0=1,Tn=0.001,xmax=10){	
	adj<-get.adjacency(G,sparse=F)
	v<-c(1:length(V(G)))
	neighbor<-neighborhood(G,order=1)
	nei_matrix<-matrix(-1,length(v),max(degree(G)+1))
	for(i in 1:length(v)){
		nei_matrix[i,1:length(neighbor[[i]])]<-neighbor[[i]]
		}
	s0<-global_move(k,nei_matrix) ##start subgraph 
	e<-energy_edge_density(s0,adj)
	trace<-c()
	pointer<-c()	
	x=1
	while(x<=xmax){
		t1<-temperature(x,xmax,T0,Tn) ##calculate current temperature
		for(r in 1:M){
			ind<-sample(c(0,1),1,prob=c(alpha,1-alpha))
			if(ind){ #localmove
				vector_id<-as.vector(nei_matrix[s0,])
				vector_id<-vector_id[vector_id>0]
				vector_id<-unique(setdiff(vector_id,s0))
				new_node<-sample(vector_id,1) #add one node from the neighbors
				s<-c(s0,new_node) 
				arti_points<-articulation.points(graph.adjacency(adj[s,s]))
				if(length(arti_points)==0){
					int<-sample(c(1:(k+1)),1)
				}
				else{
					int<-sample(c(1:(k+1))[-arti_points],1)
				}
				s_new<-s[-int] #delete one node that is not a cut vertex
				e_new<-energy_edge_density(s_new,adj)
				p<-min(c(1,exp((e-e_new)/t1)))			
				if(p>runif(1)){
					s0<-s_new
					e<-e_new
				}
				pointer<-c(pointer,length(trace))
				trace<-c(trace,e,s0)
			}
			else{#globalmove
				s_new<-global_move(k,nei_matrix)
				e_new<-energy_edge_density(s_new,adj)
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
	
	e<-c()
	for(i in 1:length(pointer)){
		e<-c(e,trace[pointer[i]+1])
	}
	opt<-which.max(-e)
	member<-trace[(pointer[opt]+2):(pointer[opt]+k+1)]
	Results=list(MaximumDensity=max(-e),Membership=member)
	return(Results)
}

##example
G<-erdos.renyi.game(1000,0.05)
dense_subgraph_fixed_size(G,40)