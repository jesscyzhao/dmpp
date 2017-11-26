# This is a funciton implements Direchlet Process Mixure Model simulation 
# Assume a gaussian mixture in R2. 
# Each cluster center, i.e. mean of the Gaussian, comes from a Normal with MU_0, SIGMA_0


# Total number of iterations 
N = 50
# Hyper parameters for the normal distibution that generates the cluster center 
MU_0 = cbind(0, 0)
SIGMA_0 = cbind(1.5, 1.5)
# Hyper parameter for the variance in each individual Gaussian
sigma = 0.3
# Hyper parameters for GEM 
ALPHA = 1.5



stick.breaking.gem = function(mu_0, sigma_0, sigma, alpha, n){
    # Now let's initialize some data structure to keep track of our simulations 
    # Rho is the individual piece in the "stick-breaking" process
    rho = c()
    # rhosum is the sum of rhos we have drawn 
    rhosum = 0
    # csum is the cumulative sum of rhos we have drawn, each time we append rhosum to csum; We add a 0 at the beginning incase that first rho is greater than u
    csum = c(0)
    # z is the cluster assignment for each data points 
    z = c()
    # thisz is this cluster that we draw for this data point, z = c(z, thisz)
    thisz = 0 
    # mu is the list of mu that we generated so far
    mu = c()
    # x is the data simulated 
    x = c()
    # number of clusters we have generated so far   
    num_cluster = 0 
    # number of data point we have simulated so far 
    sample_size = 0
    
    u = 1
    rhosum = 0
    palette(rainbow(10))
    par(mar = rep(2,4))
    layout(matrix(c(1,2), 2, 1, byrow=TRUE),
           heights=c(1,4)
    )
    bar_colors = c() 
    for (i in 1:n){
    
    	u = runif(1)
    	# compare rhosum to u. Think it in terms of unrevealed rhos. We land between 0 and 1 first, and then reveal rhos one at a time, until a rho covers u. 
    	while (rhosum < u){
    		V = rbeta(1, 1, alpha)
    	
    		newrho = (1-rhosum)*V
    		rho = rbind(rho, newrho)
    		rhosum = rhosum + newrho 
    		csum = c(csum, rhosum)
    		
    		bar_colors = c(bar_colors,"grey")	
    		
    		mu = rbind(mu, rnorm(cbind(1,1), mu_0, sigma_0))
    		# The length of csum represent how many cluster generated
    
    	}
    
    	# This step is actually pretty deep. In the finite case, sampling K rhos is equivalent to draw u from uniform(0,1), and see where which rho u lands on. So here the csum we generated is the rhosum we generated corresponding to each rho we draw, meaning, the new rho is always at the top of the stick. By comparing csum to u we are trying to see if u lands on the pho that is at the top of the stick, i.e. the i-th csum indicates i-th rho, which indicates i-th cluster
    
    	thisz = max(which(csum<u))
    	z = c(z, thisz)
    	thismu = mu[thisz, ]
    	newx = rnorm(cbind(1,1), thismu, cbind(sigma, sigma))
    	x = rbind(x, newx)
    	num_cluster = length(csum)
    
    
    	
    	barplot(rbind(rho,1-rhosum),
    	        beside=FALSE,
    	        horiz=TRUE,
    	        col=c(bar_colors, "white"),  # remaining mass in (0,1)
    	        xlim=c(0,1),
    	        width=0.7,
    	        main=bquote(rho~"~GEM("~.(alpha)~")")
    	)
    	
        points(u, 1, pch=25, col="red", bg="red")
    	
    	
    	plot(x,
    	     pch=".",
    	     xlim=c(-5,5),
    	     ylim=c(-5,5),
    	     main=paste("N = ", toString(N),
    	                ", #means = ", toString(length(rho)),
    	                ", #clust = ", toString(length(unique(z))),sep="")
    	)
    	
    	# plot all the instantiated means
    	points(mu,
    	       pch=15,
    	       col="black"
    	)
    	
    	# plot the data points generated from the DPMM thus far
    	points(x,
    	       pch=19,
    	       col=z
    	)
    	
    	
    	
    }
    
    print(csum)
    print(num_cluster)
    
}

stick.breaking.gem(MU_0, SIGMA_0, sigma, ALPHA, N)

