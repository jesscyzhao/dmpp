## Simulate CRP 

# Total number of iterations 
N = 100
# Hyper parameters for the normal distibution that generates the cluster center 
mu_0 = cbind(0, 0)
sigma_0 = cbind(1.5, 1.5)
# Hyper parameter for the variance in each individual Gaussian
sigma = 0.3
# Hyper parameters for GEM 
ALPHA = 2


cpr.generate.data = function(mu_0 = cbind(0, 0), sigma_0=cbind(1.5, 1.5), sigma=0.3, alpha=2){

    table_assingment = list(1)
    table_id = 1
    tables = c(table_id)
    get.seating.prob = function(table_assingment, alpha){
        head_counts = sapply(table_assingment, length) 
        rate = c(head_counts, alpha)/sum(head_counts, alpha)
        return(rate)
    }
    
    for (i in 2:N){
        seating_prob = get.seating.prob(table_assingment, ALPHA)
        table_i = sample(c(tables, table_id+1), 1,  prob=seating_prob)
        if(table_i == (table_id+1)){
            table_id = table_id + 1
            tables = c(tables, table_id)
            table_assingment[[table_i]] = i
        }
        else{
        table_assingment[[table_i]] = c(table_assingment[[table_i]], i) 
        }
    }
    
    mu = lapply(replicate(length(table_assingment), list), function(x){rnorm(cbind(1,1), mu_0, sigma_0)})
    
    simulate.one.cluster = function(i, table_assingment, mu, sigma){
        num_data = length(table_assingment[[i]])
        this_mu = mu[[i]]
        data_points = t(replicate(num_data, rnorm(cbind(1,1), this_mu, sigma)))
        this_mu_matrix = t(matrix(rep(this_mu, num_data), ncol=num_data))
        return(cbind(table_assingment[[i]], rep(i, num_data), this_mu_matrix, data_points))
    }
    
    sim_data_list = lapply(1:length(table_assingment), simulate.one.cluster, table_assingment = table_assingment,mu= mu, sigma=sigma) 
    
    sim_data = do.call(rbind, sim_data_list)
    
    colnames(sim_data) = c("id", "cluster", "mu_x","mu_y", "data_x", "data_y")
    
    return(sim_data[order(sim_data[,"id"]), ])
}

data = cpr.generate.data()
