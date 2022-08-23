def montecarlosir(infection,recovery,Nnodes,time,Nsim):
    infection = float(infection)
    recovery = float(recovery)
    simsus = np.zeros((Nnodes,time,Nsim))
    siminf = np.zeros((Nnodes,time,Nsim))
    simrec = np.zeros((Nnodes,time,Nsim))
    tempadjs = temp_adj_mats(Nnodes,time)
    susprobmc = np.zeros((Nnodes,time))
    infprobmc = np.zeros((Nnodes,time))
    recprobmc = np.zeros((Nnodes,time))
    for i in range(Nsim):
        susceptible = np.ones(Nnodes)
        infected = np.zeros(Nnodes)
        recovered = np.zeros(Nnodes)
        init_node = np.random.randint(0,Nnodes)
        infected[init_node] = 1
        susceptible[init_node] = 0
        simsus[:,0,i] = susceptible
        siminf[:,0,i] = infected
        for j in range(1,time):
            adj_mat = tempadjs[j]
            # INFECTION PROCESS
            for k in range(Nnodes):
                if infected[k] == 1:
                    for l in range(Nnodes):
                        if (infected[l] == recovered[l] == 0) & bool(adj_mat[k,l])==True :
                            inf_event = np.random.binomial(1,infection)
                            if inf_event == 0:
                                susceptible[l] = 1
                                infected[l] = 0
                            elif inf_event == 1:
                                susceptible[l] = 0
                                infected[l] = 1
                    #RECOVERY PROCESS
                    rec_event = np.random.binomial(1,recovery)
                    if rec_event == 1:
                        infected[k] = 0
                        recovered[k] = 1
                    elif rec_event == 0:
                        infected[k] = 1
                        recovered[k] = 0
                    
            simsus[:,j,i] = susceptible
            siminf[:,j,i] = infected
            simrec[:,j,i] = recovered
            
        for i in range(Nnodes):
            for j in range(time):
                susprobmc[i,j] = np.sum(simsus[i,j,:])/Nsim
                infprobmc[i,j] = np.sum(siminf[i,j,:])/Nsim
                recprobmc[i,j] = np.sum(simrec[i,j,:])/Nsim
    return susprobmc, infprobmc, recprobmc
                
            
