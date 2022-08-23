def models(model_type, beta, mu, Nnodes, Tmax ):
    
    infection = float(beta)# force given infection parameter to float
    recovery  = float(mu)  # force given recovery parameter to float
    timemats = temp_adj_mats(Nnodes,Tmax)
    s0 = (1-infection)*np.ones(Nnodes)
    i0 = infection*np.ones(Nnodes)
    r0 = np.zeros(Nnodes)
    ss0 = (1-infection)**2*np.ones([Nnodes,Nnodes])
    si0 = infection*(1-infection)*np.ones([Nnodes,Nnodes])
    assert (infection >= 0 and infection <= 1), "infection probability incorrect"
    assert (recovery >= 0 and recovery <= 1), "recovery probability incorrect"
    

    
    sus_over_time = np.zeros((Nnodes, Tmax+1))
    inf_over_time = np.zeros((Nnodes, Tmax+1))
    rec_over_time = np.zeros((Nnodes, Tmax+1))
    
    sus_over_time[:,0]= s0
    inf_over_time[:,0] = i0

    time_iter = range(Tmax) #time starts at 0
    
    if model_type == 'tpb':
        time_idx = 0
        adj_mat = timemats[time_idx]
        for time in time_iter:
            psi = np.ones(Nnodes)
            phi = np.ones([Nnodes,Nnodes])
            for i in range(Nnodes):
                for j in range(Nnodes):
                        psi[i] = psi[i]*(1-infection*adj_mat[i,j]*si0[i,j]/s0[i])
                        phi[i,j] = psi[i]/(1-infection*adj_mat[i,j]*si0[i,j]/s0[i])
                        ss0[i,j] = phi[i,j]*phi[j,i]*ss0[i,j]
                        si0[i,j] = (1-recovery)*(1-infection*adj_mat[i,j])*phi[i,j]*si0[i,j] + phi[i,j]*(1-phi[j,i])*ss0[i,j]

            
            s0 = psi*s0
            i0 = (1. - recovery)*i0 + (1. - psi)*s0
            r0 = 1 - s0 - i0


            sus_over_time[:, time_idx+1] = s0
            inf_over_time[:, time_idx+1 ] = i0
            rec_over_time[:, time_idx+1 ] = r0

            time_idx = (time_idx + 1) % Tmax
            
        
        if model_type == 'cb':
            time_idx = 0
            adj_mat = timemats[time_idx]
            theta = np.ones([Nnodes,Nnodes])
            for time in time_iter:
                for i in range(Nnodes):
                    for j in range(Nnodes):
                            theta[i,j] -= (1-infection*adj_mat[j,i])
                            
                    s0[i] = theta*s0[i]
                    i0[i] = (1. - recovery)*i0[i] + (1. - psi)*sus_over_time[:,time_idx][i]
                    r0[i] = 1 - s0 - i0


                sus_over_time[:, time_idx+1] = s0
                inf_over_time[:, time_idx+1 ] = i0
                rec_over_time[:, time_idx+1 ] = r0

                time_idx = (time_idx + 1) % Tmax


    return sus_over_time, inf_over_time, rec_over_time 
