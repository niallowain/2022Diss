def adj_mat_gen(dim):
    adj_mat_temp = np.random.randint(0,2,(dim,dim))
    adj_mat_temp[np.diag_indices_from(adj_mat_temp)] = 0 
    return adj_mat_temp

def temp_adj_mats(dim,time):
    times_adj_mat = np.zeros((time,dim,dim))
    for i in range(time):
        times_adj_mat[i,:,:] = adj_mat_gen(dim)
    return times_adj_mat
