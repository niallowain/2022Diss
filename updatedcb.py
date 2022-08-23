import numpy as np
import networkx as nx


def integrate_contact_based_model(alpha,beta,edgelist,outbreak_origin = None,Tmax = None,individual_based = True,directed = True,verbose = True ):

    edgelist = np.array(edgelist, dtype=int)
    alpha = float(alpha)
    beta = float(beta)
    verbose = bool(verbose)
    individual_based = bool(individual_based)
    directed = bool(directed)
   
    nodes = set( edgelist[:,1] ) | set( edgelist[:,2] )
    Nnodes = len(nodes)
    assert (nodes == set(range(Nnodes))),"nodes must be named continuously from 0 to 'number of nodes'-1"
    assert (alpha >= 0 and alpha <= 1), "infection probability incorrect"
    assert (beta >= 0 and beta <= 1), "recovery probability incorrect"

    if isinstance(outbreak_origin, int):
    
        assert (outbreak_origin >= 0 and outbreak_origin <  Nnodes), "0 <= outbreak origin < 'number of nodes'"
    out = _integrate_fixed_source(alpha,
    beta,
    outbreak_origin, 
    edgelist,
    Tmax=Tmax,
    verbose=verbose,
    individual_based=individual_based,
    directed=directed)
    
    susceptible, infected, recovered = out

    return susceptible, infected, recovered


def _integrate_fixed_source(alpha,
    beta,
    outbreak_origin,
    edgelist,
    Tmax=None,
    verbose=False,
    individual_based=True,
    directed = True):
    
    edgelist = np.array(edgelist, dtype=int)

    if not directed:
        edgelist_directed = np.zeros(( len(edgelist)*2, 3 ))
        for idx in range( len(edgelist) ):
            t, u, v = edgelist[idx]
            edgelist_directed[2*idx, :] = t, u, v
            edgelist_directed[2*idx + 1, :] = t, v, u
        edgelist = edgelist_directed.astype(int)

    times = np.unique(edgelist[:,0])
    T = max(times) + 1 # last time step + 1
    if Tmax is None:
        Tmax = T
    else:
        Tmax += 1

    G = nx.DiGraph()
    G.add_edges_from([(u,v) for t,u,v in edgelist])
    assert (G.remove_edges_from(list(nx.selfloop_edges(G))), "self loops are not allowed")
    static_out_neighbours = G.succ #for speed
    static_in_neighbours = G.pred #for speed
    Nedges = G.number_of_edges()
    Nnodes = G.number_of_nodes()
    edge_to_idx = {edge: idx for idx, edge in enumerate(G.edges())}
    reciprical_link_exists = {edge_to_idx[ u, v ]: True if G.has_edge(v,u) else False for u,v in G.edges()} 
    active_edges_dct = {t: [] for t in range(Tmax)} #for speed
    temporal_in_neighbours = {time: {} for time in times}
    for time in times:
        snapshot = edgelist[ edgelist[:,0] == time ] # contacts at time instance
        temporal_in_neighbours[ time ] = { target: set() for target in np.unique( snapshot[:,2] ) }
        for _, source, target in edgelist[ edgelist[:,0] == time ]:
            temporal_in_neighbours[time][target].add(source)
        active_edges_dct[time] = np.array( [edge_to_idx[ u, v ] for _, u, v in snapshot ] )
    active_targets = {t: set(temporal_in_neighbours[t].keys()) if t in times else set() for t in range(T)} 
    
    sus_over_time = np.zeros((Nnodes, Tmax+1))
    inf_over_time = np.zeros((Nnodes, Tmax+1))
    rec_over_time = np.zeros((Nnodes, Tmax+1))

    inf_over_time = np.zeros((Nnodes, Tmax+1))

    inf_over_time = np.zeros((Nnodes, Tmax+1))

    sus     = np.ones(Nedges) 
    sus_new = np.ones(Nedges) 
    theta   = np.ones(Nedges) 
    phi     = np.zeros(Nedges)
    
    init_idx = np.array([edge_to_idx[node, cavity] for (node, cavity) in G.edges() if node == outbreak_origin], dtype=int)
    sus[init_idx] = 0.
    sus_new[init_idx] = 0.
    phi[init_idx] = 1.
    
    susceptible = np.ones(Nnodes)
    susceptible[outbreak_origin] = 0.
    infected = np.zeros(Nnodes) 
    infected[outbreak_origin] = 1.
    recovered = np.zeros(Nnodes)
    
    sus_over_time[:,0] = susceptible
    inf_over_time[:,0] = infected
    
    if verbose: 
        try:
            from tqdm import tqdm_notebook as tqdm
        except ImportError:
            print ("package 'tqdm' not found")
            print ("continue with verbose = False...")
            time_iter = range(Tmax)
        else:
            time_iter = tqdm(range(Tmax)) #requires tqdm package installed
    else:
        time_iter = range(Tmax) #time starts at 0
    
    time_idx = 0
    for time in time_iter:
        for target in active_targets[time_idx]:
            if target != outbreak_origin:
                susceptible[target] = 1.
                for source in static_in_neighbours[ target ]:
                    edge_idx = edge_to_idx[ source, target ]
                    if source in temporal_in_neighbours[time_idx][target]:
                        if not individual_based:
                            theta[edge_idx] -= alpha * phi[edge_idx]
                        else:
                            theta[edge_idx] *= 1. - alpha * infected[source]
                    susceptible[target] *= theta[ edge_idx ]

                if not individual_based:
                    for cavity in static_out_neighbours[target]:
                        edge_idx = edge_to_idx[ target, cavity ]
                        sus_new[ edge_idx ] = susceptible[ target ]
                        if reciprical_link_exists[ edge_idx ]:
                            reciprical_edge_idx = edge_to_idx[ cavity, target ]
                            sus_new[ edge_idx ] /= theta[ reciprical_edge_idx ]

        if not individual_based:    
            active_edges = np.zeros(Nedges, dtype=bool)
            active_edges[ active_edges_dct[time_idx] ] = True
            phi *= (1. - beta) * (1. - alpha * active_edges)
            phi +=  sus - sus_new
            sus = sus_new.copy()
        
        recovered += infected * beta
        infected = 1. - recovered - susceptible
        
        sus_over_time[:, time + 1] = susceptible
        inf_over_time[:, time + 1] = infected
        rec_over_time[:, time + 1] = recovered

        time_idx = (time_idx + 1) % T
        
    return sus_over_time, inf_over_time, rec_over_time

