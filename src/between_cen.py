#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import networkx as nx
import numpy as np
from networkx.algorithms import centrality
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys
a = pd.read_csv(str(sys.argv[1]),header=None,sep="\t")
base=str(sys.argv[2])
G = nx.from_pandas_edgelist(a,0,1)
b_cen = centrality.betweenness_centrality(G)
pd.DataFrame.from_records([[key, rank+1,np.exp(b_cen[key]*2)*15] for rank, key in enumerate(sorted(b_cen,key=b_cen.get, reverse=True))],columns=['Gene','rank','size_by_betcen']).to_csv(base+"/size_by_np.csv",sep="\t",index=False)
pd.DataFrame.from_records([[rank+1, key,b_cen[key]] for rank, key in enumerate(sorted(b_cen,key=b_cen.get, reverse=True))],columns=['Rank','Gene','Centrality']).to_csv(base+"/Betweenness_Centrality.tsv",sep="\t",index=False)
