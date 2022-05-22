import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
base = str(sys.argv[2])
p = pd.read_csv(str(sys.argv[1]),sep="\t",index_col=0)
print(p)
plt.figure(figsize=(12, 10))
vmin = -0.6
vmax = 0.6
norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
plt.scatter(p['NPscore'],p['Centrality'], c = p['log2FC'], s= 300, alpha=0.8, cmap='bwr',norm=norm)
plt.xlabel("NPscore",fontsize=15)
plt.ylabel("Centrality",fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15)
for st in ["Prkaca","Gapdh","Lrrk2","Mtor","Pkm","Hif3a","Hk2","Ppargc1a"]:
    plt.text(p['NPscore'][st],p['Centrality'][st],s=st,fontsize=18)
plt.colorbar().set_label("log2FC")
plt.savefig(base+"/Cen_NP_figure.png")
