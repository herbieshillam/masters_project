# gene_count.py

import numpy as np
import pandas as pd

def gene_count(final_links):
    # Counts the number of links for each gene in glasso
    
    # final_links = table of links produced from links()
    
    gene_no_dup = final_links['Link from'].drop_duplicates('first')
    gene_no_dup = gene_no_dup.append(final_links['Link to'].drop_duplicates('first'))
    gene_no_dup = gene_no_dup.drop_duplicates('first')
    gene_no_dup = gene_no_dup.reset_index(drop=True)
    
    counter = pd.DataFrame(data=np.zeros(len(gene_no_dup)))
    
    for i in range(len(gene_no_dup)):
        for j in range(len(final_links['Link from'])):
            for k in range(len(final_links.iloc[0])):
                if gene_no_dup.iloc[i] == final_links.iloc[j][k]:
                    counter.iloc[i] += 1
    gene_no_dup = gene_no_dup.to_frame()
    gene_no_dup.columns = ['Gene ID']
    gene_no_dup['Link quantity'] = counter
    gene_no_dup = gene_no_dup.sort_values(by=['Link quantity'],ascending=False)
    
    return gene_no_dup