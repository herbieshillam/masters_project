# links.py

import pandas as pd

def links(Theta,data): 
    # Produces DataFrame of links produced using glasso
    
    # NOTE: data must be in DataFrame form, first col: gene ID, last col: p-value/kurtosis
    
    result = [] # initialise empty list
    for i in range(len(Theta)):
        for j in range(len(Theta)):
            if i < j and abs(Theta[i][j]) > 0.0001 and data.iloc[i][0] != data.iloc[j][0]:
                result.append([data.iloc[i][0],data.iloc[j][0]])
                final_links = pd.DataFrame(data=result)
                final_links.columns = ['Link from', 'Link to']
    return final_links