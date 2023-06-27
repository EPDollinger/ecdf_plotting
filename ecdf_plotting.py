import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Union
import scanpy as sc

def plot_ecdf(
        adata: AnnData, 
        group: str, 
        genes: Union[str, list],
        hue_group: Union[str, type(None)] = None,
        col_wrap: Union[int, type(None)] = None,
        height: float = 5.0, 
        plot_on_x: bool = False,
        return_df: bool = False,  
        **kwargs):
    '''
    Takes adata, a grouping, genes, and returns a df and displot of the gene expression in the groups.
    Parameters
    ----------
    adata - AnnData object
    group - str, obs column name in AnnData object to do columns by.
    genes - [str, list], genes expression to plot.
    hue_group - [str, None], obs column name. If plotting only one gene allow hue to be a categorical variable.
    col_wrap - [int, None], how many columns to plot. Defaults to None.
    height - float, height of plot. Defaults to 5.0.
    plot_on_x - bool, whether to plot expression on x axis. Defaults to False (plot expression on y axis).
    return_df - bool, whether to return melted dataframe for further plotting. Defaults to False.
    '''
    if isinstance(genes, str):
        genes = [genes]
    if hue_group is not None:
        data_to_extract = np.hstack((genes, group, hue_group))
        id_vars = [group, hue_group]
        kwargs_for_plotting = {'col':group, 'hue':hue_group, 'row':'Gene'}
    
    else:
        data_to_extract = np.hstack((genes, group))
        id_vars = [group]
        kwargs_for_plotting = {'col':'Gene', 'hue':group}
    df = sc.get.obs_df(adata, data_to_extract.tolist())
    melteddf = df.melt(id_vars=id_vars, var_name='Gene', value_name='Expression')
    if plot_on_x:
        g = sns.displot(data=melteddf, x='Expression', kind='ecdf', col_wrap=col_wrap, height=height, **kwargs_for_plotting, **kwargs)
        plt.yticks([0,0.25,0.5,0.75,1.0])
    else:
        g = sns.displot(data=melteddf, y='Expression',  kind='ecdf', col_wrap=col_wrap, height=height, **kwargs_for_plotting, **kwargs)
        plt.xticks([0,0.25,0.5,0.75,1.0])
    
    
    if return_df:
        return g, melteddf
    return g