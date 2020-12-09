import rfpimp 
import pandas as pd
import numpy as np

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import os
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import pickle
import numba

from bokeh.io import show, output_file
from bokeh.io import show, output_notebook, push_notebook
from bokeh.plotting import figure, output_file, show, ColumnDataSource, save
from bokeh.models import LinearColorMapper, CategoricalColorMapper, CDSView, GroupFilter, HoverTool, Label
from bokeh.palettes import Category20
from bokeh.transform import factor_cmap
from bokeh.layouts import gridplot
from bokeh.models.glyphs import Text
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, BoxZoomTool, ResetTool, Legend
from bokeh.palettes import Spectral4
from bokeh.models.graphs import from_networkx, EdgesAndLinkedNodes, NodesAndLinkedEdges

import networkx as nx
import matplotlib.pyplot as plt

from definitions import *

import argparse

def load_data(path, fname):
    return pd.read_csv(os.path.join(path, fname))

def remove_from_list(in_list, to_be_removed):
    return [l for l in in_list if l not in to_be_removed]

def build_rfpimp_graph(rfpimp_output, gtype='directed', meta=['name', 'bet_cent', 'deg_cent', 'weight']):
    rfpimp_df = rfpimp_output.copy()
    
    #Remove dependence score
    if 'Dependence' in rfpimp_df.columns:
        print("Dependence column found, removing it")
        rfpimp_df.drop('Dependence', axis=1, inplace=True)
    rfpimp_df = rfpimp_all.astype(float)
    
    #zero diagonal to cancel self loops
    for i in range(rfpimp_df.shape[0]):
        rfpimp_df.iloc[i,i] = 0            
        
    if gtype == 'undirected':
        graph = nx.from_numpy_matrix(rfpimp_all.values)
    elif gtype == 'directed':
        graph = nx.from_numpy_matrix(rfpimp_all.values, create_using=nx.DiGraph())
    
    
    degrees = [len(list(graph.neighbors(n))) for n in graph.nodes()]
    bet_cen = nx.betweenness_centrality(graph)
    
    if 'name' in meta:
        nx.set_node_attributes(graph, dict(enumerate(rfpimp_all.columns)), 'name')
    if 'deg_cent' in meta:
        deg_cent = nx.degree_centrality(graph) #Returns a dictionary
        nx.set_node_attributes(graph, deg_cent, 'deg_cent')
    if 'bet_cent' in meta:
        bet_cen = nx.betweenness_centrality(graph)
        nx.set_node_attributes(graph, bet_cen, 'bet_cent')
        
    return graph


class MetaboliteConditionCorrs:
    """
    Takes in the correlations of a certain filtering conditions, computes
    the tsne and clustering. Keeps track of:
    - Original correlations data
    - tsne x and y coordinates
    - cluster names
    """
    def __init__(self, name, data_df, cols_to_drop, random_seed, tsne_params, dbscan_params, metabo_special_groups_dict):
        self.name = name
        self.original_data = data_df.copy()
        self.original_data.drop(cols_to_drop, axis=1, inplace=True)
        self.col_names = list(self.original_data.columns)
        
        self.metabo_special_groups_dict = metabo_special_groups_dict
        self.random_seed = random_seed
        
        self.pearson_plot_df = None
        self.spearman_plot_df = None
        self.kindall_plot_df = None
        self.copula_plot_df = None
        self.rfpimp_plot_df = None
        
        self.pearson_df = self.compute_correlation('pearson')
        print('Pearson correlation completed')
        
        self.spearman_df = self.compute_correlation('spearman')
        print('Spearman correlation completed')
        
        self.kindall = self.compute_correlation('kindall')
        print('Kindall correlation completed')
        
        self.copula = self.compute_correlation('copula')
        print('Copula correlation completed')
        
        self.rfpimp_feature_dependence = rfpimp.feature_dependence_matrix(self.original_data)
        print('Features InterDependence completed')
        
        self.plot_df = None
        self.gc = None
        
        self.compute_all_tsne(random_seed, tsne_params)
        print('TSNE completed')
        
        self.compute_clustering(dbscan_params)
        print('Clustering completed')
        
        #Add special groups dict to self.plot_df
        self.add_special_features_to_all()
        
        print(f"Condition {name} processed.")
        
    def compute_all_tsne(self, random_seed, tsne_params):
        self.pearson_plot_df = MetaboliteConditionCorrs.compute_tsne(self.pearson_df, 
                                                             ppl=tsne_params['ppl'],
                                                             n_iters=tsne_params['n_iters'],
                                                             random_state=self.random_seed,
                                                             lr=tsne_params['lr'],
                                                             col_names=self.col_names)
        
        self.spearman_plot_df = MetaboliteConditionCorrs.compute_tsne(self.spearman_df, 
                                                             ppl=tsne_params['ppl'],
                                                             n_iters=tsne_params['n_iters'],
                                                             random_state=self.random_seed,
                                                             lr=tsne_params['lr'],
                                                             col_names=self.col_names)
        
        self.kindall_plot_df = MetaboliteConditionCorrs.compute_tsne(self.kindall, 
                                                             ppl=tsne_params['ppl'],
                                                             n_iters=tsne_params['n_iters'],
                                                             random_state=self.random_seed,
                                                             lr=tsne_params['lr'],
                                                             col_names=self.col_names)
        
        self.copula_plot_df = MetaboliteConditionCorrs.compute_tsne(self.copula, 
                                                             ppl=tsne_params['ppl'],
                                                             n_iters=tsne_params['n_iters'],
                                                             random_state=self.random_seed,
                                                             lr=tsne_params['lr'],
                                                             col_names=self.col_names)
        
        self.rfpimp_plot_df = MetaboliteConditionCorrs.compute_tsne(self.rfpimp_feature_dependence, 
                                                             ppl=tsne_params['ppl'],
                                                             n_iters=tsne_params['n_iters'],
                                                             random_state=self.random_seed,
                                                             lr=tsne_params['lr'],
                                                             col_names=self.col_names)
        
    def compute_clustering(self, dbscan_params):
        
        self.pearson_plot_df['cluster'] = self.compute_dbscan(self.pearson_plot_df,
                                                              eps = dbscan_params['eps'],
                                                              min_samples=dbscan_params['min_samples'])
        
        self.spearman_plot_df['cluster'] = self.compute_dbscan(self.spearman_plot_df,
                                                              eps = dbscan_params['eps'],
                                                              min_samples=dbscan_params['min_samples'])
        
        self.kindall_plot_df['cluster'] = self.compute_dbscan(self.kindall_plot_df,
                                                              eps = dbscan_params['eps'],
                                                              min_samples=dbscan_params['min_samples'])
        
        self.copula_plot_df['cluster'] = self.compute_dbscan(self.copula_plot_df,
                                                              eps = dbscan_params['eps'],
                                                              min_samples=dbscan_params['min_samples'])
        
        self.rfpimp_plot_df['cluster'] = self.compute_dbscan(self.rfpimp_plot_df,
                                                              eps = dbscan_params['eps'],
                                                              min_samples=dbscan_params['min_samples'])
        
    def get_clusters(self, cluster_id, df_name):
        if df_name == 'pearson':
            df_to_use = self.pearson_plot_df
        elif df_name == 'spearman':
            df_to_use = self.spearman_plot_df
        elif df_name == 'kindall':
            df_to_use = self.kindall_plot_df
        elif df_name == 'copula':
            df_to_use = self.copula_plot_df
        elif df_name == 'rfpimp':
            df_to_use = self.rfpimp_plot_df
        else:
            print("df_name not recognized")
            return None
        
        clusters_dict = {}
        for cluster_id in df_to_use['cluster'].unique():
            clusters_dict[cluster_id] = df_to_use.index[df_to_use['cluster'] == cluster_id].tolist()
            
        return clusters_dict
    
    def compute_tsne(df, ppl, n_iters, random_state, lr, col_names):
        X_embedded = TSNE(n_components=2, perplexity=ppl, n_iter=n_iters, random_state=random_state,
                      learning_rate=lr).fit_transform(df)
        
        plot_df = pd.DataFrame()
        plot_df['name'] = col_names
        plot_df['x'] = [a for (a, b) in X_embedded]
        plot_df['y'] = [b for (a, b) in X_embedded]
        
        return plot_df
    
    def compute_dbscan(self, df, eps, min_samples):
        clustering = DBSCAN(eps=eps,
                        min_samples=min_samples).fit(np.array([[x, y] for x, y in zip(df['x'], df['y'])]))
        
        return clustering.labels_.astype(str)
    
    def compute_correlation(self, corr_type):
        if corr_type in ['pearson', 'spearman', 'kindall']:
            return self.original_data.corr('pearson')
        elif corr_type == 'copula':
            self.gc = GaussianMultivariate()
            self.gc.fit(self.original_data.values)
            return pd.DataFrame(self.gc.covariance, 
                                columns=self.original_data.columns, 
                                index=self.original_data.columns)
        
    def add_special_groups(self, plot_df):
        for k,v in self.metabo_special_groups_dict.items():
            plot_df[k] = 0
            plot_df.loc[plot_df['name'].isin(v), k] = 1
            
        return plot_df
    
    def add_special_features_to_all(self):
        self.pearson_plot_df = self.add_special_groups(self.pearson_plot_df)
        self.spearman_plot_df = self.add_special_groups(self.spearman_plot_df)
        self.kindall_plot_df = self.add_special_groups(self.kindall_plot_df)
        self.copula_plot_df = self.add_special_groups(self.copula_plot_df)
        self.rfpimp_plot_df = self.add_special_groups(self.rfpimp_plot_df)

def swap_dict_k_vs(d):
    if d is not None:
        ret_dict = {}
        for k, vs in d.items():
            for v in vs:
                ret_dict[v] = k
        return ret_dict
    return None

@numba.njit
def maxify_diag(in_matrix, result_mat):
    for i in range(in_matrix.shape[0]):
        for j in range(in_matrix.shape[1]):
            #sum
            #result_mat[i,j] =  in_matrix[i,j] + in_matrix[j,i]
            
            #max
            result_mat[i,j] =  max([in_matrix[i,j] + in_matrix[j,i]])
            
            
    return result_mat


def build_networkx(network_map_df, gtype='directed', meta=['name', 'bet_cent', 'deg_cent', 'weight'], edge_cutoff=None):
    df = network_map_df.copy()
    
    #We need all nodes connections to be floats
    df = df.astype(float)
    
    #zero diagonal to cancel self loops
    for i in range(df.shape[0]):
        df.iloc[i,i] = 0      
    
    if edge_cutoff is not None:
        df = df.mask(df < edge_cutoff)
        df = df.fillna(0)
        
    if gtype == 'undirected':
        graph = nx.from_numpy_matrix(df.values)
    elif gtype == 'directed':
        graph = nx.from_numpy_matrix(df.values, create_using=nx.DiGraph())
        
    degrees = [len(list(graph.neighbors(n))) for n in graph.nodes()]
    bet_cen = nx.betweenness_centrality(graph)
    
    if 'name' in meta:
        nx.set_node_attributes(graph, dict(enumerate(rfpimp_all.columns)), 'name')
    if 'deg_cent' in meta:
        deg_cent = nx.degree_centrality(graph) #Returns a dictionary
        nx.set_node_attributes(graph, deg_cent, 'deg_cent')
    if 'bet_cent' in meta:
        bet_cen = nx.betweenness_centrality(graph)
        nx.set_node_attributes(graph, bet_cen, 'bet_cent')
        
    return graph

def build_network_and_bkdatasource(network_map_df, groups=None, color_dict=None, shape_dict=None, gtype='directed', 
                                   meta=['name', 'bet_cent', 'deg_cent', 'weight'],
                                   default_color='gainsboro', default_shape='o', 
                                   edge_cutoff=None, seed=None):
    
    graph = build_networkx(network_map_df, gtype=gtype, meta=meta, edge_cutoff=edge_cutoff)
    nodes_pos = nx.spring_layout(graph, seed=seed)
    
    color_dict = swap_dict_k_vs(color_dict)
    shape_dict = swap_dict_k_vs(shape_dict)
    
    data_source_dicts = []
    
    for i, (n, d) in enumerate(graph.nodes(data=True)):
        #set default color and shape, and if we find this entry in the groups, we will update the default value
        node_col = default_color
        node_shape = default_shape
        
        if groups is not None:
            for group in groups:
                if (d['name'] in group) and (color_dict is not None) and (d['name'] in color_dict.keys()):
                    node_col = color_dict[d['name']] 

                if (d['name'] in group) and (shape_dict is not None) and (d['name'] in shape_dict.keys()):
                    node_shape = shape_dict[d['name']]
                    
        sd = {}
        pos = nodes_pos[i]

        sd['x'] = pos[0]
        sd['y'] = pos[1]
        sd['name'] = d['name']
        sd['color'] = node_col
        sd['shape'] = node_shape
        data_source_dicts.append(sd)
            
    return graph , pd.DataFrame(data_source_dicts), nodes_pos


def plot_netwok_bk(source_df, g, nodes_pos):
    source = ColumnDataSource(data=source_df)
    circles = CDSView(source=source, filters=[GroupFilter(column_name=f'shape', group='o')])
    triangles = CDSView(source=source, filters=[GroupFilter(column_name=f'shape', group='^')])
    squares = CDSView(source=source, filters=[GroupFilter(column_name=f'shape', group='s')])
    diamonds = CDSView(source=source, filters=[GroupFilter(column_name=f'shape', group='d')])
    stars = CDSView(source=source, filters=[GroupFilter(column_name=f'shape', group='*')])
    
    TOOLTIPS = [("metabo_name", "@name")]
    compare_all_TOOLS = "box_select,lasso_select,help,save, wheel_zoom, box_zoom, pan, reset"
    
    fig=figure(title='Metabolite RF Permutation Importance Network',
           plot_width=1200, 
           plot_height=600, #tooltips=TOOLTIPS, 
           tools=compare_all_TOOLS, toolbar_location="above")
    
    dummy_spearator0 = fig.circle(0, 0, size=15, alpha=1.0, color='white')
    dummy_c1 = fig.circle(0, 0, size=15, alpha=1.0, line_color='red', color='white', line_width=2)
    dummy_c2 = fig.circle(0, 0, size=15, alpha=1.0, line_color='green', color='white', line_width=2)
    dummy_c3 = fig.circle(0, 0, size=15, alpha=1.0, line_color='gainsboro', color='white', line_width=2)
    dummy_hider0 = fig.circle(0, 0, size=15, alpha=1.0, color='white')
    dummy_hider1 = fig.circle(0, 0, size=15, alpha=1.0, color='white')

    r = from_networkx(g, nodes_pos, scale=1, center=(0,0))
    r.node_renderer.glyph = Circle(size=4, fill_color='#000000')
    r.node_renderer.hover_glyph = Circle(size=15, fill_color='#abdda4')
    node_hover_tool = HoverTool(tooltips=[("name", "@name")])
    fig.add_tools(node_hover_tool)

    r.edge_renderer.glyph = MultiLine(line_alpha=0.1, line_width=1)  # zero line alpha
    r.edge_renderer.hover_glyph = MultiLine(line_color='#abdda4', line_width=2)

    r.inspection_policy = NodesAndLinkedEdges()
    fig.renderers.append(r)

    
    
    dummy_spearator = fig.circle(0, 0, size=1, alpha=1.0, color='white')
    
    ts = fig.triangle('x', 'y', size=20, source=source, view=triangles,
                 color='color', line_width=2, line_color='black',
                 nonselection_alpha=0.01)

    sq = fig.square('x', 'y', size=12, source=source, view=squares,
                 color='color', line_width=2, line_color='black',
               nonselection_alpha=0.01)

    ds = fig.diamond('x', 'y', size=15, source=source, view=diamonds,
                 color='color', line_width=2, line_color='black',
                nonselection_alpha=0.01)
    
    c = fig.circle('x', 'y', size=8, source=source, view=circles, alpha=0.6,
                   color='color',
                   nonselection_alpha=0.01)
    
    st = fig.circle_cross('x', 'y', size=12, source=source, view=stars,
                 color='color', line_width=2, line_color='black',
                nonselection_alpha=0.01)
    
    

    fig.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    fig.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks

    fig.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    fig.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks

    fig.xaxis.major_label_text_color = None  # turn off x-axis tick labels leaving space
    fig.yaxis.major_label_text_color = None  # turn off y-axis tick labels leaving space 
    
    
    
    legend1 = Legend(items=[("Colors:", [dummy_spearator0]),
                            ("Omega 3 or 6" , [dummy_c1]), 
                            ("T-SNE Green Cluster", [dummy_c2]),
                            #('Everything else', [dummy_c3]),
                            ("   ", [dummy_hider0]),
                            ("Shapes:", [dummy_hider1]),
                            #("Non-significant" , [c]), 
                            #("From the ANOVA model:", [dummy_spearator]),
                            ("ageXstate interaction metabolites", [ts]), 
                            ("treatmentXstate interaction metabolites", [sq]), 
                            ("timeXstate interaction metabolites", [ds]),
                            ("sexXstate interaction metabolites", [st])] ) 
                    #location=(10,1))#, orientation="horizontal")
    

    fig.add_layout(legend1, 'right')
    
    return fig


def parse_args():
    """
    Parse the arguments from the command line
    All arguments have default values
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--data_path',  # Specified in the config file
                        help='Folder containing all the datasets',
                        default='../data',
                        type=str)
    
    parser.add_argument('--f1name',  # Specified in the config file
                        help='Name of the normalized dataset file. Disregarded if we reload precomuted file.',
                        default='data_normalized_all_BL_DC.csv',
                        type=str)
    
    parser.add_argument('--f2name',  # Specified in the config file
                        help='Name of the metadata file. Disregarded if we reload precomuted file.',
                        default='patient_IDs_metadata.csv',
                        type=str)
    
    parser.add_argument('--load_precomputed',  # Specified in the config file
                        help='If the file has already been computed, just reuse it',
                        default=True,
                        type=bool)
    
    parser.add_argument('--precomputed_file',  # Specified in the config file
                        help='Name of the metadata file. Note that this argument is still used even if we recompute from scratch, in which case this is the name of the saved file.',
                        default='condition_all.pickle',
                        type=str)
    
    parser.add_argument('--edge_thresh',  # Specified in the config file
                        help='Threshold of predictibility below which we filter out the edge',
                        default=0.015,
                        type=float)
    
    parser.add_argument('--random_seed',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=123,
                        type=int)
    
    parser.add_argument('--tsne_ppl',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=7,
                        type=int)
    
    parser.add_argument('--tsne_lr',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=200,
                        type=int)
    
    parser.add_argument('--dbscan_eps',  # Specified in the config file
                        help='random seed. It can help with reproducibility',
                        default=8,
                        type=int)
    
    parser.add_argument('--dbscan_min_samples',  # Specified in the config file
                        help='dbscan_ seed. It can help with reproducibility',
                        default=4,
                        type=int)
    return parser.parse_args()

    
if __name__ == "__main__":
    parser = parse_args()
    
    edge_thresh = parser.edge_thresh #0.015
    tsne_params = {'ppl':parser.tsne_ppl, 'n_iters':5000, 'lr':parser.tsne_lr}
    dbscan_params = {'eps':parser.dbscan_eps, 'min_samples':parser.dbscan_eps}
    metabo_special_groups_dict = shapes_dict
    cols_to_drop=['samples', 'SAMPLE_NAME', 'AGE', 'SEX', 'CONDITION']
    
    data_path = parser.data_path  #'../Clean/Data/ml_preprocessed'
    np.random.seed(parser.random_seed)
    
    data_path = parser.data_path #'../Data/original'
    file_1 = parser.f1name#"data_normalized_all_BL_DC.csv"
    file_2 = parser.f2name#"patient_IDs_metadata.csv"
    
    #First time run, set it to false, it will take a bit of time
    load_precomputed = parser.load_precomputed
    
    if load_precomputed:
        file = open(os.path.join(parser.data_path, parser.precomputed_file),'rb')
        condition_all = pickle.load(file)
        file.close()
        
    else:
        first_analysis_df = pd.read_csv(os.path.join(data_path,file_1), sep='\t')
        meta_1st_analysis_df = pd.read_csv(os.path.join(data_path,file_2), sep='\t')

        #Rename the 1st column to be samples
        first_analysis_df['samples'] = first_analysis_df['metabolites']
        meta_1st_analysis_df['CONDITION'] = meta_1st_analysis_df['TIME']
        
        filtered_df = first_analysis_df.copy()
        filtered_df = pd.merge(filtered_df, meta_1st_analysis_df[['SAMPLE_NAME','AGE', 'SEX', 'CONDITION', 'TREATMENT']], left_on='samples', right_on='SAMPLE_NAME')
        filtered_df.drop(['metabolites'], axis=1, inplace=True)
        
        condition_all = MetaboliteConditionCorrs(name='all', 
                                            data_df=filtered_df, 
                                            cols_to_drop=cols_to_drop,
                                            random_seed=parser.random_seed,
                                            tsne_params=tsne_params,
                                            dbscan_params=dbscan_params,
                                            metabo_special_groups_dict=metabo_special_groups_dict)
    
        with open(os.path.join(parser.data_path, parser.precomputed_file), 'wb') as f:
            pickle.dump(condition_all, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    rfpimp_all = condition_all.rfpimp_feature_dependence.copy()
    #Drop the extra score row
    rfpimp_all.drop('Dependence', axis=1, inplace=True)
    rfpimp_all = rfpimp_all.astype(float)
    
    #zero the diagonals, otherwise it will be interpreted as self loops
    for i in range(rfpimp_all.shape[0]):
        rfpimp_all.iloc[i,i] = 0
    
    result_df = rfpimp_all.copy()
    result_matrix = np.zeros_like(rfpimp_all.values, dtype=np.float32)
    rfpimp_vals = rfpimp_all.values.astype(np.float32)
    
    #Since the matrix is not diagonal but we are not interested in the directed form of the graph, we just
    #take the maximum of both directions and copy it to the diagonal
    result_matrix = maxify_diag(rfpimp_vals, result_matrix)
    result_df.iloc[:,:] = result_matrix
    
    #Build the graph
    colors_dict = {'green':green_mets, 'red':omega_3_mets}
    groups = [green_mets, metabo_agePstatus, metabo_treat, metabo_timePstatus, metabo_sex]
    
    gg, dff, nodes_poss = build_network_and_bkdatasource(result_df, groups, colors_dict, shapes_dict, 
                                                     edge_cutoff=0.015, gtype='undirected', seed=456)
    fig = plot_netwok_bk(dff, gg, nodes_poss)

    output_file("rfpimp_network.html")
    save(fig)
