## Dynamic Slicing and Attractor Analysis of Biological Networks

import networkx as nx
import pandas as pd
from functools import reduce
from swiplserver import PrologMQI, PrologThread
import sys
import os
import glob
from itertools import chain, combinations
import json

class Target:
    def __init__(self,pr,ab):
        self.present = pr
        self.absent = ab
    
    def __str__(self):
        return f"present: {self.present} absent: {self.absent}"    
        
targets = []
interesting_entities = ["tyr_b", "bh4_b"]
def powerset(s):
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
for x in powerset(interesting_entities):
    if (len(x)>0):
        targets.append(Target([k for k in x],[k for k in interesting_entities if k not in x]))

        # checks whether a node is in an attractor
def check_node(node):
    cycle = list(nx.find_cycle(G,node))
    tmp = map(lambda x : x[0]==node or x[1]==node, cycle) 
    return reduce(lambda b1, b2: b1 or b2, tmp)

# Compute the list of entities that are present in the attractor reachable from "node"
def compute_attractor(node):
    cycle = list(nx.find_cycle(G,node))
    tmp1 = map(lambda x: x[0].split(';') + x[1].split(';'), cycle)
    tmp2 = reduce(lambda x, y: x+y,tmp1)
    res = list(dict.fromkeys(tmp2))
    return res

def count_states (d):
    tmp = 0
    for v in d.values():
        tmp += len(v)
    return tmp

# Find computations (LTS traces) that lead to the selected target
def target_computations(target):
    all_nodes = list(G.nodes)
    filtered = [k for k in all_nodes if check_node(k)] # Filter out intermediate nodes that are not in an attractor
    attractors_map = dict()
    for f in filtered:
        attractors_map[f] = compute_attractor(f) # Create a "state -> attractor" map
    for s in target.present:
        filtered = [k for k in filtered if s in attractors_map[k]]   # Filter out states that do not contain s, which is required to be present in the target
    for s in target.absent:
        filtered = [k for k in filtered if s not in attractors_map[k]] # Filter out states that contain s, which is required to be absent from the target

    # Filter out states in target attractors that do not contain any entity from target.present
    # (slicing would give an empty result on these states)
    filtered2 = [k for k in filtered if len(list(set(k.split(';')) & set(target.present)))>0]

    # Remove from attractors_map the states that are not part of the target set
    keys_to_delete = list() 
    for k in attractors_map.keys():
        if k not in filtered: keys_to_delete.append(k) 
    for k in keys_to_delete:
        del attractors_map[k] 
    
    contexts = list()   # List of contexts that lead to the target
    contexts_dict = {}  # For each context, list the states in the corresponding attractor
    filtered_splitted = [f.split(';') for f in filtered2]
    for f in filtered_splitted:
        prefix = f[0:10] # The first 10 elements in the state represent the context [UPDATE 10 if the context size changes]
        pure_state = f[10:] # The remaining elements represent the actual state
        if (not prefix in contexts):
            contexts.append(prefix)
            contexts_dict[','.join(prefix)] = [','.join(pure_state)]
        else:
            contexts_dict[','.join(prefix)].append(','.join(pure_state))
    print("TARGET --> " + str(target))
    print("Found " + str(len(contexts)) + " contexts that lead to the target")
    print("Found " + str(count_states(contexts_dict) + (len(filtered)-len(filtered2))) + " states in reachable attractors")
    print("Of these, " + str(count_states(contexts_dict)) + " contain entities from the target")
    print()
    return contexts, contexts_dict


from IPython.display import display

# Create a global DataFrame to accumulate rows
table_data_df = pd.DataFrame()

def generate_table(target_desc, intersection_pos_set, union_pos_set, intersection_neg_set, union_neg_set, show_neg=True):
    global table_data_df
        
    # Get all unique entity names
    all_strings = sorted(set(table_data_df.columns.get_level_values(0)).union(union_pos_set, union_neg_set))
    
    # Create the column header by alternating the POS and NEG subcolumns
    columns = [(s, cat) for s in all_strings for cat in ('POS', 'NEG')]
    
    # Ensure that the DataFrame contains all required columns
    if table_data_df.empty:
        table_data_df = pd.DataFrame(columns=pd.MultiIndex.from_tuples(columns))
    else:
        missing_columns = set(columns) - set(table_data_df.columns)
        for col in missing_columns:
            table_data_df[col] = 0
        table_data_df = table_data_df[columns]  # Reorder the columns
    
    # Create a new row
    new_row = {col: 0 for col in table_data_df.columns}  # Initialize all columns to 0
    for s in all_strings:
        if s in intersection_pos_set:
            new_row[(s, 'POS')] = 2
        elif s in union_pos_set:
            new_row[(s, 'POS')] = 1
        if s in intersection_neg_set:
            new_row[(s, 'NEG')] = 2
        elif s in union_neg_set:
            new_row[(s, 'NEG')] = 1
    
    # Add the new row to the DataFrame and fill any NaN values with 0
    table_data_df = pd.concat([table_data_df, pd.DataFrame([new_row], index=[target_desc])]).fillna(0)
    
    # Function used to apply the style
    def highlight_cells(val):
        if val == 2:
            return 'background-color: black; color: white'
        elif val == 1:
            return 'background-color: gray; color: white'
        return ''
  
    # Apply the style
    styled_df = table_data_df.style.applymap(highlight_cells)
    display(styled_df)


def save_table_to_file(filename, show_neg=True):
    """
    Save the DataFrame to a CSV file.
    """
    global table_data_df
    if show_neg:
        df_to_save = table_data_df
    else:
        # Save only the POS columns with a single-level header
        df_to_save = table_data_df.xs('POS', level=1, axis=1)
        df_to_save.columns = df_to_save.columns.get_level_values(0)
    # if df.empty:
    #     print("No data to save.")
    # else:
        
        df_to_save.to_csv(filename, index=True)


def slicing_analysis(version):
    target_tot = len(targets)
    target_count = 0
    # enzymes = ["tat","dhpr","pah_fe2","gstz1","hpdfe","hgdfe","fahmgca"]
    # pat_enzymes_count = {enzyme: [0, 0] for enzyme in enzymes}
    sliced_computations = []
    print("TARGETS TO ANALYZE: " + str(target_tot))
    print()

    #outfile = open("Results/"+dot_name+"/out.txt","w")
    if not os.path.exists("Results/Final/slicing/"+dot_name):
        os.makedirs("Results/Final/slicing/"+dot_name)

    #outfile = open("Results/slicing/"+dot_name+"/out_"+str(targets[0].present)+".txt","w")
    for target in targets:

        outfile = open("Results/Final/slicing/"+dot_name+"/out_"+str(target.present)+".txt","w")
        target_count = target_count+1
        print("TARGET COUNT: " + str(target_count) + "/" + str(target_tot))
        contexts, contexts_dict = target_computations(target)

        prolog_target = ','.join(target.present)

        cont = 1
        tot = str(count_states(contexts_dict))
        union_set = set()
        intersection_set = set()
        first_time = True
        for ctx in contexts:
            prolog_context = ','.join(ctx)
            prolog_target_states = contexts_dict[','.join(ctx)]
            for i,state in enumerate(prolog_target_states):
                print("TEST CASE: " + str(cont) + "/" + str(tot))
                cont=cont+1
                print("CONTEXT: " + prolog_context + "      ATTRACTOR STATE: " + str(i+1) + "/" + str(len(prolog_target_states)))
                print("STATE: " + state)
            
                param_file = open("params.pl",'w') 
                if (version=="original"):
                    #param_file.write("myenvironment('[nul =  {}.nul,\n bh2_p_b = {bh2_b}.bh2_p_b,\n bh2_p_n = {bh2}.bh2_p_n,\n bh4_p_b = {bh4_b}.bh4_p_b,\n bh4_p_n = {bh4}.bh4_p_n,\n phe_p_b = {phe_b}.phe_p_b,\n phe_p_n = {phe}.phe_p_n,\n tat_p_b = {tat_b}.tat_p_b,\n tat_p_n = {tat}.tat_p_n,\n plp_p_b = {plp_b}.plp_p_b,\n plp_p_n = {plp}.plp_p_n,\n alphakg_p_b = {alphakg_b}.alphakg_p_b,\n alphakg_p_n = {alphakg}.alphakg_p_n,\n dhpr_p_b = {dhpr_b}.dhpr_p_b,\n dhpr_p_n = {dhpr}.dhpr_p_n,\n pah_p_b = {pah_fe2_b}.pah_p_b,\n pah_p_n = {pah_fe2}.pah_p_n,\n gstz1_p_b = {gstz1_b}.gstz1_p_b,\n gstz1_p_n = {gstz1}.gstz1_p_n,\n hpd_p_b = {hpdfe_b}.hpd_p_b,\n hpd_p_n = {hpdfe}.hpd_p_n,\n hgdfe_p_b = {hgdfe_b}.hgdfe_p_b, \n fah_p_b = {fahmgca_b}.fah_p_b,\n fah_p_n = {fahmgca}.fah_p_n \n ]').\n \n")
                    #Try differentiating entities from the context at the second iteration
                    #param_file.write("myenvironment('[nul =  {}.nul,\n bh2_p_b = {bh2_b}.bh2_p_b1,\n bh2_p_b1 = {xbh2_b}.bh2_p_b1,\n bh2_p_n = {bh2}.bh2_p_n1,\n  bh2_p_n1 = {xbh2}.bh2_p_n1,\n  bh4_p_b = {bh4_b}.bh4_p_b1,\n  bh4_p_b1 = {xbh4_b}.bh4_p_b1,\n  bh4_p_n = {bh4}.bh4_p_n1,\n  bh4_p_n1 = {xbh4}.bh4_p_n1,\n  phe_p_b = {phe_b}.phe_p_b1,\n  phe_p_b1 = {xphe_b}.phe_p_b1,\n  phe_p_n = {phe}.phe_p_n1,\n  phe_p_n1 = {xphe}.phe_p_n1,\n  tat_p_b = {tat_b}.tat_p_b1,\n  tat_p_b1 = {xtat_b}.tat_p_b1,\n  tat_p_n = {tat}.tat_p_n1,\n  tat_p_n1 = {xtat}.tat_p_n1,\n  plp_p_b = {plp_b}.plp_p_b1,\n  plp_p_b1 = {xplp_b}.plp_p_b1,\n  plp_p_n = {plp}.plp_p_n1,\n  plp_p_n1 = {xplp}.plp_p_n1,\n  alphakg_p_b = {alphakg_b}.alphakg_p_b1,\n  alphakg_p_b1 = {xalphakg_b}.alphakg_p_b1,\n  alphakg_p_n = {alphakg}.alphakg_p_n1,\n  alphakg_p_n1 = {xalphakg}.alphakg_p_n1,\n  dhpr_p_b = {dhpr_b}.dhpr_p_b1,\n  dhpr_p_b1 = {xdhpr_b}.dhpr_p_b1,\n  dhpr_p_n = {dhpr}.dhpr_p_n1,\n  dhpr_p_n1 = {xdhpr}.dhpr_p_n1,\n  pah_p_b = {pah_fe2_b}.pah_p_b1,\n  pah_p_b1 = {xpah_fe2_b}.pah_p_b1,\n  pah_p_n = {pah_fe2}.pah_p_n1,\n pah_p_n1 = {xpah_fe2}.pah_p_n1,\n  gstz1_p_b = {gstz1_b}.gstz1_p_b1,\n  gstz1_p_b1 = {xgstz1_b}.gstz1_p_b1,\n  gstz1_p_n = {gstz1}.gstz1_p_n1,\n  gstz1_p_n1 = {xgstz1}.gstz1_p_n1,\n  hpd_p_b = {hpdfe_b}.hpd_p_b1,\n  hpd_p_b1 = {xhpdfe_b}.hpd_p_b1,\n hpd_p_n = {hpdfe}.hpd_p_n1,\n  hpd_p_n1 = {xhpdfe}.hpd_p_n1,\n  hgdfe_p_b = {hgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_b1 = {xhgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_n = {hgdfe}.hgdfe_p_n1,\n  hgdfe_p_n1 = {xhgdfe}.hgdfe_p_n1,\n    fah_p_b = {fahmgca_b}.fah_p_b1,\n  fah_p_b1 = {xfahmgca_b}.fah_p_b1,\n  fah_p_n = {fahmgca}.fah_p_n1,\n fah_p_n1 = {xfahmgca}.fah_p_n1 ,\n  nadh_p_n = {nadh}.nadh_p_n1,\n  nadh_p_n1 = {xnadh}.nadh_p_n1,\n gsh_p_n = {gsh}.gsh_p_n1,\n  gsh_p_n1 = {xgsh}.gsh_p_n1,\n ala_p_n = {ala}.ala_p_n1,\n  ala_p_n1 = {xala}.ala_p_n1,\n arat_p_n = {arat}.arat_p_n1,\n  arat_p_n1 = {xarat}.arat_p_n1   \n]').\n \n") 
                    #param_file.write("myenvironment('[nul =  {}.nul,\n bh2_p_b = {bh2_b}.bh2_p_b1,\n bh2_p_b1 = {xbh2_b}.bh2_p_b1,\n bh2_p_n = {bh2}.bh2_p_n1,\n  bh2_p_n1 = {xbh2}.bh2_p_n1,\n  bh4_p_b = {bh4_b}.bh4_p_b1,\n  bh4_p_b1 = {xbh4_b}.bh4_p_b1,\n  bh4_p_n = {bh4}.bh4_p_n1,\n  bh4_p_n1 = {xbh4}.bh4_p_n1,\n  phe_p_b = {phe_b}.phe_p_b1,\n  phe_p_b1 = {xphe_b}.phe_p_b1,\n  phe_p_n = {phe}.phe_p_n1,\n  phe_p_n1 = {xphe}.phe_p_n1,\n tat_p_nul = {tat_nul}.tat_p_nul1,\n tat_p_nul1 = {xtat_nul}.tat_p_nul1,\n tat_p_b = {tat_b}.tat_p_b1,\n  tat_p_b1 = {xtat_b}.tat_p_b1,\n  tat_p_n = {tat}.tat_p_n1,\n  tat_p_n1 = {xtat}.tat_p_n1,\n  plp_p_b = {plp_b}.plp_p_b1,\n  plp_p_b1 = {xplp_b}.plp_p_b1,\n  plp_p_n = {plp}.plp_p_n1,\n  plp_p_n1 = {xplp}.plp_p_n1,\n  alphakg_p_b = {alphakg_b}.alphakg_p_b1,\n  alphakg_p_b1 = {xalphakg_b}.alphakg_p_b1,\n  alphakg_p_n = {alphakg}.alphakg_p_n1,\n  alphakg_p_n1 = {xalphakg}.alphakg_p_n1,\n dhpr_p_nul = {dhpr_nul}.dhpr_p_nul1,\n dhpr_p_nul1 = {xdhpr_nul}.dhpr_p_nul1,\n dhpr_p_b = {dhpr_b}.dhpr_p_b1,\n  dhpr_p_b1 = {xdhpr_b}.dhpr_p_b1,\n  dhpr_p_n = {dhpr}.dhpr_p_n1,\n  dhpr_p_n1 = {xdhpr}.dhpr_p_n1,\n pah_fe2_p_nul = {pah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_nul1 = {xpah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_b = {pah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_b1 = {xpah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_n = {pah_fe2}.pah_fe2_p_n1,\n pah_fe2_p_n1 = {xpah_fe2}.pah_fe2_p_n1,\n gstz1_p_nul = {gstz1_nul}.gstz1_p_nul1,\n gstz1_p_nul1 = {xgstz1_nul}.gstz1_p_nul1,\n gstz1_p_b = {gstz1_b}.gstz1_p_b1,\n  gstz1_p_b1 = {xgstz1_b}.gstz1_p_b1,\n  gstz1_p_n = {gstz1}.gstz1_p_n1,\n  gstz1_p_n1 = {xgstz1}.gstz1_p_n1,\n hpdfe_p_nul = {hpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_nul1= {xhpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_b = {hpdfe_b}.hpdfe_p_b1,\n  hpdfe_p_b1 = {xhpdfe_b}.hpdfe_p_b1,\n hpdfe_p_n = {hpdfe}.hpdfe_p_n1,\n  hpdfe_p_n1 = {xhpdfe}.hpdfe_p_n1, \n hgdfe_p_nul = {hgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_nul1= {xhgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_b = {hgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_b1 = {xhgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_n = {hgdfe}.hgdfe_p_n1,\n  hgdfe_p_n1 = {xhgdfe}.hgdfe_p_n1,\n fahmgca_p_nul = {fahmgca_nul}.fahmgca_p_nul1,\n fahmgca_p_nul1= {xfahmgca_nul}.fahmgca_p_nul1, \n fahmgca_p_b = {fahmgca_b}.fahmgca_p_b1,\n fahmgca_p_b1 = {xfahmgca_b}.fahmgca_p_b1,\n fahmgca_p_n = {fahmgca}.fahmgca_p_n1,\n fahmgca_p_n1 = {xfahmgca}.fahmgca_p_n1 ,\n  nadh_p_n = {nadh}.nadh_p_n1,\n  nadh_p_n1 = {xnadh}.nadh_p_n1,\n gsh_p_n = {gsh}.gsh_p_n1,\n  gsh_p_n1 = {xgsh}.gsh_p_n1,\n ala_p_n = {ala}.ala_p_n1,\n  ala_p_n1 = {xala}.ala_p_n1,\n arat_p_n = {arat}.arat_p_n1,\n  arat_p_n1 = {xarat}.arat_p_n1   \n]').\n \n")  
                    #,\n  nadh_p_n = {nadh}.nadh_p_n1,\n  nadh_p_n1 = {xnadh}.nadh_p_n1,\n gsh_p_n = {gsh}.gsh_p_n1,\n  gsh_p_n1 = {xgsh}.gsh_p_n1 
                    #Try differentiating entities from the context at the first iteration
                    #param_file.write("myenvironment('[nul =  {}.nul,\n bh2_p_b = {xbh2_b}.bh2_p_b,\n bh2_p_n = {xbh2}.bh2_p_n,\n bh4_p_b = {xbh4_b}.bh4_p_b,\n bh4_p_n = {xbh4}.bh4_p_n,\n phe_p_b = {xphe_b}.phe_p_b,\n phe_p_n = {xphe}.phe_p_n,\n tat_p_b = {xtat_b}.tat_p_b,\n tat_p_n = {xtat}.tat_p_n,\n plp_p_b = {xplp_b}.plp_p_b,\n plp_p_n = {xplp}.plp_p_n,\n alphakg_p_b = {xalphakg_b}.alphakg_p_b,\n alphakg_p_n = {xalphakg}.alphakg_p_n,\n dhpr_p_b = {xdhpr_b}.dhpr_p_b,\n dhpr_p_n = {xdhpr}.dhpr_p_n,\n pah_p_b = {xpah_fe2_b}.pah_p_b,\n pah_p_n = {xpah_fe2}.pah_p_n,\n gstz1_p_b = {xgstz1_b}.gstz1_p_b,\n gstz1_p_n = {xgstz1}.gstz1_p_n,\n hpd_p_b = {xhpdfe_b}.hpd_p_b,\n hpd_p_n = {xhpdfe}.hpd_p_n,\n hgdfe_p_b = {xhgdfe_b}.hgdfe_p_b, \n fah_p_b = {xfahmgca_b}.fah_p_b,\n fah_p_n = {xfahmgca}.fah_p_n \n ]').\n \n")
                    param_file.write("myenvironment('[bh4_p_b = {bh4_b}.bh4_p_b1,\n  bh4_p_b1 = {xbh4_b}.bh4_p_b1,\n  bh4_p_n = {bh4}.bh4_p_n1,\n  bh4_p_n1 = {xbh4}.bh4_p_n1,\n  phe_p_b = {phe_b}.phe_p_b1,\n  phe_p_b1 = {xphe_b}.phe_p_b1,\n  phe_p_n = {phe}.phe_p_n1,\n  phe_p_n1 = {xphe}.phe_p_n1,\n tat_p_nul = {tat_nul}.tat_p_nul1,\n tat_p_nul1 = {xtat_nul}.tat_p_nul1,\n tat_p_b = {tat_b}.tat_p_b1,\n  tat_p_b1 = {xtat_b}.tat_p_b1,\n  tat_p_n = {tat}.tat_p_n1,\n  tat_p_n1 = {xtat}.tat_p_n1,\n  plp_p_b = {plp_b}.plp_p_b1,\n  plp_p_b1 = {xplp_b}.plp_p_b1,\n  plp_p_n = {plp}.plp_p_n1,\n  plp_p_n1 = {xplp}.plp_p_n1,\n dhpr_p_nul = {dhpr_nul}.dhpr_p_nul1,\n dhpr_p_nul1 = {xdhpr_nul}.dhpr_p_nul1,\n dhpr_p_b = {dhpr_b}.dhpr_p_b1,\n  dhpr_p_b1 = {xdhpr_b}.dhpr_p_b1,\n  dhpr_p_n = {dhpr}.dhpr_p_n1,\n  dhpr_p_n1 = {xdhpr}.dhpr_p_n1,\n pah_fe2_p_nul = {pah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_nul1 = {xpah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_b = {pah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_b1 = {xpah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_n = {pah_fe2}.pah_fe2_p_n1,\n pah_fe2_p_n1 = {xpah_fe2}.pah_fe2_p_n1,\n gstz1_p_nul = {gstz1_nul}.gstz1_p_nul1,\n gstz1_p_nul1 = {xgstz1_nul}.gstz1_p_nul1,\n gstz1_p_b = {gstz1_b}.gstz1_p_b1,\n  gstz1_p_b1 = {xgstz1_b}.gstz1_p_b1,\n  gstz1_p_n = {gstz1}.gstz1_p_n1,\n  gstz1_p_n1 = {xgstz1}.gstz1_p_n1,\n hpdfe_p_nul = {hpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_nul1= {xhpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_b = {hpdfe_b}.hpdfe_p_b1,\n  hpdfe_p_b1 = {xhpdfe_b}.hpdfe_p_b1,\n hpdfe_p_n = {hpdfe}.hpdfe_p_n1,\n  hpdfe_p_n1 = {xhpdfe}.hpdfe_p_n1,\n hgdfe_p_nul = {hgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_nul1= {xhgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_b = {hgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_b1 = {xhgdfe_b}.hgdfe_p_b1,\n hgdfe_p_n = {hgdfe}.hgdfe_p_n1,\n  hgdfe_p_n1 = {xhgdfe}.hgdfe_p_n1,\n fahmgca_p_nul = {fahmgca_nul}.fahmgca_p_nul1,\n fahmgca_p_nul1= {xfahmgca_nul}.fahmgca_p_nul1,\n fahmgca_p_b = {fahmgca_b}.fahmgca_p_b1,\n fahmgca_p_b1 = {xfahmgca_b}.fahmgca_p_b1,\n fahmgca_p_n = {fahmgca}.fahmgca_p_n1,\n fahmgca_p_n1 = {xfahmgca}.fahmgca_p_n1 \n]').\n \n")  

                    param_file.write("myentities([h,o2,h2o,bh2, alphakg, nadh, gsh, ala, arat]).\n \n")

                    #"mycontext(\"[bh2_p,bh4_p,phe_p,tat_p,plp_p,alphakg_p,dhpr_p,pah_fe2_p,gstz1_p,hpdfe_p,hgdfe_p,fahmgca_p]\").\n")
                    #param_file.write("myenvironment('[ x1 = ({tgfb}.x11 + {}.x0),\n x2 = ({il23}.x21 + {}.x0),\n x3 = ({il12}.x31 + {}.x0),\n x4 = ({il18}.x41 + {}.x0),\n x5 = ({il4e}.x51 + {}.x0),\n x6 = ({il27}.x61 + {}.x0),\n x7 = ({il6e}.x71 + {}.x0),\n x8 = ({ifnge}.x81 + {}.x0),\n x9 = ({tcr}.x91 + {}.x0),\n x11 = {tgfb}.x11,\n x21 = {il23}.x21,\n x31 = {il12}.x31,\n x41 = {il18}.x41,\n x51 = {il4e}.x51,\n x61 = {il27}.x61,\n x71 = {il6e}.x71,\n x81 = {ifnge}.x81,\n x91 = {tcr}.x91,\n x0 = {}.x0\n]').\n \n myentities([]).\n \n")
                elif (version=="positive"):
                    param_file.write("myenvironment('[ x1 = ({tgfb}.x11 + {neg_tgfb}.x12),\n x2 = ({il23}.x21 + {neg_il23}.x22),\n x3 = ({il12}.x31 + {neg_il12}.x32),\n x4 = ({il18}.x41 + {neg_il18}.x42),\n x5 = ({il4e}.x51 + {neg_il4e}.x52),\n x6 = ({il27}.x61 + {neg_il27}.x62),\n x7 = ({il6e}.x71 + {neg_il6e}.x72),\n x8 = ({ifnge}.x81 + {neg_ifnge}.x82),\n x9 = ({tcr}.x91 + {neg_tcr}.x92),\n x11 = {tgfb}.x11,\n x21 = {il23}.x21,\n x31 = {il12}.x31,\n x41 = {il18}.x41,\n x51 = {il4e}.x51,\n x61 = {il27}.x61,\n x71 = {il6e}.x71,\n x81 = {ifnge}.x81,\n x91 = {tcr}.x91,\n x12 = {neg_tgfb}.x12,\n x22 = {neg_il23}.x22,\n x32 = {neg_il12}.x32,\n x42 = {neg_il18}.x42,\n x52 = {neg_il4e}.x52,\n x62 = {neg_il27}.x62,\n x72 = {neg_il6e}.x72,\n x82 = {neg_ifnge}.x82,\n x92 = {neg_tcr}.x92,\n x0 = {}.x0\n]').\n \n") 
                    param_file.write("myentities([neg_foxp3,neg_gata3,neg_ifng,neg_ifngr,neg_il12r,neg_il17,neg_il18r,neg_il2,neg_il21,neg_il21r,neg_il23r,neg_il2r,neg_il4,neg_il4r,neg_il6,neg_il6r,neg_irak,neg_jak1,neg_nfat,neg_nfkb,neg_rorgt,neg_socs1,neg_stat1,neg_stat3,neg_stat4,neg_stat5,neg_stat6,neg_tbet,neg_tgfbr]).\n \n")
                else:
                    print("Invalid version.")
                    return
                param_file.write('mycontext("[' + prolog_context + ']").\n\n')  # x0,
                param_file.write('mytarget([' + state + ']).\n')
                param_file.write('mypos([' + prolog_target + ']).\n')
                param_file.write('myneg([]).')
                param_file.flush()
                param_file.close()
                
                with PrologMQI() as mqi:
                    with mqi.create_thread() as prolog_thread:
                        
                        prolog_thread.query('["BioReSolvePositive.pl"]')
                        result = prolog_thread.query('main_do(ordslice,EKs,ListReactNumbR,ListEpos,ListCS).')
                        print(result)
                        tmp_set = set(chain.from_iterable(result[0]['ListEpos'])) 
                        sliced_computations.append(tmp_set)

                        # for enzyme in enzymes:
                        #     low = enzyme+"_b"
                        #     null = enzyme+"_nul"
                        #     if low in tmp_set:
                        #         pat_enzymes_count[enzyme][0] += 1
                        #     if null in tmp_set:
                        #         pat_enzymes_count[enzyme][1] += 1

                        union_set.update(tmp_set)
                        if (first_time):
                            intersection_set = tmp_set
                            first_time = False
                        else:
                            intersection_set = intersection_set.intersection(tmp_set)
        
        
        for f in glob.glob("tmp-slice*.txt"):
            os.remove(f)
        for f in glob.glob("tmp-legenda*.txt"):
            os.remove(f)
        for f in glob.glob("tmp-slicingrun*.txt"):
            os.remove(f)
        
        print()
        union_set=sorted(union_set)
        intersection_set=sorted(intersection_set)
        print("SET OF ENTITIES IN SLICED COMPUTATIONS FOR TARGET " + str(target) + ":")
        outfile.write("SET OF ENTITIES IN SLICED COMPUTATIONS FOR TARGET " + str(target) + ":\n")
        print("UNION: " + str(union_set))
        outfile.write("         UNION: " + str(union_set)+"\n")
        print("INTERSECTION: " + str(intersection_set))
        outfile.write("         INTERSECTION: " + str(intersection_set)+"\n\n")
        print()
        outfile.flush()
    
        intersection_neg_set = {s[4:] for s in intersection_set if s.startswith("neg")}
        intersection_pos_set = {s for s in intersection_set if not s.startswith("neg")}

        union_neg_set = {s[4:] for s in union_set if s.startswith("neg")}
        union_pos_set = {s for s in union_set if not s.startswith("neg")}
    
        generate_table(str(target.present),intersection_pos_set, union_pos_set, intersection_neg_set, union_neg_set)

        outfile.close()
        sliced_computations_tot = [list(s) for s in sliced_computations]  # JSON cannot serialize sets directly
        with open("Results/Final/slicing/"+dot_name+'/sliced_computations_to_'+str(target.present)+'.json', 'w', encoding='utf-8') as f:
            json.dump(sliced_computations_tot, f, indent=4)
        
    # save_table_to_file("Results/"+dot_name+"/table_results.csv", show_neg=False)

    # save_table_to_file("Results/slicing/"+dot_name+"/table_results_"+str(target.present)+".csv", show_neg=False)

def dynamic_negative_slicing_analysis(version):
    target_tot = len(targets)
    target_count = 0
    
    print("TARGETS TO ANALYZE: " + str(target_tot))
    print()

    #outfile = open("out.txt","w")
    if not os.path.exists("Results/Final/slicing/negative/"+dot_name):
        os.makedirs("Results/Final/slicing/negative/"+dot_name)

    outfile = open("Results/Final/slicing/negative/"+dot_name+"/out_"+str(targets[0].present)+".txt","w")

    for target in targets:
        sliced_computations_pos = []
        sliced_computations_neg = []
        target_count = target_count+1
        print("TARGET COUNT: " + str(target_count) + "/" + str(target_tot))
        contexts, contexts_dict = target_computations(target)

        prolog_target = ','.join(target.present)

        cont = 1
        tot = str(count_states(contexts_dict))
        union_pos_set = set()
        intersection_pos_set = set()
        union_neg_set = set()
        intersection_neg_set = set()
        first_time = True
        for ctx in contexts:
            prolog_context = ','.join(ctx)
            prolog_target_states = contexts_dict[','.join(ctx)]
            for i,state in enumerate(prolog_target_states):
                print("TEST CASE: " + str(cont) + "/" + str(tot))
                cont=cont+1
                print("CONTEXT: " + prolog_context + "      ATTRACTOR STATE: " + str(i+1) + "/" + str(len(prolog_target_states)))
                print("STATE: " + state)
                param_file = open("params.pl",'w')
                if (version=="original"):
                    param_file.write("myenvironment('[bh4_p_b = {bh4_b}.bh4_p_b1,\n  bh4_p_b1 = {xbh4_b}.bh4_p_b1,\n  bh4_p_n = {bh4}.bh4_p_n1,\n  bh4_p_n1 = {xbh4}.bh4_p_n1,\n  phe_p_b = {phe_b}.phe_p_b1,\n  phe_p_b1 = {xphe_b}.phe_p_b1,\n  phe_p_n = {phe}.phe_p_n1,\n  phe_p_n1 = {xphe}.phe_p_n1,\n tat_p_nul = {tat_nul}.tat_p_nul1,\n tat_p_nul1 = {xtat_nul}.tat_p_nul1,\n tat_p_b = {tat_b}.tat_p_b1,\n  tat_p_b1 = {xtat_b}.tat_p_b1,\n  tat_p_n = {tat}.tat_p_n1,\n  tat_p_n1 = {xtat}.tat_p_n1,\n  plp_p_b = {plp_b}.plp_p_b1,\n  plp_p_b1 = {xplp_b}.plp_p_b1,\n  plp_p_n = {plp}.plp_p_n1,\n  plp_p_n1 = {xplp}.plp_p_n1,\n dhpr_p_nul = {dhpr_nul}.dhpr_p_nul1,\n dhpr_p_nul1 = {xdhpr_nul}.dhpr_p_nul1,\n dhpr_p_b = {dhpr_b}.dhpr_p_b1,\n  dhpr_p_b1 = {xdhpr_b}.dhpr_p_b1,\n  dhpr_p_n = {dhpr}.dhpr_p_n1,\n  dhpr_p_n1 = {xdhpr}.dhpr_p_n1,\n pah_fe2_p_nul = {pah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_nul1 = {xpah_fe2_nul}.pah_fe2_p_nul1,\n pah_fe2_p_b = {pah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_b1 = {xpah_fe2_b}.pah_fe2_p_b1,\n  pah_fe2_p_n = {pah_fe2}.pah_fe2_p_n1,\n pah_fe2_p_n1 = {xpah_fe2}.pah_fe2_p_n1,\n gstz1_p_nul = {gstz1_nul}.gstz1_p_nul1,\n gstz1_p_nul1 = {xgstz1_nul}.gstz1_p_nul1,\n gstz1_p_b = {gstz1_b}.gstz1_p_b1,\n  gstz1_p_b1 = {xgstz1_b}.gstz1_p_b1,\n  gstz1_p_n = {gstz1}.gstz1_p_n1,\n  gstz1_p_n1 = {xgstz1}.gstz1_p_n1,\n hpdfe_p_nul = {hpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_nul1= {xhpdfe_nul}.hpdfe_p_nul1,\n hpdfe_p_b = {hpdfe_b}.hpdfe_p_b1,\n  hpdfe_p_b1 = {xhpdfe_b}.hpdfe_p_b1,\n hpdfe_p_n = {hpdfe}.hpdfe_p_n1,\n  hpdfe_p_n1 = {xhpdfe}.hpdfe_p_n1,\n hgdfe_p_nul = {hgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_nul1= {xhgdfe_nul}.hgdfe_p_nul1,\n hgdfe_p_b = {hgdfe_b}.hgdfe_p_b1,\n  hgdfe_p_b1 = {xhgdfe_b}.hgdfe_p_b1,\n hgdfe_p_n = {hgdfe}.hgdfe_p_n1,\n  hgdfe_p_n1 = {xhgdfe}.hgdfe_p_n1,\n fahmgca_p_nul = {fahmgca_nul}.fahmgca_p_nul1,\n fahmgca_p_nul1= {xfahmgca_nul}.fahmgca_p_nul1,\n fahmgca_p_b = {fahmgca_b}.fahmgca_p_b1,\n fahmgca_p_b1 = {xfahmgca_b}.fahmgca_p_b1,\n fahmgca_p_n = {fahmgca}.fahmgca_p_n1,\n fahmgca_p_n1 = {xfahmgca}.fahmgca_p_n1 \n]').\n \n")  
                    param_file.write("myentities([h,o2,h2o,bh2, alphakg, nadh, gsh, ala, arat]).\n \n")
                elif (version=="positive"):
                    param_file.write("myenvironment('[ x1 = ({tgfb}.x11 + {neg_tgfb}.x12),\n x2 = ({il23}.x21 + {neg_il23}.x22),\n x3 = ({il12}.x31 + {neg_il12}.x32),\n x4 = ({il18}.x41 + {neg_il18}.x42),\n x5 = ({il4e}.x51 + {neg_il4e}.x52),\n x6 = ({il27}.x61 + {neg_il27}.x62),\n x7 = ({il6e}.x71 + {neg_il6e}.x72),\n x8 = ({ifnge}.x81 + {neg_ifnge}.x82),\n x9 = ({tcr}.x91 + {neg_tcr}.x92),\n x11 = {tgfb}.x11,\n x21 = {il23}.x21,\n x31 = {il12}.x31,\n x41 = {il18}.x41,\n x51 = {il4e}.x51,\n x61 = {il27}.x61,\n x71 = {il6e}.x71,\n x81 = {ifnge}.x81,\n x91 = {tcr}.x91,\n x12 = {neg_tgfb}.x12,\n x22 = {neg_il23}.x22,\n x32 = {neg_il12}.x32,\n x42 = {neg_il18}.x42,\n x52 = {neg_il4e}.x52,\n x62 = {neg_il27}.x62,\n x72 = {neg_il6e}.x72,\n x82 = {neg_ifnge}.x82,\n x92 = {neg_tcr}.x92,\n x0 = {}.x0\n]').\n \n") 
                    param_file.write("myentities([neg_foxp3,neg_gata3,neg_ifng,neg_ifngr,neg_il12r,neg_il17,neg_il18r,neg_il2,neg_il21,neg_il21r,neg_il23r,neg_il2r,neg_il4,neg_il4r,neg_il6,neg_il6r,neg_irak,neg_jak1,neg_nfat,neg_nfkb,neg_rorgt,neg_socs1,neg_stat1,neg_stat3,neg_stat4,neg_stat5,neg_stat6,neg_tbet,neg_tgfbr]).\n \n")
                else:
                    print("Invalid version.")
                    return
                param_file.write('mycontext("[' + prolog_context + ']").\n\n')  
                param_file.write('mytarget([' + state + ']).\n')
                param_file.write('mypos([' + prolog_target + ']).\n')
                param_file.write('myneg([]).')
                param_file.flush()
                param_file.close()
                print()
                with PrologMQI() as mqi:
                    with mqi.create_thread() as prolog_thread:
                        prolog_thread.query('["BioReSolvePositive.pl"]')
                        
                        result = prolog_thread.query('main_do(negslice,EKs,ListReactNumbR,ListEpos,ListEneg,ListCS).')
                        tmp_pos_set = set(chain.from_iterable(result[0]['ListEpos'])) 
                        sliced_computations_pos.append(tmp_pos_set)
                        union_pos_set.update(tmp_pos_set)
                        tmp_neg_set = set(chain.from_iterable(result[0]['ListEneg'])) 
                        sliced_computations_neg.append(tmp_neg_set)
                        union_neg_set.update(tmp_neg_set)
                        if (first_time):
                            intersection_pos_set = tmp_pos_set
                            intersection_neg_set = tmp_neg_set
                            first_time = False
                        else:
                            intersection_pos_set = intersection_pos_set.intersection(tmp_pos_set)
                            intersection_neg_set = intersection_neg_set.intersection(tmp_neg_set)
        for f in glob.glob("tmp-slice*.txt"):
            os.remove(f)
        for f in glob.glob("tmp-legenda*.txt"):
            os.remove(f)
        for f in glob.glob("tmp-slicingrun*.txt"):
            os.remove(f)
        
        print()
        union_pos_set=sorted(union_pos_set)
        intersection_pos_set=sorted(intersection_pos_set)
        union_neg_set=sorted(union_neg_set)
        intersection_neg_set=sorted(intersection_neg_set)
        print("SET OF ENTITIES IN SLICED COMPUTATIONS FOR TARGET " + str(target) + ":")
        outfile.write("SET OF ENTITIES IN SLICED COMPUTATIONS FOR TARGET " + str(target) + ":\n")
        print("UNION POSITIVE: " + str(union_pos_set))
        outfile.write("         UNION POSITIVE: " + str(union_pos_set)+"\n")
        print("INTERSECTION POSITIVE: " + str(intersection_pos_set))
        outfile.write("         INTERSECTION POSITIVE: " + str(intersection_pos_set)+"\n\n")
        print("UNION NEGATIVE: " + str(union_neg_set))
        outfile.write("         UNION NEGATIVE: " + str(union_neg_set)+"\n")
        print("INTERSECTION NEGATIVE: " + str(intersection_neg_set))
        outfile.write("         INTERSECTION NEGATIVE: " + str(intersection_neg_set)+"\n\n")
        print()
        outfile.flush()
    
        generate_table(str(target.present),intersection_pos_set, union_pos_set, intersection_neg_set, union_neg_set)
        sliced_computations_tot_pos = [list(s) for s in sliced_computations_pos]  # JSON cannot serialize sets directly


        with open("Results/Final/slicing/negative/"+dot_name+'/pos_sliced_computations_to_'+str(target.present)+'_absent'+str(target.absent)+'.json', 'w', encoding='utf-8') as f:
            json.dump(sliced_computations_tot_pos, f, indent=4)
        sliced_computations_tot_neg = [list(s) for s in sliced_computations_neg]  # JSON cannot serialize sets directly
        with open("Results/Final/slicing/negative/"+dot_name+'/neg_sliced_computations_to_'+str(target.present)+'_absent'+str(target.absent)+'.json', 'w', encoding='utf-8') as f:
            json.dump(sliced_computations_tot_neg, f, indent=4)


    outfile.close()
    save_table_to_file("table_results_neg.csv")
    


dot_name= "hgd0_04-02"
G = nx.DiGraph(nx.nx_pydot.read_dot("output/"+dot_name+".dot"))
print(len(G.nodes))
print(len(G.edges))
for target in targets:
    contexts, contexts_dict = target_computations(target)
    if (len(contexts)>0):
        df = pd.DataFrame(contexts).sort_values(by=[9,8,7,6,5,4,3,2,1,0],ignore_index=True)
        print("CONTEXTS THAT LEAD TO THE TARGET:")
        print(df)
        print()
        if not os.path.exists("Results/Final/context_to_target/"+dot_name):
            os.makedirs("Results/Final/context_to_target/"+dot_name)
        df.to_csv("Results/Final/context_to_target/"+dot_name+"/contexts_to_" + str(target.present) + ".csv",index=False,header=False)

table_data_df = pd.DataFrame()
dynamic_negative_slicing_analysis("original")




