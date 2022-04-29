import requests
from hubmap_api_py_client import Client
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

api_endpoint = 'https://cells.api.hubmapconsortium.org/api/'
client = Client(api_endpoint)

#Helper functions for transforming API results into Pandas and Seaborn objects

def get_num_cells(metadata_df, dataset_uuid):
  dataset_df = metadata_df[metadata_df['dataset'] == dataset_uuid]
  return dataset_df['num_cells'][0]

def get_metadata_df_by_cells(metadata_df, attribute):
  cell_count_dict = {value:0 for value in metadata_df[attribute].unique()}
  for value in metadata_df[attribute].unique():
    attribute_df = metadata_df[metadata_df[attribute] == value]
    cell_count_dict[value] += attribute_df['num_cells'].sum()
  records = [{attribute:value, 'num_cells':cell_count_dict[value]} for value in cell_count_dict]
  return pd.DataFrame(records)

def get_metadata_df_binned(metadata_df, attribute, bin_size, by_cells=False):
  cell_count_dict = {value:0 for value in metadata_df[attribute].unique()}
  for value in metadata_df[attribute].unique():
    attribute_df = metadata_df[metadata_df[attribute] == value]
    cell_count_dict[value] += attribute_df['num_cells'].sum()
  records = [{attribute:value, 'num_cells':cell_count_dict[value]} for value in cell_count_dict]
  return pd.DataFrame(records)

def filter_df_by_modality(metadata_df, modality):
  modality_datasets = [dataset["uuid"] for dataset in list(client.select_datasets(where="modality", has = [modality]).get_list())]
  return metadata_df[metadata_df["dataset"].isin(modality_datasets)]

def get_datasets_by_organ(metadata_df, organ, modality = None):
  if modality:
    metadata_df = filter_df_by_modality(metadata_df, modality)
  organ_df = metadata_df[metadata_df["organ"] == organ]
  return list(organ_df["dataset"].unique())

def get_full_organ_name(abbreviation):
  organs_dict = {'HT':'Heart', 'SI': 'Small Intestine', "LK":"Kidney",'RK':"Kidney", "LI":"Large Intestine", "LV":"Liver", "RL":"Lung", "LL":"Lung", "PA":"Pancreas","LY":"Lymph Node", "SP":"Spleen", "TH":"Thymus"}
  return organs_dict[abbreviation]

def get_dataset_json(dataset):
  dataset_query_dict = {
        "query": {
            "bool": {
                "must": [],
                "filter": [
                    {
                        "match_all": {}
                    },
                    {
                        "exists": {
                            "field": "files.rel_path"
                        }
                    },
                    {
                        "match_phrase": {
                            "uuid": {
                                "query": dataset
                            },
                        }

                    }
                ],
                "should": [],
                "must_not": [
                    {
                        "match_phrase": {
                            "status": {
                                "query": "Error"
                            }
                        }
                    }
                ]
            }
        }
    }

  dataset_response = requests.post(
          'https://search.api.hubmapconsortium.org/search',
          json=dataset_query_dict)

  hits = dataset_response.json()['hits']['hits']
  return hits


def get_sample_metadata(dataset_uuid):
  hits = get_dataset_json(dataset_uuid)
  for hit in hits:
    for ancestor in hit['_source']['ancestors']:
        if 'organ' in ancestor.keys():
            return {'organ': get_full_organ_name(ancestor['organ'])}

def get_donor_metadata(dataset_uuid):
  donor_metadata_dict = {}
  properties = ['Age', 'Race', 'Sex', 'Body mass index', 'Blood type']
  hits = get_dataset_json(dataset_uuid)
  if len(hits) == 0:
    print(f"No hits found for uuid {dataset_uuid}")
  for hit in hits:
    if len(hit['_source']['ancestors']) == 0:
      print(f"No ancestors found for hit in dataset {dataset_uuid}")
    for ancestor in hit['_source']['ancestors']:
      donor_metadata_dict = {}
      if 'lab_donor_id' in ancestor.keys():
        try:
          if 'metadata' in ancestor.keys():
            if 'organ_donor_data' in ancestor['metadata'].keys():
              organ_donor_data = ancestor['metadata']['organ_donor_data']
            else:
              organ_donor_data = ancestor['metadata']['living_donor_data']
            donor_metadata_dict = {prop:get_donor_property(prop, organ_donor_data) for prop in properties}
        except:
          donor_metadata_dict = {}
        return donor_metadata_dict

  return donor_metadata_dict

def get_donor_property(prop, donor_metadata):
  if donor_metadata is None:
    return "Unknown"
  for attribute in donor_metadata:
    if attribute['grouping_concept_preferred_term'] == prop:
      if prop in ['Age', 'Body mass index']:
        return float(attribute['data_value'])
      elif prop in ['Sex', 'Race', 'Blood type']:
        return attribute['preferred_term']

def get_metadata_record(dataset_uuid):
  dataset_cells = client.select_cells(where="dataset", has=[dataset_uuid])
  dataset_clusters = client.select_clusters(where="dataset", has=[dataset_uuid])
  portal_url = get_portal_url(dataset_uuid)
  metadata_dict = {"dataset":dataset_uuid, "portal_url":portal_url, "num_cells":len(dataset_cells), "num_cluster":len(dataset_clusters)}
  metadata_dict.update(get_sample_metadata(dataset_uuid))
  metadata_dict.update(get_donor_metadata(dataset_uuid))
  return metadata_dict

def get_metadata_df(dataset_uuids):
  records = [get_metadata_record(uuid) for uuid in dataset_uuids]
  return pd.DataFrame(records)

def get_barplot_by_column(dataframe, x, y=None):
    if y:
        sns.barplot(data=dataframe, x=x, y=y)
    else:
        sns.countplot(data=dataframe, x=x)

def get_barplot_cells_by_column(dataframe, x, y):
    sns.barplot(data=dataframe, x=x, y=y)

def get_dataset_comparison_dataframe(dataset_one, dataset_two, gene_symbol):
    dataset_one_cells = client.select_cells(where="dataset", has=[dataset_one])
    dataset_two_cells = client.select_cells(where="dataset", has=[dataset_two])
    dataset_cells_list = list(dataset_one_cells.get_list(values_included=[gene_symbol])) + list(
        dataset_two_cells.get_list(values_included=[gene_symbol]))
    values_list = [cell['values'] for cell in dataset_cells_list]
    for i in range(len(values_list)):
        values_list[i]['dataset'] = dataset_cells_list[i]['dataset']
    values_df = pd.DataFrame(values_list)
    return values_df

def get_organ_comparison_dataframe(metadata_df, organ_one, organ_two, gene_symbol):
    organ_one_datasets = get_datasets_by_organ(metadata_df, organ_one, modality="rna")
    organ_two_datasets = get_datasets_by_organ(metadata_df, organ_two, modality="rna")
    # These steps involve bulk retrieval of many datasets, so they take a long time
    dataset_uuids = organ_one_datasets + organ_two_datasets
    dataset_cells_list = []
    for dataset in dataset_uuids:
        dataset_cells = client.select_cells(where="dataset", has=[dataset])
        dataset_cells_list.extend(list(dataset_cells.get_list(values_included=[gene_symbol])))
    values_list = [cell['values'] for cell in dataset_cells_list]
    for i in range(len(values_list)):
        values_list[i]['organ'] = dataset_cells_list[i]['organ']
    values_df = pd.DataFrame(values_list)
    return values_df

def get_coexpression_scatterplot(values_df, gene_one, gene_two):
    sns.scatterplot(data=values_df, x=gene_one, y=gene_two)

def get_comparison_histogram(values_df, gene_symbol, hue):
    sns.histplot(data=values_df, x=gene_symbol, hue=hue, log_scale=(False, True))

def make_histogram(metadata_df, x, binwidth=10):
    sns.histplot(data=metadata_df, x=x, binwidth=binwidth)

def get_dataset_comparison_table(dataset_one, dataset_two, gene_symbol):
    dataset_one_cells = client.select_cells(where="dataset", has=[dataset_one])
    dataset_two_cells = client.select_cells(where="dataset", has=[dataset_two])
    dataset_cells_list = list(dataset_one_cells.get_list(values_included=[gene_symbol])) + list(
        dataset_two_cells.get_list(values_included=[gene_symbol]))
    values_list = [cell['values'] for cell in dataset_cells_list]
    for i in range(len(values_list)):
        values_list[i]['dataset'] = dataset_cells_list[i]['dataset']
    values_df = pd.DataFrame(values_list)
    return values_df

def get_organ_comparison_table(metadata_table, organ_one, organ_two, gene_symbol):
    organ_one_datasets = get_datasets_by_organ(metadata_table, organ_one, modality="rna")
    organ_two_datasets = get_datasets_by_organ(metadata_table, organ_two, modality="rna")
    # These steps involve bulk retrieval of many datasets, so they take a long time
    dataset_uuids = organ_one_datasets + organ_two_datasets
    dataset_cells_list = []
    for dataset in dataset_uuids:
        dataset_cells = client.select_cells(where="dataset", has=[dataset])
        dataset_cells_list.extend(list(dataset_cells.get_list(values_included=[gene_symbol])))
    values_list = [cell['values'] for cell in dataset_cells_list]
    for i in range(len(values_list)):
        values_list[i]['organ'] = dataset_cells_list[i]['organ']
    values_df = pd.DataFrame(values_list)
    return values_df

def get_metadata_table_by_cells(metadata_table, attribute):
  cell_count_dict = {value:0 for value in metadata_table[attribute].unique()}
  for value in metadata_table[attribute].unique():
    attribute_df = metadata_table[metadata_table[attribute] == value]
    cell_count_dict[value] += attribute_df['num_cells'].sum()
  records = [{attribute:value, 'num_cells':cell_count_dict[value]} for value in cell_count_dict]
  return pd.DataFrame(records)

def get_metadata_table_binned(metadata_table, attribute, bin_size, by_cells=False):
  cell_count_dict = {value:0 for value in metadata_table[attribute].unique()}
  for value in metadata_table[attribute].unique():
    attribute_df = metadata_table[metadata_table[attribute] == value]
    cell_count_dict[value] += attribute_df['num_cells'].sum()
  records = [{attribute:value, 'num_cells':cell_count_dict[value]} for value in cell_count_dict]
  return pd.DataFrame(records)

def filter_table_by_modality(metadata_table, modality):
  modality_datasets = [dataset["uuid"] for dataset in list(client.select_datasets(where="modality", has = [modality]).get_list())]
  return metadata_table[metadata_table["dataset"].isin(modality_datasets)]

def get_metadata_table(dataset_uuids):
  records = [get_metadata_record(uuid) for uuid in dataset_uuids]
  return pd.DataFrame(records)

def get_metadata_attributes(metadata_table):
  return list(metadata_table.columns)

def get_distinct_values_for_attribute(metadata_table, attribute):
  return list(metadata_table[attribute].unique())

def select_from_table_by_metadata_value(metadata_table, attribute, value):
  return metadata_table[metadata_table[attribute] == value]

def get_portal_url(uuid):
  return f"https://portal.hubmapconsortium.org/browse/dataset/{uuid}/"

documentation_dict = {get_metadata_table:"""def get_metadata_table()
  This function doesn't take any inputs
  It returns a table with sample and donor metadata for each dataset in the Cells API index""",
  
  get_metadata_table_by_cells:"""def get_metadata_table_by_cells(metadata_table)
  This function takes a metadata table as input
  It returns the table with counts by cells instead of by datasets""",
  
  get_barplot_by_column:"""def get_barplot_by_column()
  This function takes a metadata table as input, along with x and y axis labels for a barplot
  It outputs a barplot""",
  
  get_coexpression_scatterplot:"""def get_coexpression_scatterplot()
  This function takes a table of cells as input, along with identifiers for two genes whose expression to compare
  It outputs a scatterplot""",

  get_comparison_histogram:"""def get_comparison_histogram()
  This function takes a table of cells, along with a gene identifier and the name of the attribute to color by
  It produces a histogram comparing the expression of a certain gene across two or more groups of cells, colored by group""",

  get_organ_comparison_table:"""def get_organ_comparison_table()
  This function takes the identifiers for two organ and a gene and returns a table of cells from those datasets with expression
  values for that gene""",

  get_dataset_comparison_table: """def get_dataset_comparison_table()
  This function takes the identifiers for two datasets and a gene and returns a table of cells from those datasets with expression
  values for that gene""",

  get_datasets_by_organ: """def get_datasets_by_organ()
  This function takes an identifier for an organ and returns a list of identifiers for datasets associated with that organ""",

  get_metadata_attributes:"""def get_metadata_attributes(metadata_table)
  This function takes a metadata table and displays all the donor, sample, and dataset metadata attributes available""",

  get_distinct_values_for_attribute:"""def get_distinct_values_for_attribute()
  This function takes a metadata table and a metadata attribute and returns all the distinct values for that attribute in that table""",

  select_from_table_by_metadata_value:"""def select_from_metadata_table_by_value():
  This function takes a metadata table, a metadata attribute, and a value for that attribute and returns the subset of the table
  Where that attribute matches that value, i.e. where donor sex = Female""",

  get_portal_url:"""def get_portal_url():
  This function takes a dataset uuid and returns a url to view that dataet in the HuBMAP data portal"""}

def get_function_help_text(function):
  if function not in documentation_dict:
      return "No documentation available for that function"
  else:
      return documentation_dict[function]