import pandas as pd

from google.cloud import bigquery


def query_set_of_clusters(biosamples: list[str]) -> list[str]:
    biosamples_str = python_list_to_sql_str(biosamples)
    query = f'''
    SELECT DISTINCT erd_group AS cluster
    FROM `ncbi-pathogen-detect.pdbrowser.isolates`
    WHERE biosample_acc IN ({biosamples_str})
    '''
    df = execute_query(query)
    clusters = df['cluster'].to_list()
    return clusters

def query_isolates(clusters: list[str], biosamples: list[str]) -> pd.DataFrame:
    clusters_str = python_list_to_sql_str(clusters)
    biosamples_str = python_list_to_sql_str(biosamples)
    query = f'''
    SELECT
        isolate_identifiers[0] AS isolate_id,
        biosample_acc AS biosample,
        target_acc,
        erd_group AS cluster,
        Run AS sra_id,
        isolation_source,
        geo_loc_name,
        collection_date,
        creation_date,
        taxgroup_name,
        scientific_name,
        bioproject_acc,
        STRING_AGG(a.element, ",") AS beta_lactamase_genes
    FROM `ncbi-pathogen-detect.pdbrowser.isolates`
    LEFT JOIN UNNEST(AMR_genotypes) AS a 
    WHERE erd_group IN ({clusters_str})
    OR biosample_acc IN ({biosamples_str}) 
    GROUP BY
        isolate_id,
        biosample,
        target_acc,
        cluster,
        sra_id,
        isolation_source,
        geo_loc_name,
        collection_date,
        creation_date,
        taxgroup_name,
        scientific_name,
        bioproject_acc
    ORDER BY isolate_id;
    '''
    df = execute_query(query)
    return df


def query_clusters(biosamples: list[str]) -> pd.DataFrame:
    biosamples_str = python_list_to_sql_str(biosamples)
    query = f'''
    SELECT DISTINCT
        cluster_isolates.erd_group AS cluster,
        cluster_size.num AS total_count,
        taxgroup_name,
        earliest_added,
        latest_added,
        earliest_year_collected,
        latest_year_collected
    FROM
    (
        SELECT
            erd_group,
            taxgroup_name
        FROM `ncbi-pathogen-detect.pdbrowser.isolates`
        WHERE biosample_acc IN ({biosamples_str})
        GROUP BY
            erd_group,
            taxgroup_name
    ) AS cluster_isolates
    LEFT JOIN
    (
        SELECT
            erd_group,
            COUNT(*) AS num,
            MIN(SUBSTRING(creation_date, 0, 10)) AS earliest_added,
            MAX(SUBSTRING(creation_date, 0, 10)) AS latest_added,
            MIN(SUBSTRING(collection_date, 0, 4)) AS earliest_year_collected,
            max(substring(collection_date, 0, 4)) as latest_year_collected
        FROM `ncbi-pathogen-detect.pdbrowser.isolates`
        GROUP BY erd_group
    ) AS cluster_size
    ON cluster_isolates.erd_group = cluster_size.erd_group
    WHERE cluster_size.num IS NOT NULL
    ORDER BY total_count DESC
    '''
    df = execute_query(query)
    return df

def python_list_to_sql_str(python_list: list[str]) -> str:
    quoted_list = [f'"{l}"' for l in python_list]
    sql_str = ', '.join(quoted_list)
    return sql_str

def execute_query(query: str) -> pd.DataFrame:
    client = bigquery.Client()
    df = client.query_and_wait(query).to_dataframe()
    return df
