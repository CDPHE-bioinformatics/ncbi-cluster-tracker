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
        bioproject_acc
    FROM `ncbi-pathogen-detect.pdbrowser.isolates`
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
    df['collection_date'] = df['collection_date'].astype('string')
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
            MAX(SUBSTRING(collection_date, 0, 4)) as latest_year_collected
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

def isolates_df_from_browser_df(
    browser_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Parse DataFrame from --browser-tsv into `isolates_df` DataFrame format.
    The output format should match that from query_isolates().
    """
    browser_df['isolate_id'] = browser_df['Isolate identifiers'].str.split(',').str[0]
    rename_cols = {
        'isolate_id': 'isolate_id',
        'BioSample': 'biosample',
        'Isolate': 'target_acc',
        'SNP cluster': 'cluster',
        'Run': 'sra_id',
        'Isolation source': 'isolation_source',
        'Location': 'geo_loc_name',
        'Collection date': 'collection_date',
        'Create date': 'creation_date',
        '#Organism group': 'taxgroup_name',
        'Scientific name': 'scientific_name',
        'BioProject': 'bioproject_acc'
    }
    df = browser_df[rename_cols.keys()].rename(columns=rename_cols)
    df['collection_date'] = df['collection_date'].astype('string')
    df.to_csv('test_out.csv', index=False)
    return df 

def cluster_df_from_isolates_df(
    isolates_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Return DataFrame with counts and date ranges for each cluster in
    `isolates_df`.
    """
    isolates_df = isolates_df.copy()
    cluster_isolates = (
        isolates_df
        .groupby(['cluster', 'taxgroup_name'], as_index=False)
        .first()
    )
    isolates_df['creation_date'] = pd.to_datetime(
        isolates_df['creation_date'], errors='coerce'
    )
    isolates_df['collection_date'] = isolates_df['collection_date'].astype('string')
    isolates_df['collection_year'] = isolates_df['collection_date'].str[:4]
    cluster_size = (
        isolates_df
        .groupby('cluster')
        .agg(
            total_count=('cluster', 'count'),
            earliest_added=(
                'creation_date',
                lambda x: x.dropna().min().strftime('%Y-%m-%d'),
            ),
            latest_added=(
                'creation_date',
                lambda x: x.dropna().max().strftime('%Y-%m-%d'),
            ),
            earliest_year_collected=(
                'collection_year',
                lambda x: x.dropna().min(),
            ),
            latest_year_collected=(
                'collection_year',
                lambda x: x.dropna().max(),
            )
        )
        .reset_index()
    )

    df = (
        cluster_isolates
        .merge(cluster_size, on='cluster', how='left')
        .dropna(subset=['total_count'])
        .sort_values(by='total_count', ascending=False)
    )
    keep_cols = [
        'cluster',
        'total_count',
        'taxgroup_name',
        'earliest_added',
        'latest_added',
        'earliest_year_collected',
        'latest_year_collected'
    ]
    df = df[keep_cols] 
    df = df.astype({
            'earliest_year_collected': 'string',
            'latest_year_collected': 'string',
        })
    return df