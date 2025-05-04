# ncbi-cluster-tracker

**ncbi-cluster-tracker** is a tool for creating static, shareable HTML reports for tracking SNP clusters within the [NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) system. Given an input sample sheet CSV containing BioSample IDs for isolates of interest (referred to as "internal isolates"), the tool creates a report with a high-level overview of all of the clusters associated with the internal isolates. For each cluster the output report links to the corresponding NCBI Pathogen Detection tree and displays additional visualizations such as a pairwise SNP distance matrix heatmap. Any additional metadata or alternate IDs provided in the sample sheet are used to further annotate the internal isolates within the report. If provided a directory to previous outputs, ncbi-cluster-tracker will compare to the previous report and indicate any new isolates, new clusters, and changes to isolate counts for existing clusters.

### Clusters tab
The Clusters tab contains a table of clusters and their associated isolate counts and any changes to isolate counts from the previous report. Columns in the table can be sorted and filtered. A timeline chart is displayed below the table showing when each isolate was added to each cluster. Details about a particular isolate can be viewed by hovering over a point.

![Clusters tab](assets/clusters.png)

### Isolates tab
The Isolates tab contains a table with details about the internal and external isolates, such as which cluster (if any) they belong to, whether they are new, their BioSample metadata, and any custom metadata provided in the sample sheet. SQL can be used directly in the report to query this table and the Clusters table.

![Isolates tab](assets/isolates.png)

### Cluster details tab
The Cluster details tab displays information about each cluster. Clusters can be searched and selected in the dropdown menu. Each page displays the count (and change in count) of internal and external isolates in the cluster and links are provided to the SNP Tree Viewer in the NCBI Pathogen Detection website. A downloadable custom labels file is provided to annotate the tree in the SNP Tree viewer with any alternate sample IDs and other metadata.

![Cluster details overview](assets/cluster_details_overview.png)

This section of the Cluster details page shows a SNP distance matrix. If there are a large number of isolates in the cluster, the matrix is filtered to only show the nearest external isolates.

![Cluster details SNP distance matrix](assets/cluster_details_matrix.png)

This section of the Cluster details page displays a histogram showing the isolate counts over time (by record creation date).

![Cluster details histogram](assets/cluster_details_histogram.png)

## Requirements
There are two options to provide the required cluster metadata and isolate metadata from Pathogen Detection:

 1. The default option is to use NCBI's public Google BigQuery [`pdbrowser`](https://www.ncbi.nlm.nih.gov/pathogens/docs/gcp/) dataset to automatically download the most recent isolate and cluster metadata. However, this option requires a Google Cloud Platform account and may not reflect the most up-to-date information in the Pathogen Detection system since there is up to a one day delay for data to be published to the BigQuery dataset. 
2. Alternatively, users can manually provide this data by downloading search results as a TSV or CSV file from the Pathogen Detection Isolates Browser website and passing the file into the program using `--browser-file`. To create this file, copy-and-paste the BioSample IDs from the sample sheet into the [Isolates Browser](https://www.ncbi.nlm.nih.gov/pathogens/) search bar and press Enter. Make sure all columns are displayed (select "Choose columns"), then click the "Download" button to export the Matched Isolates table as a CSV or TSV. This file will contain the isolate metadata for internal isolates and their associated clusters. To get the isolate metadata for external isolates in these clusters, a second search needs to be run. Copy the "SNP cluster" column from the downloaded file and paste into to the search bar along with the existing BioSample IDs from the previous search. Run the search again and download the second TSV or CSV as before, and pass the path to this file to ncbi-cluster-tracker with the `--browser-file` option.

## Installation
1. Install as command-line program with `pip`. Requires Python version >=3.12.
```
pip install ncbi-cluster-tracker
```

2. Ensure logged in to Google Cloud Platform account (if using BigQuery)

```
gcloud auth login
```

## Usage

```
$ ncbi-cluster-tracker --help
usage: ncbi-cluster-tracker [-h] [--out-dir OUT_DIR] [--retry | --no-retry] [--browser-file BROWSER_FILE]
                            [--compare-dir COMPARE_DIR | --no-compare]
                            sample_sheet

positional arguments:
  sample_sheet          Path to sample sheet CSV with required "biosample" column and any additional metadata columns.
                        Use "id" column for alternate isolate IDs.

options:
  -h, --help            show this help message and exit
  --out-dir OUT_DIR, -o OUT_DIR
                        Path to directory to store outputs. Defaults to "./outputs/" if not specified.
  --retry, --no-retry   Do not query BigQuery or NCBI, assumes data has already been downloaded to --out-dir or
                        directory with most recent timestamp.
  --browser-file BROWSER_FILE
                        Path to isolates TSV or CSV downloaded from the Pathogen Detection Isolates Browser with
                        information for all internal and external isolates. When specified, data in file will be used
                        instead of querying the BigQuery dataset.
  --compare-dir COMPARE_DIR
                        Path to previous output directory to detect and report new isolates. Defaults to directory
                        inside --out-dir with most recent timestamp if not specified.
  --no-compare          Do not compare to most recent output directory, all clusters and isolates will be considered
                        "new".
```
