# TSG-Anomaly-Detection

This repository contains code supporting the text "Multiple Network Embedding for Anomaly Detection in Time Series of Graphs"

Guodong Chen, Jes\'{u}s Arroyo, Avanti Athreya, Joshua Cape, Joshua T. Vogelstein, Youngser Park, Chris White,
Jonathan Larson, Weiwei Yang, and Carey E. Priebe

https://arxiv.org/abs/2008.10055

## Abstract
This paper considers the graph signal processing problem of anomaly detection in time series of graphs. We examine two related, complementary inference tasks: the detection of anomalous graphs within a time series, and the detection of temporally anomalous vertices. We approach these tasks via the adaptation of statistically principled methods for joint graph inference, specifically 
multiple adjacency spectral embedding (MASE). We demonstrate that our method is effective for our inference tasks. Moreover, we assess the performance of our method in terms of the underlying nature of detectable anomalies. We further provide the theoretical justification for our method and insight into its use. Applied to the Enron communication graph and a large-scale commercial search engine time series of graphs, our approaches demonstrate their applicability and identify the anomalous vertices beyond just large degree change.

## Results

| Figure(s) | R file | HTML file |
|:--------- |:-----------|:-----------|
| 1  | figure1.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure1.html |
| 2   | figure2.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure2.html |
| 3    | figure3.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure3.html |
| 4   | figure4.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure4.html |
| 5    | figure5.R |
| 6, 8, 9, 10, 11, 12, 13 | figure6_8_9_10_11_12_13.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure6_8_9_10_11_12_13.html |
| 7, 14, 15, 16, 17 | figure7_14_15_16_17_table1_table2.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure7_14_15_16_17_table1_table2.html |

| Table | R file | HTML file |
|:----- |:-----------|:-----------|
| 1   | figure7_14_15_16_17_table1_table2.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure7_14_15_16_17_table1_table2.html |
| 2   | figure7_14_15_16_17_table1_table2.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure7_14_15_16_17_table1_table2.html |
| 3   | figure7_14_15_16_17_table1_table2.R | https://www.cis.jhu.edu/~parky/AnomalyDetection/figure7_14_15_16_17_table1_table2.html |

## Code

| Utility functions | File |
|:----- |:-----------|
| 1   | utils.R |
| 2   | enron_utils.R |
| 3   | ms_utils.R |
| 4   | qcc.R |

## Data

### Enron
https://www.cis.jhu.edu/~parky/Enron/

### MSB
The commercial search engine data used in Section 6 are available from the authors upon reasonable request and with the permission of Microsoft Corporation.
