# MethylSeek: DNA Methylation Analysis Tool

MethylSeek is a comprehensive R-based tool designed for analyzing DNA methylation patterns using Illumina Human Methylation EPIC array data. This tool facilitates preprocessing of raw methylation data, performs differential methylation analysis using `limma`, and generates visualizations such as volcano plots, MA plots, and heatmaps to explore and interpret methylation changes associated with diseases or conditions.

## Features

- **Data Preprocessing:** Includes quality control, normalization, and filtering of DNA methylation data.
- **Differential Methylation Analysis:** Utilizes `limma` for statistical analysis to identify differentially methylated regions (DMRs).
- **Visualization:** Generates informative plots (volcano plots, MA plots, heatmaps) for visualizing methylation patterns and differential analysis results.
- **Flexible:** Customizable for different experimental designs and datasets.

## Requirements

- R 
- R packages: minfi, limma, ggplot2, IlluminaHumanMethylationEPICanno.ilm10b2.hg19

## Usage

1. **Clone the repository:**
   ```bash
   git clone https://github.com/username/MethylSeek-DNA-Methylation-Analysis-Tool.git
   cd MethylSeek-DNA-Methylation-Analysis-Tool

   Run the analysis script:

bash
Copy code
Rscript methylation_analysis.R /path/to/your/idat/files /path/to/output/directory
Replace /path/to/your/idat/files with the directory containing IDAT files and /path/to/output/directory with the desired output directory.

Output:

Plots and analysis results will be saved in the specified output directory.
Check the README.md file in the output directory for a summary of results.
License
This project is licensed under the MIT License. See the LICENSE file for details.

Contributing
Contributions are welcome! Please fork the repository and create a pull request for any improvements or features.

Authors
Balaji
Contact
For questions or feedback, please contact balajimadhvn@gmail.com.
