# AlphaMut![Logo](src/LOGO.png)

A web application and analysis pipeline for searching and visualizing folded proteins of human genes taking account different mutations.

## Description

This project provides a platform for users to search for human genes and view detailed information about their folded protein structures. Utilizing interactive search features and a curated database of folded proteins, users can explore gene data, protein structures, and collaborate on further research.

This is only a prototype proof of concept. All the foldings were done in my laptop using [local colabfold](https://github.com/YoshitakaMo/localcolabfold). I would love to set this up in a proper server and running foldings constantly using proper hardware.

If you enjoyed this project consider checking out its sister project, [AlphaSplice](https://github.com/DrAnomalocaris/AlphaSplice), A similar project which is focused on splice variants rather than mutations.

As this is a prototype, there is only a handful of genes available. If you are particularly interested in a specific gene, feel free to contact me to have said gene processed.

## Features

Front end features:

* **Interactive Search Bar** : Autocomplete suggestions for genes and transcripts.
* **Available Folded Proteins List** : Quick access to a list of genes with available folded protein data.
* **Gene Information Pages** : Detailed pages displaying gene-specific data and protein structures.
* **Responsive Design** : User-friendly interface compatible with various devices.
* **Collaboration Opportunities** : Contact information for users interested in contributing.

Back end features:

* **Pipeline:** Analysis pipeline in snakemake format
* **Ease of use:** To process a specific gene, all that is needed is to request the output analysis file using snakemake, and snakemake will get the gene data and process it.

## Usage

* **Search for a Gene** : Use the search bar at the top of the page to find a gene or transcript. Autocomplete suggestions will help you find the correct gene.
* **Gene Pages** : Explore the gene-specific pages to view protein structures and associated data.
* **View Folded Proteins** : Scroll down on the homepage to see a list of available folded proteins. Click on any gene to view detailed information.

## Road map

* [ ] Add support for different species (as I use [clinvar](https://www.ncbi.nlm.nih.gov/clinvar/intro/), and it only has human data, adding another species would require modifying a lot of code)
* [ ] Add a daemon to tally requested genes, to allow for the homepage to display the most popular genes and to prioritize gene folding for most requested genes that are not been processed yet.
* [ ] Develop an API or creation of a large file to allow users to download bulk data.
* [ ] Change to locally running MSA, as this runs on a laptop, I do not have enough storage to download the required files for running local MSA.
* [ ] Implement a proper API which uses SQL in the backend, as this prototype is designed to be hosted in github pages, there were some limitations. However a proper implementation will have a proper SQL based backend.

## Contact

If you're interested in collaborating or have any questions, feel free to reach out:

* **Email** : [lahat.albert@gmail.com](mailto:lahat.albert@gmail.com)
* **GitHub** : [DrAnomalocaris](https://github.com/DrAnomalocaris/)
* **LinkedIn** : [Dr. Albert Lahat](https://www.linkedin.com/in/dr-albert-lahat-85b30578/)

## License

This project is licensed under the [CC-NC ](https://creativecommons.org/licenses/by-nc/4.0/deed.en)License.
