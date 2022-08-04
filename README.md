# HOST-PATHOGEN PROTEIN-PROTEIN INTERACTION PREDICTOR (HPIPred)
This tool allows for the prediction of putative host-pathogen protein-protein interactions based on numerical encoding of physicochemical descriptors


### GETTING STARTED

#### Pre-requisites

This tool has been tested on Ubuntu (18.04) and MacOS Catalina.
Please, make sure to install local blast, R and Zenity.

##### Ubuntu:
Install blast locally: ```bash sudo apt-get install ncbi-blast+ ```

Install R: ```bash sudo apt install r-base ```

Install Zenity: ```bash sudo apt install zenity ```


##### MacOS:

To install blast locally follow the instructions in: https://www.ncbi.nlm.nih.gov/books/NBK52640/

Install R package in: https://cran.r-project.org/bin/macosx/

Install Zenity*: brew install zenity

* Please note that you need homebrew installed in your computer. To install homebrew for mac follow the instructions in: https://brew.sh


#### Installation

Please, download the repository as a zip file and uncompress it

Download dataset/ and database/ directories* from Zenodo into the main directory where the repository is allocated in your system. Then, unzip the files.


*Dataset and database directories are available under: https://zenodo.org/record/4668840#.YHm-B-0zYVs

### RUNNING EXAMPLE

In order to run the tool, go to /scripts inside the main directory and run ```./main_script.bash ```

After this, Zenity interface will be prompted and the user will be asked to choose:
1) A host organism (either from the available database in the repository, from uniprot via its taxon ID or from a custom file) 
2) A pathogen organism
3) From a set of physicochemical descriptors (either custom or default) to perform the predictive analysis
4) The False discovery rate (this is a measure of astringency, the lower the chosen value, the more astringent the analysis will be)
5) The percentage of physicochemical models that have to agree on a prediction in order to report it in the consensus interactome

### LICENSE

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
