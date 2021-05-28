# HOST-PATHOGEN PROTEIN-PROTEIN INTERACTION PREDICTOR 
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

...Under construction...

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



### AUTHORS/CONTRIBUTORS

...Under construction...


### LICENSE

...Under construction...
