# Code base
## Tool for EHR & mutation profile based patient clustering & visualization  

This tool was developed in partial fulfillment of the requirements for the course “Medical Informatics” at the University Medical Center Göttingen.

Thesis title : Prototypical Development of a Data Exploration Tool for Cancer Patients based on Patient Similarities

Author: Jan Janosch Schneider

First Supervisor: Prof. Dr. Anne-Christin Hauschild

Second Supervisor: Prof. Dr. Ulrich Sax


## Installation

To use the tool, clone this repository and install the requirements.txt using pip install -r requirements.txt.

## Usage

### Running the tool

To run the tool, execute the 'dash_visualization.py' file.
Once running, go to the IP adress that is shown in the terminal.

### Data integration

To upload your data, the files needs be of a certain format.
The tool requires a .csv file with the columns "age", "age at diagnosis", "sex", "icds","mutations","cohort".
The age values must be numerical, the sex can either be "M"/"Male" or "F/Female".
The ICD-codes have to be space seperated and represent a valid ICD-10 code.
The currently embedded version of the ICD-taxonomy is ICD-10 GM 2022.
The mutations have to be HUGO gene symbols and space seperated as well.
The cohort can be any string you wish to assign to the uploaded cohort.

To use the cBioPortal API for integration of public data sets, an active internet connection is needed.
However, using the "create_local_resources.py" file you can download as many studies from the cBioPortal and save them under resources/local_studies to have them available offline.

### Current limitations

**Wait until the browser tab does not say *Updating..* to avoid deadlocks in the program.**
