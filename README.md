# *Eureka!* custom scripts
A collection of miscellaneous scripts connected to the Eureka! data reduction 
pipeline for JWST observations. Excluding some exceptions such as the changed
source files in `eureka_nirspec_changes`, I'll try to stick to a common naming
scheme for custom plotting routines, data evaluation runs etc., meaning that
these files will start with `eucus_*`. This will make it easier for me to keep
track of them.

- `EUREKA_FRESH_PROJECT.sh`: This bash-script can be placed in the *Eureka!* 
top-directory to automate the creation of new data reduction projects.


## `eureka_custom_plots`
A collection of customised plotting routines for the data output of various
stages in *Eureka!*.
- `data/`: Storage for data input into custom routines
- `plots/*`: Storage for custom plots
- All custom routines start with an **Overview** header that describes the
plotting script


## `eureka_environment`
The two files in here should help with the creation of a working Python 
environment for *Eureka!* (which should probably be created within a venv).
- `eureka_packages.txt`: Provides a list of all installed Python packages 
and their version numbers from the working environment for me.
- `eureka_packages.yml`: Environment-file that can be used to create a working
*Eureka!* environment with conda


## `eureka_nirspec_changes`
Two NIRSpec related files that need to be changed when looking at real JWST 
data. Their location is `code_source/eureka/src/eureka/S3_data_reduction/*`
- `nirspec.py`: The `flag_bg()` function needs to be adjusted (see also [GitHub 
Issue #193](https://github.com/kevin218/Eureka/issues/193)). It would 
be interesting to note once this gets incorporated into the source code. For 
now, it sets the background masking back form 0's to the routine used in 
`nircam.py`
- `s3_reduce.py`: Added an additional log-message noting that the file above 
has been changed
