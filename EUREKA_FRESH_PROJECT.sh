# USER VARIABLES (PLEASE EDIT)
PROJECT_NAME="trappist1"
INSTRUMENT="nirspec" # TODO: Add functionality for other instruments
RAW_DATA_DIR="RELATIVE/TO/FILE"
PROJECT_DIR="RELATIVE/TO/FILE"


#######################################################
#                                                     #
#          LOOK BELOW AT YOUR OWN RISK                #
#                                                     #
#######################################################

# STATIC
PROJECT_TOPDIR=`pwd`
TEMPLATE_DIR="$PROJECT_TOPDIR/RAW/demos"

# CHECK IF TOP DIRECTORY IS CORRECT
echo; echo "You are working from $PROJECT_TOPDIR and creating the project '$PROJECT_NAME'. Is that correct [y/n]?"; read input
if [ $input != "y" ]; then echo; exit; fi; echo

# SANITY CHECK FOR PATHS (FIX THIS!)
[ ! -d "$PROJECT_TOPDIR/$RAW_DATA_DIR" ] && echo && echo "$RAW_DATA_DIR is not valid as raw data directory" && echo && exit
[ ! -d "$PROJECT_DIR" ] && echo && echo "$PROJECT_DIR is not valid as project directory" && echo && exit
[ ! -d "$TEMPLATE_DIR" ] && echo && echo "$TEMPLATE_DIR is not valid as template directory" && echo && exit


# COPY THE NECESSARY TEMPLATE FILES
cd $TEMPLATE_DIR

# The run-script, S4, S5 and S6 are the same for all files and rename correctly
cp run_eureka.py S4*.ecf S5* S6*.ecf $PROJECT_TOPDIR/$PROJECT_DIR
cd $PROJECT_TOPDIR/$PROJECT_DIR
for file in *template.*; do mv "$file" "${file/template/$PROJECT_NAME}"; done
sed -i "s/eventlabel = 'event'/eventlabel = '$PROJECT_NAME'/" run_eureka.py

# Copy instrument-specific files
cd $TEMPLATE_DIR; 
if [ $INSTRUMENT="nirspec" ]; then 
	cp S1_nirx* *nirspec* $PROJECT_TOPDIR/$PROJECT_DIR 
	cd $PROJECT_TOPDIR/$PROJECT_DIR
	mv S1_nirx_template.ecf S1_$PROJECT_NAME.ecf
	mv S2_nirspec_fs_template.ecf S2_$PROJECT_NAME.ecf
	mv S3_nirspec_fs_template.ecf S3_$PROJECT_NAME.ecf
fi

# CHANGE PROJECT TOP DIRECTORY IN ALL FILES
cd $PROJECT_TOPDIR/$PROJECT_DIR; INDEX=0
for file in S*_$PROJECT_NAME.ecf; 
	do
	
	# Generate Input and Output Pointers
	IDX_IN=$INDEX
	IDX_OUT=$((INDEX + 1))
	
	# Replace strings
	sed -i "/topdir/ s%/home/User%$PROJECT_TOPDIR%" $file
	sed -i "/inputdir/ s%/.*%/$PROJECT_DIR/Stage$IDX_IN%" $file
	sed -i "/outputdir/ s%/.*%/$PROJECT_DIR/Stage$IDX_OUT%" $file
	
	# Need exception for Stage 1
	if [ $INDEX == 0 ]; then
	sed -i "/inputdir/ s%/.*%/$RAW_DATA_DIR%" $file
	fi
	
	# Additional change for LC fit stage
	if [ $INDEX == 4 ]; then
	LC_PARAM="S5_fit_par_$PROJECT_NAME.epf"
	sed -i "/fit_par/ s%S.*%$LC_PARAM%" $file
	fi
	
	INDEX=$IDX_OUT
	done

# WHOLE PROCESS IS DONE
echo "Created new Eureka! project"; echo

