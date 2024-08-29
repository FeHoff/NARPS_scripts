#!/bin/bash

##
#				EPI preprocessing using FSL, FIX, AROMA
#
#       Felix Hoffstaedter
#       Data and Platforms
#       INM-7 - Brain and Behaviour
#       September 2017 - Research Center Juelich
#
# 			use: execute plineEPI_RS [ <4DrawEPI> <subject_name> <outputdir> <TR> <fieldmapdir>]
##
EPIRAW=$1				  # raw data as Nifti
subject=$2				# subject name
OUT=$3						# output directoty
TR=$4							# TR
CAT=$5					# fieldmap folder if present

QC=$OUT/QC				# fixed QC folder: subject ohne "RS/name"
CWD=$(pwd)

AROMADIR=/data/BnB2/TOOLS/AROMA-0.3beta
FIXDIR=/data/BnB2/TOOLS/fix

QCnew=yes #yes				# change to "yes" to redo quality control images
NOcleanup=no #yes		# change to "yes" to keep intermediate files
JURECA=no #yes			# change to "yes" if on JURECA

if [ $JURECA = yes ]; then
	source /data/inm1/mapping/software_func.source
fi

# WhII_Standard.RData derived from more traditional early parallel scanning in the
# Whitehall imaging study, using no EPI acceleration: TR=3s, Resolution=3x3x3mm,
# Session=10mins, no spatial smoothing, 100s FWHM highpass temporal filtering.

# UKBiobank.RData derived from fairly HCP-like scanning in the UK Biobank imaging study:
# 40 subjects, TR=0.735s, Resolution=2.4x2.4x2.4mm, Session=6mins, no spatial smoothing,
# 100s FWHM highpass temporal filtering.

trainData=$FIXDIR/training_files/UKBiobank.RData
# trainData=$FIXDIR/training_files/WhII_MB6.RData

if [ ! -d $OUT ]; then	# create output folder if missing
	mkdir $OUT
fi

# functions:			prepro calls [ EPI AROMA FIX QC_EPI ]

function prepro	{

	# $1 = subject name
	# define directories & logfile
	CATDIR=$CAT
	# FUNCDIR=$OUT/$1/struc
	FUNCDIR=$OUT/
	LOCDIR=$FUNCDIR
	# LOCDIR=/dev/shm/$(basename $OUT)/$subject
	subLog=$OUT/${subject}_EPI.log

	# checks whether basic images for functional preprocessing are available
	if [ -f $EPIRAW ] && [ -f $CATDIR/mri/rmi${1}_affine.nii ]&& [ ! -f $FUNCDIR/AROMA_func_data_MNI152.nii ]; then
			date >> $subLog
			# if [ ! -f $FUNCDIR/AROMA_func_data_MNI152.nii.gz ] && [ ! -f $FUNCDIR/smoothed_func_data.nii.gz ]; then
			# 	echo "Remove incomplete func/ of subject $subject" >> $subLog
			# 	rm -R $FUNCDIR
			# fi
			echo "Start preprocessing of $EPIRAW of subject $subject" >> $subLog
			T1_CAT $subject
			# QC_T1_CAT $subject
			EPI $subject $EPIRAW
			CLN
			AROMA $subject
			CLN
			# FIX $subject
			# CLN
			# QC_EPI $subject
	else
		  echo "--- $EPIRAW or $subject/CAT results missing for subject $subject "
			echo "--- $EPIRAW or $subject/CAT results missing for subject $subject " >> $subLog
	fi

}

### log commands
run() {
  echo $@ >> $subLog
  $@
}


function T1_CAT	{

	if [ -f $FUNCDIR/cat4fsl/T1_MNI152_2mm.nii.gz ] && [ -f $FUNCDIR/cat4fsl/T1_to_MNI_warp.nii.gz ]; then
		date >> $subLog; echo "- FSL warp field from CATed T1 for $subject found"	>> $subLog
	else
		if [ -f $CATDIR/mri/rmi${1}_affine.nii ]; then
			date >> $subLog; echo "FSL warp field from CATed T1 for $subject"	>> $subLog
			run mkdir -p $FUNCDIR/cat4fsl
			run imcp $CATDIR/mri/rmi${1}_affine.nii $FUNCDIR/cat4fsl/T1.nii
			run imcp $CATDIR/mri/rp0${1}_affine.nii $FUNCDIR/cat4fsl/T1_pve.nii
			run imcp $CATDIR/mri/wp0${1}.nii $FUNCDIR/cat4fsl/pveseg_MNI152_1p5mm.nii

			run fslmaths $FUNCDIR/cat4fsl/T1.nii -mas $FUNCDIR/cat4fsl/T1_pve.nii $FUNCDIR/cat4fsl/T1_brain
			run flirt -in $FUNCDIR/cat4fsl/T1_brain -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain -omat $FUNCDIR/cat4fsl/T1_to_MNI_affine.mat
			run fnirt --in=$FUNCDIR/cat4fsl/T1 --ref=$FSLDIR/data/standard/MNI152_T1_2mm --fout=$FUNCDIR/cat4fsl/T1_to_MNI_warp --aff=$FUNCDIR/cat4fsl/T1_to_MNI_affine.mat --iout=$FUNCDIR/cat4fsl/T1_MNI152_2mm --logout=$FUNCDIR/cat4fsl/T1_to_MNI_nonlin.txt --config=$FSLDIR/etc/flirtsch/T1_2_MNI152_2mm.cnf
			run fslmaths $FUNCDIR/cat4fsl/T1_MNI152_2mm -mas $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil $FUNCDIR/cat4fsl/T1_MNI152_2mm_brain
		fi
	fi

}


function QC_T1_CAT	{

	if [ -f $FUNCDIR/cat4fsl/T1_MNI152_2mm.nii.gz ]; then
		if [ ! -f $QC/struct/$subject_T1cat-MNI152.png ] || [ $QCnew = yes ]; then
			mkdir -p $QC/struct/
			date >> $subLog; echo "QC: T1cat to MNI152"	>> $subLog
			# QC: T1 to MNI warp # creation of png images for quality control
			slicer $FUNCDIR/cat4fsl/T1_MNI152_2mm $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 $FUNCDIR/cat4fsl/sla.png -x 0.45 $FUNCDIR/cat4fsl/slb.png -x 0.55 $FUNCDIR/cat4fsl/slc.png -x 0.65 $FUNCDIR/cat4fsl/sld.png -y 0.35 $FUNCDIR/cat4fsl/sle.png -y 0.45 $FUNCDIR/cat4fsl/slf.png -y 0.55 $FUNCDIR/cat4fsl/slg.png -y 0.65 $FUNCDIR/cat4fsl/slh.png -z 0.35 $FUNCDIR/cat4fsl/sli.png -z 0.45 $FUNCDIR/cat4fsl/slj.png -z 0.55 $FUNCDIR/cat4fsl/slk.png -z 0.65 $FUNCDIR/cat4fsl/sll.png
			pngappend $FUNCDIR/cat4fsl/sla.png + $FUNCDIR/cat4fsl/slb.png + $FUNCDIR/cat4fsl/slc.png + $FUNCDIR/cat4fsl/sld.png + $FUNCDIR/cat4fsl/sle.png + $FUNCDIR/cat4fsl/slf.png + $FUNCDIR/cat4fsl/slg.png + $FUNCDIR/cat4fsl/slh.png + $FUNCDIR/cat4fsl/sli.png + $FUNCDIR/cat4fsl/slj.png + $FUNCDIR/cat4fsl/slk.png + $FUNCDIR/cat4fsl/sll.png $FUNCDIR/cat4fsl/T1_to_MNI_warped1.png
			slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain $FUNCDIR/cat4fsl/T1_MNI152_2mm -s 2 -x 0.35 $FUNCDIR/cat4fsl/sla.png -x 0.45 $FUNCDIR/cat4fsl/slb.png -x 0.55 $FUNCDIR/cat4fsl/slc.png -x 0.65 $FUNCDIR/cat4fsl/sld.png -y 0.35 $FUNCDIR/cat4fsl/sle.png -y 0.45 $FUNCDIR/cat4fsl/slf.png -y 0.55 $FUNCDIR/cat4fsl/slg.png -y 0.65 $FUNCDIR/cat4fsl/slh.png -z 0.35 $FUNCDIR/cat4fsl/sli.png -z 0.45 $FUNCDIR/cat4fsl/slj.png -z 0.55 $FUNCDIR/cat4fsl/slk.png -z 0.65 $FUNCDIR/cat4fsl/sll.png
			pngappend $FUNCDIR/cat4fsl/sla.png + $FUNCDIR/cat4fsl/slb.png + $FUNCDIR/cat4fsl/slc.png + $FUNCDIR/cat4fsl/sld.png + $FUNCDIR/cat4fsl/sle.png + $FUNCDIR/cat4fsl/slf.png + $FUNCDIR/cat4fsl/slg.png + $FUNCDIR/cat4fsl/slh.png + $FUNCDIR/cat4fsl/sli.png + $FUNCDIR/cat4fsl/slj.png + $FUNCDIR/cat4fsl/slk.png + $FUNCDIR/cat4fsl/sll.png $FUNCDIR/cat4fsl/T1_to_MNI_warped2.png
			pngappend $FUNCDIR/cat4fsl/T1_to_MNI_warped1.png - $FUNCDIR/cat4fsl/T1_to_MNI_warped2.png $QC/struct/$subject_T1cat-MNI152.png
			rm -f $FUNCDIR/cat4fsl/sl?.png $FUNCDIR/cat4fsl/T1_to_MNI_warped1.png $FUNCDIR/cat4fsl/T1_to_MNI_warped2.png
		fi
	fi

}


function EPI	{

  if [ -f $FUNCDIR/filtered_func_data_MNI152.nii.gz ] && [ -f $FUNCDIR/conf/meanWM.txt ]; then
    date >> $subLog; echo "- EPI preprocessing for $subject already performed"	>> $subLog
  else
    if [ -f $EPIRAW ]; then
		    if [ -f $FUNCDIR/reg/example_func.nii.gz ]; then
			    date >> $subLog; echo "	- EPIs ready for preprocessing"	>> $subLog
		    else
			    mkdir -p $FUNCDIR/reg
			    run fslmaths $EPIRAW $FUNCDIR/prefiltered_func_data -odt float
			    date >> $subLog; echo "	Initialising functional data"	>> $subLog
			    NUMVOL=$(fslnvols $FUNCDIR/prefiltered_func_data)
			    run fslroi $FUNCDIR/prefiltered_func_data $FUNCDIR/reg/example_func $(expr $NUMVOL / 2) 1
		    fi

				# EPI to T1
		    if [ -f $FUNCDIR/reg/example_func2standard.mat ] && [ -f $FUNCDIR/reg/example_func2standard_warp.nii.gz ]; then
			    date >> $subLog; echo "	- EPI registration to T1 already performed"	>> $subLog
		    else
			    date >> $subLog; echo "	EPI registration to T1"	>> $subLog
			    # registration EPI to T1, including Fieldmap based unwarping if applicable
					# if [ -d $FMDIR ]; then
						# # --> gradient distortion correction (if fieldmap applicable)
						# 			    run bet $FMDIR/FIELD_MAP_MAG_${FM_ID} $FUNCDIR/reg/FIELD_MAP_MAG_${FM_ID}_brain
						# 			    run fslmaths $FUNCDIR/reg/FIELD_MAP_MAG_${FM_ID}_brain -ero $FUNCDIR/reg/FIELD_MAP_MAG_${FM_ID}_brain
						# 			    run fsl_prepare_fieldmap SIEMENS $FMDIR/FIELD_MAP_PHASE_${FM_ID} $FUNCDIR/reg/FIELD_MAP_MAG_${FM_ID}_brain $FUNCDIR/reg/FIELD_MAP_RADS_${FM_ID} 2.46
						# 			    run epi_reg --epi=$FUNCDIR/reg/example_func --t1=$FUNCDIR/cat4fsl/T1.nii --t1brain=$FUNCDIR/cat4fsl/T1_brain --out=$FUNCDIR/reg/example_func2highres --fmap=$FUNCDIR/reg/FIELD_MAP_RADS_${FM_ID} --fmapmag=$FMDIR/FIELD_MAP_MAG_${FM_ID} --fmapmagbrain=$FUNCDIR/reg/FIELD_MAP_MAG_${FM_ID}_brain --echospacing=0.00069 --pedir=-y
					# else
						# --> registration of middle reference EPI_image onto T1 (needs whole head image and brainextracted image), T1 = biasfield_corrected (maybe rather ---wm_segm than t1_brain (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT/UserGuide)); Output is 1) an image, multiplication with it registers EPIs to T1 (subjectspace) 2) .mat file
					  run epi_reg --epi=$FUNCDIR/reg/example_func --t1=$FUNCDIR/cat4fsl/T1.nii --t1brain=$FUNCDIR/cat4fsl/T1_brain --out=$FUNCDIR/reg/example_func2highres
			  	# fi

			    run convert_xfm -inverse -omat $FUNCDIR/reg/highres2example_func.mat $FUNCDIR/reg/example_func2highres.mat
				 	run convert_xfm -omat $FUNCDIR/reg/example_func2standard.mat -concat $FUNCDIR/cat4fsl/T1_to_MNI_affine.mat $FUNCDIR/reg/example_func2highres.mat
		     	run convertwarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain --premat=$FUNCDIR/reg/example_func2highres.mat --warp1=$FUNCDIR/cat4fsl/T1_to_MNI_warp --out=$FUNCDIR/reg/example_func2standard_warp
			    run applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -i $FUNCDIR/reg/example_func -o $FUNCDIR/reg/example_func2standard -w $FUNCDIR/reg/example_func2standard_warp
					run convert_xfm -inverse -omat $FUNCDIR/reg/standard2func.mat $FUNCDIR/reg/example_func2standard.mat
		    fi

				# # slice timing correction
		    # if [ ! -f $FUNCDIR/prefiltered_func_data_fs.nii.gz ] && [ ! -f $FUNCDIR/mc/prefiltered_func_data_mcf.par ]; then
			  #   date >> $subLog; echo "	Slice time correction"	>> $subLog
			  #   # slicetimer -i $FUNCDIR/mc/prefiltered_func_data_mcf --out=$FUNCDIR/prefiltered_func_data_st -r $TR
			  #   ### Filter-Shift - new slicetimer method (Parker et al. 2017)
			  #   run $(pwd)/filtershift -i $FUNCDIR/prefiltered_func_data -o $FUNCDIR/prefiltered_func_data_fs --TR=$TR
		    # fi

				###					with or without Slice timing ?!?
				if [ ! -f $FUNCDIR/mc/prefiltered_func_data_mcf.par ]; then
			    date >> $subLog; echo "	Motion correction"	>> $subLog
			    mkdir -p $FUNCDIR/mc
			    # motion correction ### MAYBE without slice time correction
					run mcflirt -in $FUNCDIR/prefiltered_func_data -out $FUNCDIR/mc/prefiltered_func_data_mcf -mats -plots -reffile $FUNCDIR/reg/example_func -rmsrel -rmsabs -spline_final
					#	RUN WITH filtershift
			    # run mcflirt -in $FUNCDIR/prefiltered_func_data_fs -out $FUNCDIR/mc/prefiltered_func_data_mcf -mats -plots -reffile $FUNCDIR/reg/example_func -rmsrel -rmsabs -spline_final
		    fi

		    if [ ! -f $FUNCDIR/filtered_func_data.nii.gz ] ; then
					date >> $subLog; echo "	Intensity normalization, highpass filtering & mean creation"	>> $subLog
			    # create mean functional image	- for FIX
			    run fslmaths $FUNCDIR/mc/prefiltered_func_data_mcf -Tmean $FUNCDIR/mean_func
			    # create brain mask of 4D data	- for AROMA & FIX
			    run bet $FUNCDIR/mean_func $FUNCDIR/mask -f 0.3 -n -m
			    run immv $FUNCDIR/mask_mask $FUNCDIR/mask
					# intensity normalization
	    		MEDINTENSITY=$(fslstats $FUNCDIR/mc/prefiltered_func_data_mcf -k $FUNCDIR/mask -p 50)
			    BRTHRESH=$(echo $MEDINTENSITY \* 0.75 | bc -l)
			    INTTHRESH=$(echo 10000 / $MEDINTENSITY | bc -l)
			    run fslmaths $FUNCDIR/mc/prefiltered_func_data_mcf -mas $FUNCDIR/mask $FUNCDIR/mc/prefiltered_func_data_mcf
			    run fslmaths $FUNCDIR/mc/prefiltered_func_data_mcf -mul $INTTHRESH $FUNCDIR/prefiltered_func_data_intnorm
					run imcp $FUNCDIR/prefiltered_func_data_intnorm $FUNCDIR/filtered_func_data
			    run fslmaths $FUNCDIR/filtered_func_data -Tmean $FUNCDIR/mean_func
					# smoothing for AROMA
					run susan $FUNCDIR/prefiltered_func_data_intnorm $BRTHRESH 2.12314225053 3 1 1 $FUNCDIR/mean_func $BRTHRESH $FUNCDIR/smoothed_func_data
				fi

				if [ ! -f $FUNCDIR/filtered_func_data_MNI152.nii.gz ]; then
			   	date >> $subLog; echo "	EPI to MNI152 warping "	>> $subLog
					FSLOUTPUTTYPE=NIFTI
			    run applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -i $FUNCDIR/prefiltered_func_data_intnorm -o $FUNCDIR/filtered_func_data_MNI152 -w $FUNCDIR/reg/example_func2standard_warp
					fslmaths $FUNCDIR/filtered_func_data_MNI152 -Tmean $FUNCDIR/filtered_mean_func_data_MNI152
					FSLOUTPUTTYPE=NIFTI_GZ
		    fi

				# create GM, WM & CSF masks
		    if [ ! -f $FUNCDIR/conf/WM.nii.gz ]; then
			    mkdir -p $FUNCDIR/conf
					run flirt -in $FUNCDIR/cat4fsl/pveseg_MNI152_1p5mm -ref $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz -usesqform -applyxfm -out $FUNCDIR/cat4fsl/T1_pve_MNI152_2mm
			    run fslmaths $FUNCDIR/cat4fsl/T1_pve_MNI152_2mm -thr 0.9 -uthr 1.1 -bin $FUNCDIR/conf/CSF;
			    run fslmaths $FUNCDIR/cat4fsl/T1_pve_MNI152_2mm -thr 1.9 -uthr 2.1 -bin $FUNCDIR/conf/GM;
			    run fslmaths $FUNCDIR/cat4fsl/T1_pve_MNI152_2mm -thr 2.9 -bin $FUNCDIR/conf/WM;
		    fi

				# calculates mean GS, GM, WM & CSF (for regression)
		    if [ -f $FUNCDIR/filtered_func_data_MNI152.nii.gz ] && [ ! -f $FUNCDIR/conf/meanWM.txt ]; then
			    run fslmeants -i $FUNCDIR/filtered_func_data_MNI152 -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask -o $FUNCDIR/conf/meanGS.txt;
			    run fslmeants -i $FUNCDIR/filtered_func_data_MNI152 -m $FUNCDIR/conf/CSF -o $FUNCDIR/conf/meanCSF.txt;
			    run fslmeants -i $FUNCDIR/filtered_func_data_MNI152 -m $FUNCDIR/conf/GM -o $FUNCDIR/conf/meanGM.txt;
			    run fslmeants -i $FUNCDIR/filtered_func_data_MNI152 -m $FUNCDIR/conf/WM -o $FUNCDIR/conf/meanWM.txt;
			    date >> $subLog; echo "	# EPI preprocessing finished"	>> $subLog
		    fi
    else
	    date >> $subLog; echo "	# EPI data missing"	>> $subLog
	    return
    fi
  fi

}


function AROMA	{

	if [ -f $FUNCDIR/AROMA_func_data_MNI152.nii.gz ] && [ -f $FUNCDIR/conf/AROMA_meanWM.txt ]; then
		date >> $subLog; echo "AROMA already performed"	>> $subLog
	else
		date >> $subLog; echo "start AROMA"	>> $subLog
		if [ -f $FUNCDIR/smoothed_func_data.nii.gz ]; then
			if [ ! -f $FUNCDIR/AROMA/denoised_func_data_nonaggr.nii.gz ] && [ ! -f $FUNCDIR/AROMA.zip ]; then
			 	date >> $subLog; echo " AROMA EPI cleaning of $subject"	>> $subLog
				# AROMA: ICA-based Automatic Removal Of Motion Artifacts
				# -- files needed for AROMA --
				# smoothed_func_data.nii.gz			smoothed but NOT highpass filtered 4D data
				# example_func2highres.mat      FLIRT transform from functional to structural space
				# T1_to_MNI_nonlin_field.nii.gz       warp-file of the structural data to MNI152 space
				# mc/prefiltered_func_data_mcf.par    motion parameters created by mcflirt (in mc subdirectory)
				# run python /data/inm1/mapping/software/2016b/installed/AROMA-0.3beta/ICA_AROMA.py -i $FUNCDIR/smoothed_func_data.nii.gz -o $FUNCDIR/AROMA -a $FUNCDIR/reg/example_func2highres.mat -w $FUNCDIR/cat4fsl/T1_to_MNI_warp.nii.gz -mc $FUNCDIR/mc/prefiltered_func_data_mcf.par -m $FUNCDIR/mask.nii.gz -tr $TR
				run python $AROMADIR/ICA_AROMA.py -i $FUNCDIR/smoothed_func_data.nii.gz -o $LOCDIR/AROMA -a $FUNCDIR/reg/example_func2highres.mat -w $FUNCDIR/cat4fsl/T1_to_MNI_warp.nii.gz -mc $FUNCDIR/mc/prefiltered_func_data_mcf.par -m $FUNCDIR/mask.nii.gz -tr $TR
				run imcp $LOCDIR/AROMA/denoised_func_data_nonaggr $FUNCDIR/AROMA_func_data
				date >> $subLog; echo "	AROMA_EPI to MNI152 warping & MELODIC zipping "	>> $subLog
				FSLOUTPUTTYPE=NIFTI
				run applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -i $LOCDIR/AROMA/denoised_func_data_nonaggr -o $FUNCDIR/AROMA_func_data_MNI152 -w $FUNCDIR/reg/example_func2standard_warp
				fslmaths $FUNCDIR/AROMA_func_data_MNI152 -Tmean $FUNCDIR/AROMA_mean_func_data_MNI152
				FSLOUTPUTTYPE=NIFTI_GZ
				# ZIP melodic output as it kills the file quota
				cd $LOCDIR
				run zip -rmq0 $FUNCDIR/AROMA.zip AROMA
				run rm -R AROMA
				# run rm -R $LOCDIR
				cd $CWD
			fi

			# calculates mean GS, GM, WM & CSF for noise regression
			if [ -f $FUNCDIR/AROMA_func_data_MNI152.nii.gz ] && [ ! -f $FUNCDIR/conf/AROMA_meanWM.txt ]; then
				run fslmeants -i $FUNCDIR/AROMA_func_data_MNI152 -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask -o $FUNCDIR/conf/AROMA_meanGS.txt;
				run fslmeants -i $FUNCDIR/AROMA_func_data_MNI152 -m $FUNCDIR/conf/CSF -o $FUNCDIR/conf/AROMA_meanCSF.txt;
				run fslmeants -i $FUNCDIR/AROMA_func_data_MNI152 -m $FUNCDIR/conf/GM -o $FUNCDIR/conf/AROMA_meanGM.txt;
				run fslmeants -i $FUNCDIR/AROMA_func_data_MNI152 -m $FUNCDIR/conf/WM -o $FUNCDIR/conf/AROMA_meanWM.txt;
				date >> $subLog; echo "AROMA finished"	>> $subLog
			fi
		else
			date >> $subLog; echo "# EPIs not prepared fo AROMA #"	>> $subLog
		fi
	fi

}


function FIX	{

	if [ -f $FUNCDIR/filtered_func_data_clean.nii.gz ] && [ -f $FUNCDIR/conf/FIX_meanWM.txt ]; then
		date >> $subLog; echo "FIX already performed"	>> $subLog
	else
		if [ -f $FUNCDIR/filtered_func_data.nii.gz ] && [ ! -f $FUNCDIR/filtered_func_data_clean.nii.gz ]; then
			# if [ -d $FUNCDIR/filtered_func_data.ica ]; then
			# 	# run mkdir $LOCDIR
			# 	run rm -r $FUNCDIR/filtered_func_data.ica
			# fi
			# 	date >> $subLog; echo "start FIX"	>> $subLog
			# 	# prepare DATA for FIX / copy files to local DIR (preferable to RAMdisk)
			# 	# run mkdir $LOCDIR
			# 	# run cp -r $FUNCDIR/reg $LOCDIR/
			# 	# run cp -r $FUNCDIR/mc $LOCDIR/
			# 	# run imcp $FUNCDIR/cat4fsl/T1_brain $LOCDIR/reg/highres
			# 	# run imcp $FUNCDIR/cat4fsl/T1_pve $LOCDIR/reg/highres_pveseg
			# 	# run imcp $FUNCDIR/mask $LOCDIR/mask
			# 	# run imcp $FUNCDIR/mean_func $LOCDIR/mean_func
			# 	# run imcp $FUNCDIR/filtered_func_data $LOCDIR/filtered_func_data
			# 	# echo "	Melodic for FIX computed locally"	>> $subLog
			# 	run melodic -i $LOCDIR/filtered_func_data -o $LOCDIR/filtered_func_data.ica -v --nobet --bgthreshold=1 --tr=$TR -d 0 --mmthresh=0.5 --report

			if [ -f $FUNCDIR/filtered_func_data.ica.zip ] && [ ! -d $FUNCDIR/filtered_func_data.ica ]; then
				# run mkdir $LOCDIR
				run unzip $FUNCDIR/filtered_func_data.ica.zip -d $LOCDIR
			elif [ ! -d $FUNCDIR/filtered_func_data.ica ]; then
				date >> $subLog; echo "start FIX"	>> $subLog
				# prepare DATA for FIX / copy files to local DIR (preferable to RAMdisk)
				# run mkdir $LOCDIR
				# run cp -r $FUNCDIR/reg $LOCDIR/
				# run cp -r $FUNCDIR/mc $LOCDIR/
				run imcp $FUNCDIR/cat4fsl/T1_brain $LOCDIR/reg/highres
				run imcp $FUNCDIR/cat4fsl/T1_pve $LOCDIR/reg/highres_pveseg
				# run imcp $FUNCDIR/mask $LOCDIR/mask
				# run imcp $FUNCDIR/mean_func $LOCDIR/mean_func
				# run imcp $FUNCDIR/filtered_func_data $LOCDIR/filtered_func_data
				echo "	Melodic for FIX "	>> $subLog
				run melodic -i $LOCDIR/filtered_func_data -o $LOCDIR/filtered_func_data.ica -v --nobet --bgthreshold=1 --tr=$TR -d 0 --mmthresh=0.5 --report
			fi
			### FIX - FMRIB's ICA-based Xnoiseifier
			# -- files needed for FIX --
			# filtered_func_data.nii.gz          preprocessed 4D data
		  # filtered_func_data.ica             melodic (command-line program) full output directory
		  # mc/prefiltered_func_data_mcf.par   motion parameters created by mcflirt (in mc subdirectory)
		  # mask.nii.gz                        valid mask relating to the 4D data
		  # mean_func.nii.gz                   temporal mean of 4D data
		  # reg/example_func.nii.gz            example image from 4D data
		  # reg/highres.nii.gz                 brain-extracted structural
		  # reg/highres2example_func.mat       FLIRT transform from structural to functional space
		  ##### design.fsf                         FEAT/MELODIC setup file; if present, this controls the
		  #                                      default temporal filtering of motion parameters

			echo "FIX-ing EPIs of $subject"	>> $subLog
			run $FIXDIR/fix $LOCDIR $trainData 20 -m
			### without motion regression if needed
			# fix $LOCDIR $trainData 20

			# # ZIP FIX stuff
			cd $LOCDIR
			run zip -rmq0 $FUNCDIR/filtered_func_data.ica.zip filtered_func_data.ica
			run rm -R filtered_func_data.ica
			cd $CWD
			# run mkdir $FUNCDIR/fix
			# train=$(basename $trainData)
			# run cp $LOCDIR/fix4melview_${train%%.*}_thr20.txt $FUNCDIR/fix/fix4melview_${train%%.*}_thr20.txt
			# run cp $LOCDIR/fix/features.mat $FUNCDIR/fix/features.mat
			# run cp $LOCDIR/fix/features.csv $FUNCDIR/fix/features.csv

			date >> $subLog; echo "	FIX_EPI to MNI152 warping & write mean time series "	>> $subLog
			run imcp $LOCDIR/filtered_func_data_clean $FUNCDIR/FIX_func_data
			FSLOUTPUTTYPE=NIFTI
			run applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -i $LOCDIR/filtered_func_data_clean -o $FUNCDIR/FIX_func_data_MNI152 -w $LOCDIR/reg/example_func2standard_warp
			FSLOUTPUTTYPE=NIFTI_GZ
			### careful NOT to remove FUNCDIR !!!
			# run rm -R $LOCDIR

			# calculates mean GS, GM, WM & CSF (for noise regression)
			run fslmeants -i $FUNCDIR/FIX_func_data_MNI152 -m $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask -o $FUNCDIR/conf/FIX_meanGS.txt;
			run fslmeants -i $FUNCDIR/FIX_func_data_MNI152 -m $FUNCDIR/conf/CSF -o $FUNCDIR/conf/FIX_meanCSF.txt;
			run fslmeants -i $FUNCDIR/FIX_func_data_MNI152 -m $FUNCDIR/conf/GM -o $FUNCDIR/conf/FIX_meanGM.txt;
			run fslmeants -i $FUNCDIR/FIX_func_data_MNI152 -m $FUNCDIR/conf/WM -o $FUNCDIR/conf/FIX_meanWM.txt;
			date >> $subLog; echo "FIX finished"	>> $subLog
		else
			date >> $subLog; echo "# EPIs not prepared fo FIX #"	>> $subLog
		fi
	fi

}

function CLN	{

  # cleanup obsolete files
  if [ ! $NOcleanup = yes ]; then
	 	if [ -f $FUNCDIR/filtered_func_data_MNI152.nii.gz ]; then
    	run imrm $FUNCDIR/prefiltered_func_data $FUNCDIR/prefiltered_func_data_fs $FUNCDIR/mc/prefiltered_func_data_mcf $FUNCDIR/smoothed_func_data_usan_size
		fi

		if [ -f $FUNCDIR/AROMA_func_data_MNI152.nii.gz ]; then
    	run imrm $FUNCDIR/smoothed_func_data
		fi

		# if [ -f $FUNCDIR/FIX_func_data_MNI152.nii.gz ]; then
		# 	run imrm $FUNCDIR/filtered_func_data
		# fi
	fi

}

function QC_EPI	{

	if [ ! -f $QC/func/$subject_EPI-MNI152_warp.png ] || [ $QCnew = yes ]; then
		if [ -f $FUNCDIR/reg/example_func2highres.nii.gz ]; then
			mkdir -p $QC/func
			date >> $subLog; echo "	QC: T1 EPI registration to T1"	>> $subLog
			# QC: EPI to MNI warp # creation of png images for quality control
			slicer $FUNCDIR/reg/example_func2highres $FUNCDIR/cat4fsl/T1.nii -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/example_func2highres1.png
			slicer $FUNCDIR/cat4fsl/T1.nii $FUNCDIR/reg/example_func2highres -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/example_func2highres2.png
			pngappend $FUNCDIR/example_func2highres1.png - $FUNCDIR/example_func2highres2.png $QC/func/$subject_EPI-T1.png
			rm -f $FUNCDIR/sl?.png $FUNCDIR/example_func2highres1.png $FUNCDIR/example_func2highres2.png
			# QC: movement plots for quality check
			fsl_tsplot -i $FUNCDIR/mc/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o $FUNCDIR/rot.png
			fsl_tsplot -i $FUNCDIR/mc/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o $FUNCDIR/trans.png
			fsl_tsplot -i $FUNCDIR/mc/prefiltered_func_data_mcf_abs.rms,$FUNCDIR/mc/prefiltered_func_data_mcf_rel.rms -t 'MCFLIRT estimated mean displacement (mm)' -u 1 -w 640 -h 144 -a absolute,relative -o $FUNCDIR/disp.png
			pngappend $FUNCDIR/rot.png - $FUNCDIR/trans.png - $FUNCDIR/disp.png $QC/func/$subject_EPI_motion.png
			rm -f $FUNCDIR/rot.png $FUNCDIR/trans.png $FUNCDIR/disp.png

			# QC: EPI to MNI warp # creation of png images for quality control
			slicer $FUNCDIR/filtered_func_data_MNI152 $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped1.png
			slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain $FUNCDIR/filtered_func_data_MNI152 -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped2.png
			pngappend $FUNCDIR/EPI_to_MNI_warped1.png - $FUNCDIR/EPI_to_MNI_warped2.png $QC/func/$subject_EPI-MNI152_warp.png
			rm -f $FUNCDIR/sl?.png $FUNCDIR/EPI_to_MNI_warped1.png $FUNCDIR/EPI_to_MNI_warped2.png
		fi
	fi

	if [ ! -f $QC/func/$subject_AROMA-MNI152.png ] || [ $QCnew = yes ]; then
		if [ -f $FUNCDIR/AROMA_func_data_MNI152.nii.gz ]; then
			# QC: AROMA to MNI warp # creation of png images for quality control
			slicer $FUNCDIR/AROMA_func_data_MNI152 $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped1.png
			slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain $FUNCDIR/AROMA_func_data_MNI152 -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped2.png
			pngappend $FUNCDIR/EPI_to_MNI_warped1.png - $FUNCDIR/EPI_to_MNI_warped2.png $QC/func/$subject_AROMA-MNI152.png
			rm -f $FUNCDIR/sl?.png $FUNCDIR/EPI_to_MNI_warped1.png $FUNCDIR/EPI_to_MNI_warped2.png
		fi
	fi

	if [ ! -f $QC/func/$subject_FIX-MNI_warp.png ] || [ $QCnew = yes ]; then
		if [ -f $FUNCDIR/FIX_func_data_MNI.nii.gz ]; then
			# QC: FIX to MNI warp # creation of png images for quality control
			slicer $FUNCDIR/FIX_func_data_MNI $FSLDIR/data/standard/MNI152_T1_2mm_brain -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped1.png
			slicer $FSLDIR/data/standard/MNI152_T1_2mm_brain $FUNCDIR/FIX_func_data_MNI -s 2 -x 0.35 $FUNCDIR/sla.png -x 0.45 $FUNCDIR/slb.png -x 0.55 $FUNCDIR/slc.png -x 0.65 $FUNCDIR/sld.png -y 0.35 $FUNCDIR/sle.png -y 0.45 $FUNCDIR/slf.png -y 0.55 $FUNCDIR/slg.png -y 0.65 $FUNCDIR/slh.png -z 0.35 $FUNCDIR/sli.png -z 0.45 $FUNCDIR/slj.png -z 0.55 $FUNCDIR/slk.png -z 0.65 $FUNCDIR/sll.png
			pngappend $FUNCDIR/sla.png + $FUNCDIR/slb.png + $FUNCDIR/slc.png + $FUNCDIR/sld.png + $FUNCDIR/sle.png + $FUNCDIR/slf.png + $FUNCDIR/slg.png + $FUNCDIR/slh.png + $FUNCDIR/sli.png + $FUNCDIR/slj.png + $FUNCDIR/slk.png + $FUNCDIR/sll.png $FUNCDIR/EPI_to_MNI_warped2.png
			pngappend $FUNCDIR/EPI_to_MNI_warped1.png - $FUNCDIR/EPI_to_MNI_warped2.png $QC/func/$subject_FIX-MNI_warp.png
			rm -f $FUNCDIR/sl?.png $FUNCDIR/EPI_to_MNI_warped1.png $FUNCDIR/EPI_to_MNI_warped2.png
		fi
	fi

}

date
echo "Start EPI preprocessing of: $subject"
for subject in $subject;	do
			prepro $subject
done
