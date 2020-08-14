# SPANDAK+
-----------------
The full SPANDAK+ pipeline can be found in the SPANDAK+ directory.

![Flow Graph](spandak+.png)

The pipeline as a whole is broken into two parts.

#1. Extract2cdd_auto.py
-----------------

This first step of the pipeline requires three arguments:
[1] path to database containing .fil paths and corresponding .raw paths (currently accepts csv)
[2] path to csv for a single .fil from SPANDAK
[3] desired grade given by SPANDAK (choose 'B' or 'C')

First, raw voltages are extracted and place in directories specifying the start and end time of extraction, as well as the band across which they were extracted. This is done using the extract_blocks.py script written by Greg Hell (https://github.com/greghell/extractor). These directories are placed in a surrounding directory labeled 'bursts'. Next, the raw voltages are spliced together using the splicer_raw.py script, also written by Greg Hell (https://github.com/greghell/extractor). Finally, the spliced raw files are coherently dedispersed using DSPSR and place in 'fits' directories within each 'time_band' directory.

Example Command: python extract2cdd_auto.py /datax/scratch/jfaber/SPANDAK_extension/database.csv /datax/scratch/jfaber/SPANDAK_extension/57991_49905_DIAG_FRB121102_0011.csv B

#2. Polfluxrm_auto.py
-----------------

Once step 1 is complete, the second step can be initiated, and requires two arguments:
[1] path to directory containing fits files
[2] path to directory containing calibration files for the observation

Currently, step two is only run on individaul bursts—--this will be generalized in the future. The fits files are first converted to and stored as numpy arrays and pngs in case the pulse needs to be visually verified. The, using the calibration files, polarization and flux calibration are performed, followed by fitting of the rotation measure with the script RMfit_curve.py written by Vishal Gajjar. The RM fit results are output as a csv file '[pulse_id].calib.rmfit.csv', and the PA value, including error is printed in the terminal.

Example Command: python polfluxrm_auto.py /datax/scratch/jfaber/SPANDAK_extension/bursts/26.2_26.6_3.8_9/fits /datax/scratch/jfaber/SPANDAK_extension/calib_files
