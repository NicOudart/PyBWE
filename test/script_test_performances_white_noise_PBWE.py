################################################################################
#HEADER

#This script allows you to test the performances of the PBWE package facing
#different levels of white-noise, to check for regressions.
#It will generate a report containing the test results.

#The test consists in measuring the performances of the PBWE on a radar signal
#containing 2 echoes from 2 targets, separated by distances below
#or close to the radar's resolution.

#The following test parameters are user-defined:
#   -The list of distances between targets to be tested (m).
#   -The list of SNR to be tested (dB).
#   -The number of noise cases

#The following tests will be performed:
#   -The percentage of times 2 echoes are detected (with a user-defined
#    threshold).
#   -The error on the distance between targets, estimated from the surface echo.
#   -The error on the each echo's amplitude.

#For the test to be reproducible, the same random seeds are always used.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3
#    GHz measures a 1001 frequencies spectrum when sounding.
#   -This radar can transmit and receive signals in two linear polarizations
#    named 0 and 1. This leads to 4 polirization channels: 2 co-polar named 00
#    and 11 (emission and reception in the same polarization), and 2 cross-polar
#    named 01 and 10 (emission and reception in different polarizations). For
#    simplicity, we will only consider the co-polar channels 00 and 11 in this
#    example.
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. These targets generate echoes of given complex
#    amplitudes in the radar's signal, or complex sine-waves in the measured
#    spectrum.
#   -The measured spectrum is corrupted by a white-noise of standard deviation
#    10X smaller than the complex sine-waves' amplitudes.

#The parameters used for the PBWE are the default ones.

#References: Oudart et al. (2021)

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
import pandas as pd
import os
from math import pi
from scipy.signal import hilbert,find_peaks
from importlib_metadata import version

import PyPBWE

#Test parameters:--------------------------------------------------

#List of distances between targets to be tested (m):
list_dist_targets = [0.04,0.06,0.08,0.1,0.12]

#List of SNR (dB):
list_snr_levels = [6,14,20,40,60]

#Number of noise case to be tested for each distance:
nb_noise_case = 25

#Scenario parameters:-----------------------------------------------------------

#Amplitudes of the echoes corresponding to each target for each polarimetric
#channel (00 and 11):
amp_target1_00 = 1
amp_target2_00 = 1
amp_target1_11 = 1
amp_target2_11 = -1

#Distance (m) between the 1st target and the radar:
dist_target1 = 1

#Peak detection threshold on the amplitude of echoes:
detection_level = 0.5

#Retrieve the test path:--------------------------------------------------------

#Test directory path:
test_dir_path = os.path.dirname(__file__)

#Test report path:
test_report_path = os.path.join(test_dir_path,'PyPBWE_Report_test_performances_white_noise.md')

#Initialize the test results dataframes:----------------------------------------

#For the PBWE - polar 00:
dataframe_pbwe_polar00_percentage_echoes_detection = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_distance_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_distance_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_amplitude_1_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_amplitude_1_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_amplitude_2_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar00_amplitude_2_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)

#For the PBWE - polar 11:
dataframe_pbwe_polar11_percentage_echoes_detection = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_distance_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_distance_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_amplitude_1_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_amplitude_1_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_amplitude_2_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_pbwe_polar11_amplitude_2_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)

#Perform the test:--------------------------------------------------------------

#Initialize the iterations counter:
nb_iterations = len(list_dist_targets)*len(list_snr_levels)*nb_noise_case
count_iterations = 0

#Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
freq_vect = np.linspace(0.5e9,3e9,1001)

#Iterate on the distance between targets:
for dist in list_dist_targets:

    #Distance (m) between the 2nd radar target in free-space (returning one
    #echo) and the radar:
    dist_target2 = 1 + dist

    #Generate a sum of two complex sine-waves corresponding to the targets' echoes:
    spec_vect_00 = (amp_target1_00*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_00*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))
    spec_vect_11 = (amp_target1_11*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_11*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))

    #Only keep the real part of the spectrum (In-phase component):
    spec_vect_00 = np.real(spec_vect_00)
    spec_vect_11 = np.real(spec_vect_11)

    #Calculate the power of the signal (dB):
    power_sig_00 = 10*np.log10(np.mean(spec_vect_00**2))
    power_sig_11 = 10*np.log10(np.mean(spec_vect_11**2))

    #Iterate on the SNR values:
    for snr in list_snr_levels:

        #Initialize the 2 echoes detection counter:
        nb_echoes_detection_pbwe00 = 0
        nb_echoes_detection_pbwe11 = 0

        #Initialize the lists to stock all noise cases results:
        list_distance_pbwe00 = []
        list_amplitude_target1_pbwe00 = []
        list_amplitude_target2_pbwe00 = []
        list_distance_pbwe11 = []
        list_amplitude_target1_pbwe11 = []
        list_amplitude_target2_pbwe11 = []

        #Calculate the power of the noise (dB) to obtain the right SNR:
        power_noise_00 = power_sig_00-snr
        power_noise_11 = power_sig_11-snr

        #Iterate on the noise case:
        for idx_case in range(nb_noise_case):

            #Display the test progress:
            count_iterations += 1
            print('PBWE performance test in progress ... '+str(round(100*count_iterations/nb_iterations,3))+'%')

            #Add a white noise to this spectrum signal to obtain the right SNR:
            rng_00 = np.random.RandomState(idx_case)
            rng_11 = np.random.RandomState(nb_noise_case+idx_case+1)
            wn_vect_00 = rng_00.normal(0,np.sqrt(10**(power_noise_00/10)),spec_vect_00.shape)
            wn_vect_11 = rng_11.normal(0,np.sqrt(10**(power_noise_11/10)),spec_vect_11.shape)
            spec_vect_wn_00 = spec_vect_00 + wn_vect_00
            spec_vect_wn_11 = spec_vect_11 + wn_vect_11

            #Reconstruct a complex signal with the Hilbert transform:
            spec_vect_wn_00 = np.conjugate(hilbert(spec_vect_wn_00))[::2]
            spec_vect_wn_11 = np.conjugate(hilbert(spec_vect_wn_11))[::2]

            #Assemble the 2 spectrum channels into a single matrix:
            spec_mat_wn = np.vstack((spec_vect_wn_00,spec_vect_wn_11))

            #Application of the PBWE:
            output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_wn,df=5e6,extra_factor=3,model_order=0.33,zp_factor=10,side_cut=True)

            #Testing the PBWE - polar 00:---------------------------------------

            #Detection of the 2 echoes:
            peaks_pbwe_00 = find_peaks(abs(output_pbwe[0,:]),height=detection_level)[0]

            #Check that 2 echoes are detected:
            if len(peaks_pbwe_00)==2:

                #Increment by 1 the number of echoes detection:
                nb_echoes_detection_pbwe00 += 1

                #Stock the distance between echoes:
                list_distance_pbwe00 += [0.5*abs(time_pbwe_vect[peaks_pbwe_00[0]]-time_pbwe_vect[peaks_pbwe_00[1]])*3e8]

                #Stock the amplitude of the 1st echo:
                list_amplitude_target1_pbwe00 += [abs(output_pbwe[0,peaks_pbwe_00[0]])]

                #Stock the amplitude of the 2nd echo:
                list_amplitude_target2_pbwe00 += [abs(output_pbwe[0,peaks_pbwe_00[1]])]

            #Testing the PBWE - polar 11:---------------------------------------

            #Detection of the 2 echoes:
            peaks_pbwe_11 = find_peaks(abs(output_pbwe[1,:]),height=detection_level)[0]

            #Check that 2 echoes are detected:
            if len(peaks_pbwe_11)==2:

                #Increment by 1 the number of echoes detection:
                nb_echoes_detection_pbwe11 += 1

                #Stock the distance between echoes:
                list_distance_pbwe11 += [0.5*abs(time_pbwe_vect[peaks_pbwe_11[0]]-time_pbwe_vect[peaks_pbwe_11[1]])*3e8]

                #Stock the amplitude of the 1st echo:
                list_amplitude_target1_pbwe11 += [abs(output_pbwe[1,peaks_pbwe_11[0]])]

                #Stock the amplitude of the 2nd echo:
                list_amplitude_target2_pbwe11 += [abs(output_pbwe[1,peaks_pbwe_11[1]])]

        #Tests calculations:----------------------------------------------------

        #Calculate the percentage of echoes detection:
        pc_echoes_detection_pbwe00 = 100*nb_echoes_detection_pbwe00/nb_noise_case
        pc_echoes_detection_pbwe11 = 100*nb_echoes_detection_pbwe11/nb_noise_case

        #Calculate the error on the distance between echoes:
        error_distance_pbwe00 = np.array(list_distance_pbwe00)-dist
        error_distance_pbwe11 = np.array(list_distance_pbwe11)-dist

        #Calculate the error on the amplitude of the 1st echo:
        error_amplitude_target1_pbwe00 = 100*(np.array(list_amplitude_target1_pbwe00)-amp_target1_00)/amp_target1_00
        error_amplitude_target1_pbwe11 = 100*(np.array(list_amplitude_target1_pbwe11)-amp_target1_11)/amp_target1_11

        #Calculate the error on the amplitude of the 2nd echo:
        error_amplitude_target2_pbwe00 = 100*(np.array(list_amplitude_target2_pbwe00)-amp_target2_00)/amp_target2_00
        error_amplitude_target2_pbwe11 = 100*(np.array(list_amplitude_target2_pbwe11)-amp_target2_11)/amp_target2_11

        #Add the tests results to the dataframes:-------------------------------

        #Add the PBWE - polar 00 test results to the corresponding dataframes:
        dataframe_pbwe_polar00_percentage_echoes_detection.loc[snr,dist] = pc_echoes_detection_pbwe00
        dataframe_pbwe_polar00_distance_mean_error.loc[snr,dist] = np.mean(error_distance_pbwe00)
        dataframe_pbwe_polar00_distance_std_error.loc[snr,dist] = np.std(error_distance_pbwe00)
        dataframe_pbwe_polar00_amplitude_1_mean_error.loc[snr,dist] = np.mean(error_amplitude_target1_pbwe00)
        dataframe_pbwe_polar00_amplitude_1_std_error.loc[snr,dist] = np.std(error_amplitude_target1_pbwe00)
        dataframe_pbwe_polar00_amplitude_2_mean_error.loc[snr,dist] = np.mean(error_amplitude_target2_pbwe00)
        dataframe_pbwe_polar00_amplitude_2_std_error.loc[snr,dist] = np.std(error_amplitude_target2_pbwe00)

        #Add the PBWE - polar 11 test results to the corresponding dataframes:
        dataframe_pbwe_polar11_percentage_echoes_detection.loc[snr,dist] = pc_echoes_detection_pbwe11
        dataframe_pbwe_polar11_distance_mean_error.loc[snr,dist] = np.mean(error_distance_pbwe11)
        dataframe_pbwe_polar11_distance_std_error.loc[snr,dist] = np.std(error_distance_pbwe11)
        dataframe_pbwe_polar11_amplitude_1_mean_error.loc[snr,dist] = np.mean(error_amplitude_target1_pbwe11)
        dataframe_pbwe_polar11_amplitude_1_std_error.loc[snr,dist] = np.std(error_amplitude_target1_pbwe11)
        dataframe_pbwe_polar11_amplitude_2_mean_error.loc[snr,dist] = np.mean(error_amplitude_target2_pbwe11)
        dataframe_pbwe_polar11_amplitude_2_std_error.loc[snr,dist] = np.std(error_amplitude_target2_pbwe11)

#Format columns and index names for export:-------------------------------------

data_columns = ['delta = '+str(dist)+' (m)' for dist in list_dist_targets]
data_index = ['SNR = '+str(snr)+' (dB)' for snr in list_snr_levels]

#For the PBWE - polar 00:

dataframe_pbwe_polar00_percentage_echoes_detection.index = data_index
dataframe_pbwe_polar00_distance_mean_error.index = data_index
dataframe_pbwe_polar00_distance_std_error.index = data_index
dataframe_pbwe_polar00_amplitude_1_mean_error.index = data_index
dataframe_pbwe_polar00_amplitude_1_std_error.index = data_index
dataframe_pbwe_polar00_amplitude_2_mean_error.index = data_index
dataframe_pbwe_polar00_amplitude_2_std_error.index = data_index

dataframe_pbwe_polar00_percentage_echoes_detection.columns = data_columns
dataframe_pbwe_polar00_distance_mean_error.columns = data_columns
dataframe_pbwe_polar00_distance_std_error.columns = data_columns
dataframe_pbwe_polar00_amplitude_1_mean_error.columns = data_columns
dataframe_pbwe_polar00_amplitude_1_std_error.columns = data_columns
dataframe_pbwe_polar00_amplitude_2_mean_error.columns = data_columns
dataframe_pbwe_polar00_amplitude_2_std_error.columns = data_columns

#For the PBWE - polar 11:

dataframe_pbwe_polar11_percentage_echoes_detection.index = data_index
dataframe_pbwe_polar11_distance_mean_error.index = data_index
dataframe_pbwe_polar11_distance_std_error.index = data_index
dataframe_pbwe_polar11_amplitude_1_mean_error.index = data_index
dataframe_pbwe_polar11_amplitude_1_std_error.index = data_index
dataframe_pbwe_polar11_amplitude_2_mean_error.index = data_index
dataframe_pbwe_polar11_amplitude_2_std_error.index = data_index

dataframe_pbwe_polar11_percentage_echoes_detection.columns = data_columns
dataframe_pbwe_polar11_distance_mean_error.columns = data_columns
dataframe_pbwe_polar11_distance_std_error.columns = data_columns
dataframe_pbwe_polar11_amplitude_1_mean_error.columns = data_columns
dataframe_pbwe_polar11_amplitude_1_std_error.columns = data_columns
dataframe_pbwe_polar11_amplitude_2_mean_error.columns = data_columns
dataframe_pbwe_polar11_amplitude_2_std_error.columns = data_columns

#Export the Markdown report:----------------------------------------------------

with open(test_report_path,'w') as file_report:

    #Title
    file_report.write('# PBWE performance test report\r\n')

    #Test scenario:
    file_report.write('## Test scenario\r\n')
    file_report.write('### Scenario\r\n')
    file_report.write('The test is performed on synthetic radar signals (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017):\r\n')
    file_report.write('* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.\r\n')
    file_report.write('* This radar can transmit and receive signals in two linear polarizations named 0 and 1. This leads to 4 polirization channels: 2 co-polar named 00 and 11 (emission and reception in the same polarization), and 2 cross-polar named 01 and 10 (emission and reception in different polarizations). For simplicity, we will only consider the co-polar channels 00 and 11 in this example.\r\n')
    file_report.write('* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.\r\n')
    file_report.write('* Two targets in free-space are seperated by 5 cm, slightly below the radar free-space resolution. These targets generate echoes of given complex amplitudes in the radar signal, or complex sine-waves in the measured spectrum.\r\n')
    file_report.write('* The measured spectrum is corrupted by a white-noise of standard deviation 10X smaller than the complex sine-waves amplitudes.\r\n')
    file_report.write('### Scenario parameters\r\n')
    file_report.write('* Amplitude of the 1st echo in polarimetric channel 00: '+str(amp_target1_00)+'\r\n')
    file_report.write('* Amplitude of the 2nd echo in polarimetric channel 00: '+str(amp_target2_00)+'\r\n')
    file_report.write('* Amplitude of the 1st echo in polarimetric channel 11: '+str(amp_target1_11)+'\r\n')
    file_report.write('* Amplitude of the 2nd echo in polarimetric channel 11: '+str(amp_target2_11)+'\r\n')
    file_report.write('* Distance (m) between the 1st target and the radar:  '+str(dist_target1)+'\r\n')
    file_report.write('* Peak detection threshold on the amplitude of echoes: '+str(detection_level)+'\r\n')

    #Test parameters:
    file_report.write('## Test parameters\r\n')
    file_report.write('* PyBWE version: '+version('PyBWE')+'\r\n')
    file_report.write('* Tested function: PyPBWE.PBWE\r\n')
    file_report.write('* Distances between targets **delta** (m): '+str(list_dist_targets)+'\r\n')
    file_report.write('* **SNR** levels (dB): '+str(list_snr_levels)+'\r\n')
    file_report.write('* Number of noise cases: '+str(nb_noise_case)+'\r\n')

    #Percentage of echoes detection:
    file_report.write('## Percentage of echoes detection\r\n')
    file_report.write('### Polar channel 00\r\n')
    file_report.write(dataframe_pbwe_polar00_percentage_echoes_detection.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11\r\n')
    file_report.write(dataframe_pbwe_polar11_percentage_echoes_detection.to_markdown()+'\r\n')

    #Error on the distance between targets:
    file_report.write('## Error on the distance between targets\r\n')
    file_report.write('### Polar channel 00 - Mean (m)\r\n')
    file_report.write(dataframe_pbwe_polar00_distance_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 00 - STD (m)\r\n')
    file_report.write(dataframe_pbwe_polar00_distance_std_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - Mean (m)\r\n')
    file_report.write(dataframe_pbwe_polar11_distance_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - STD (m)\r\n')
    file_report.write(dataframe_pbwe_polar11_distance_std_error.to_markdown()+'\r\n')

    #Error on the amplitude of the 1st echo:
    file_report.write('## Error on the amplitude of the 1st echo\r\n')
    file_report.write('### Polar channel 00 - Mean (%)\r\n')
    file_report.write(dataframe_pbwe_polar00_amplitude_1_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 00 - STD (%)\r\n')
    file_report.write(dataframe_pbwe_polar00_amplitude_1_std_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - Mean (%)\r\n')
    file_report.write(dataframe_pbwe_polar11_amplitude_1_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - STD (%)\r\n')
    file_report.write(dataframe_pbwe_polar11_amplitude_1_std_error.to_markdown()+'\r\n')

    #Error on the amplitude of the 2nd echo:
    file_report.write('## Error on the amplitude of the 2nd echo\r\n')
    file_report.write('### Polar channel 00 - Mean (%)\r\n')
    file_report.write(dataframe_pbwe_polar00_amplitude_2_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 00 - STD (%)\r\n')
    file_report.write(dataframe_pbwe_polar00_amplitude_2_std_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - Mean (%)\r\n')
    file_report.write(dataframe_pbwe_polar11_amplitude_2_mean_error.to_markdown()+'\r\n')
    file_report.write('### Polar channel 11 - STD (%)\r\n')
    file_report.write(dataframe_pbwe_polar11_amplitude_2_std_error.to_markdown()+'\r\n')