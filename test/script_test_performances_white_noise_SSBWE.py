################################################################################
#HEADER

#This script allows you to test the performances of the SSBWE package facing
#different levels of white-noise, to check for regressions.
#It will generate a report containing the test results.

#The test consists in measuring the performances of the SSBWE on a radar signal
#containing 2 echoes from 2 targets, separated by distances below or close to
#the radar's resolution.

#The following test parameters are user-defined:
#   -The list of distances between targets to be tested (m).
#   -The list of SNR to be tested (dB).
#   -The number of noise cases

#The State-Space model's order can either be forced to 2, or estimated using
#AIC (to do so, set the order to 0).

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
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. These targets generate echoes of given complex
#    amplitudes in the radar's signal, or complex sine-waves in the measured
#    spectrum.
#   -The measured spectrum is corrupted by a white-noise of standard deviation
#    10X smaller than the complex sine-waves' amplitudes.

#The parameters used for the SSBWE are the default ones.

#References: Oudart et al. (2021)

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
import pandas as pd
import os
from math import pi
from scipy.signal import hilbert,find_peaks

import PySSBWE

#Test parameters:---------------------------------------------------------------

#List of distances between targets to be tested (m):
list_dist_targets = [0.04,0.06,0.08,0.1,0.12]

#List of SNR (dB):
list_snr_levels = [6,14,20,40,60]

#Number of noise case to be tested for each distance:
nb_noise_case = 1000

#Order of the model (forcing 2 or 0 = AIC estimation):
param_order = 2

#Scenario parameters:-----------------------------------------------------------

#Amplitude of the 2 echoes corresponding to the 2 targets:
amp_target1 = 1
amp_target2 = 1

#Distance (m) between the 1st radar target in free-space (returning one echo)
#and the radar:
dist_target1 = 1

#Peak detection threshold on the amplitude of echoes:
detection_level = 0.5

#Retrieve the test path:--------------------------------------------------------

#Test directory path:
test_dir_path = os.path.dirname(__file__)

#Test report path:
test_report_path = os.path.join(test_dir_path,'PySSBWE_Report_test_performances_white_noise.md')

#Initialize the test results dataframes:----------------------------------------

#Columns and index names:
data_columns = ['delta = '+str(dist)+' (m)' for dist in list_dist_targets]
data_index = ['SNR = '+str(snr)+' (dB)' for snr in list_snr_levels]

#For the SSBWE - method 1:
dataframe_ssbwe_method1_percentage_echoes_detection = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_distance_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_distance_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_amplitude_1_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_amplitude_1_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_amplitude_2_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method1_amplitude_2_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)

#For the SSBWE - method 2:
dataframe_ssbwe_method2_percentage_echoes_detection = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_distance_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_distance_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_amplitude_1_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_amplitude_1_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_amplitude_2_mean_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)
dataframe_ssbwe_method2_amplitude_2_std_error = pd.DataFrame(0.0,index=list_snr_levels,columns=list_dist_targets)


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
    spec_vect = (amp_target1*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))

    #Only keep the real part of the spectrum (In-phase component):
    spec_vect = np.real(spec_vect)

    #Calculate the power of the signal (dB):
    power_sig = 10*np.log10(np.mean(spec_vect**2))

    #Iterate on the SNR values:
    for snr in list_snr_levels:

        #Initialize the 2 echoes detection counter:
        nb_echoes_detection_ssbwe1 = 0
        nb_echoes_detection_ssbwe2 = 0

        #Initialize the lists to stock all noise cases results:
        list_distance_ssbwe1 = []
        list_amplitude_target1_ssbwe1 = []
        list_amplitude_target2_ssbwe1 = []
        list_distance_ssbwe2 = []
        list_amplitude_target1_ssbwe2 = []
        list_amplitude_target2_ssbwe2 = []

        #Calculate the power of the noise (dB) to obtain the right SNR:
        power_noise = power_sig-snr
        power_noise = power_sig-snr

        #Iterate on the noise case:
        for idx_case in range(nb_noise_case):

            #Display the test progress:
            count_iterations += 1
            print('Test in progress ... '+str(round(100*count_iterations/nb_iterations,3))+'%')

            #Add a white noise to this spectrum signal to obtain the right SNR:
            rng = np.random.RandomState(idx_case)
            wn_vect = rng.normal(0,np.sqrt(10**(power_noise/10)),spec_vect.shape)
            spec_vect_wn = spec_vect + wn_vect

            #Reconstruct a complex signal with the Hilbert transform:
            spec_vect_wn = np.conjugate(hilbert(spec_vect_wn))[::2]

            #Application of the SSBWE:
            output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_wn,df=5e6,extra_factor=3,zp_factor=10,side_cut=True,order=param_order)

            #Testing the SSBWE - method 1:--------------------------------------

            #Detection of the 2 echoes:
            peaks_ssbwe_1 = find_peaks(abs(output_ssbwe_1),height=detection_level)[0]

            #Check that 2 echoes are detected:
            if len(peaks_ssbwe_1)==2:

                #Increment by 1 the number of echoes detection:
                nb_echoes_detection_ssbwe1 += 1

                #Stock the distance between echoes:
                list_distance_ssbwe1 += [0.5*abs(time_ssbwe_vect[peaks_ssbwe_1[0]]-time_ssbwe_vect[peaks_ssbwe_1[1]])*3e8]

                #Stock the amplitude of the 1st echo:
                list_amplitude_target1_ssbwe1 += [abs(output_ssbwe_1[peaks_ssbwe_1[0]])]

                #Stock the amplitude of the 2nd echo:
                list_amplitude_target2_ssbwe1 += [abs(output_ssbwe_1[peaks_ssbwe_1[1]])]

            #Testing the SSBWE - method 2:--------------------------------------

            #Detection of the 2 echoes:
            peaks_ssbwe_2 = find_peaks(abs(output_ssbwe_2),height=detection_level)[0]

            #Check that 2 echoes are detected:
            if len(peaks_ssbwe_2)==2:

                #Increment by 1 the number of echoes detection:
                nb_echoes_detection_ssbwe2 += 1

                #Stock the distance between echoes:
                list_distance_ssbwe2 += [0.5*abs(time_ssbwe_vect[peaks_ssbwe_2[0]]-time_ssbwe_vect[peaks_ssbwe_2[1]])*3e8]

                #Stock the amplitude of the 1st echo:
                list_amplitude_target1_ssbwe2 += [abs(output_ssbwe_2[peaks_ssbwe_2[0]])]

                #Stock the amplitude of the 2nd echo:
                list_amplitude_target2_ssbwe2 += [abs(output_ssbwe_2[peaks_ssbwe_2[1]])]


        #Tests calculations:----------------------------------------------------

        #Calculate the percentage of echoes detection:
        pc_echoes_detection_ssbwe1 = 100*nb_echoes_detection_ssbwe1/nb_noise_case
        pc_echoes_detection_ssbwe2 = 100*nb_echoes_detection_ssbwe2/nb_noise_case

        #Calculate the error on the distance between echoes:
        error_distance_ssbwe1 = np.array(list_distance_ssbwe1)-dist
        error_distance_ssbwe2 = np.array(list_distance_ssbwe2)-dist

        #Calculate the error on the amplitude of the 1st echo:
        error_amplitude_target1_ssbwe1 = 100*(np.array(list_amplitude_target1_ssbwe1)-amp_target1)/amp_target1
        error_amplitude_target1_ssbwe2 = 100*(np.array(list_amplitude_target1_ssbwe2)-amp_target1)/amp_target1

        #Calculate the error on the amplitude of the 2nd echo:
        error_amplitude_target2_ssbwe1 = 100*(np.array(list_amplitude_target2_ssbwe1)-amp_target2)/amp_target2
        error_amplitude_target2_ssbwe2 = 100*(np.array(list_amplitude_target2_ssbwe2)-amp_target2)/amp_target2

        #Add the tests results to the dataframes:-------------------------------

        #The SSBWE - method 1 test results to the corresponding dataframes:
        dataframe_ssbwe_method1_percentage_echoes_detection.loc[snr,dist] = pc_echoes_detection_ssbwe1
        dataframe_ssbwe_method1_distance_mean_error.loc[snr,dist] = np.mean(error_distance_ssbwe1)
        dataframe_ssbwe_method1_distance_std_error.loc[snr,dist] = np.std(error_distance_ssbwe1)
        dataframe_ssbwe_method1_amplitude_1_mean_error.loc[snr,dist] = np.mean(error_amplitude_target1_ssbwe1)
        dataframe_ssbwe_method1_amplitude_1_std_error.loc[snr,dist] = np.std(error_amplitude_target1_ssbwe1)
        dataframe_ssbwe_method1_amplitude_2_mean_error.loc[snr,dist] = np.mean(error_amplitude_target2_ssbwe1)
        dataframe_ssbwe_method1_amplitude_2_std_error.loc[snr,dist] = np.std(error_amplitude_target2_ssbwe1)

        #The SSBWE - method 2 test results to the corresponding dataframes:
        dataframe_ssbwe_method2_percentage_echoes_detection.loc[snr,dist] = pc_echoes_detection_ssbwe2
        dataframe_ssbwe_method2_distance_mean_error.loc[snr,dist] = np.mean(error_distance_ssbwe2)
        dataframe_ssbwe_method2_distance_std_error.loc[snr,dist] = np.std(error_distance_ssbwe2)
        dataframe_ssbwe_method2_amplitude_1_mean_error.loc[snr,dist] = np.mean(error_amplitude_target1_ssbwe2)
        dataframe_ssbwe_method2_amplitude_1_std_error.loc[snr,dist] = np.std(error_amplitude_target1_ssbwe2)
        dataframe_ssbwe_method2_amplitude_2_mean_error.loc[snr,dist] = np.mean(error_amplitude_target2_ssbwe2)
        dataframe_ssbwe_method2_amplitude_2_std_error.loc[snr,dist] = np.std(error_amplitude_target2_ssbwe2)

#Export the Excel report:-------------------------------------------------------

with pd.ExcelWriter(test_report_path) as writer:

    #Excel sheets corresponding to the SSBWE - method 1 tests:------------------

    dataframe_ssbwe_method1_percentage_echoes_detection.to_excel(writer, sheet_name='SSBWE_v1_echoes_detection',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_echoes_detection']
    worksheet['A1'] = 'SSBWE - method 1: Detection of the 2 echoes (%)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_distance_mean_error.to_excel(writer, sheet_name='SSBWE_v1_mean_distance_error',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_mean_distance_error']
    worksheet['A1'] = 'SSBWE - method 1: Mean error on the distance between targets (m)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_distance_std_error.to_excel(writer, sheet_name='SSBWE_v1_STD_distance_error',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_STD_distance_error']
    worksheet['A1'] = 'SSBWE - method 1: STD of the error on the distance between targets (m)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_amplitude_1_mean_error.to_excel(writer, sheet_name='SSBWE_v1_mean_amplitude_error_1',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_mean_amplitude_error_1']
    worksheet['A1'] = 'SSBWE - method 1: Mean error on the amplitude of the 1st echo (%)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_amplitude_1_std_error.to_excel(writer, sheet_name='SSBWE_v1_STD_amplitude_error_1',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_STD_amplitude_error_1']
    worksheet['A1'] = 'SSBWE - method 1: STD of the error on the amplitude of the 1st echo (%)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_amplitude_2_mean_error.to_excel(writer, sheet_name='SSBWE_v1_mean_amplitude_error_2',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_mean_amplitude_error_2']
    worksheet['A1'] = 'SSBWE - method 1: Mean error on the amplitude of the 2nd echo (%)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

    dataframe_ssbwe_method1_amplitude_2_std_error.to_excel(writer, sheet_name='SSBWE_v1_STD_amplitude_error_2',index_label='SNR (dB)',startrow = 2)
    worksheet = writer.sheets['SSBWE_v1_STD_amplitude_error_2']
    worksheet['A1'] = 'SSBWE - method 1: STD of the error on the amplitude of the 2nd echo (%)'
    worksheet['A1'].font = Font(bold=True,color="FFFFFF")
    worksheet['A1'].fill = PatternFill(start_color='000000',end_color='000000',fill_type='solid')
    worksheet.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(list_dist_targets)+1)
    worksheet['B2'] = 'Distance between targets (m)'
    worksheet['B2'].font = Font(bold=True)
    worksheet.merge_cells(start_row=2, start_column=2, end_row=2, end_column=len(list_dist_targets)+1)

#Export the Markdown report:----------------------------------------------------

with open(test_report_path,'w') as file_report:

    #Title
    file_report.write('# SSBWE performance test report\r\n')

    #Test scenario:
    file_report.write('## Test scenario\r\n')
    file_report.write('### Scenario\r\n')
    file_report.write('The test is performed on synthetic radar signals (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017):\r\n')
    file_report.write('* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.\r\n')
    file_report.write('* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.\r\n')
    file_report.write('* Two targets in free-space are seperated by 5 cm, slightly below the radar free-space resolution. These targets generate echoes of given complex amplitudes in the radar signal, or complex sine-waves in the measured spectrum.\r\n')
    file_report.write('* The measured spectrum is corrupted by a white-noise of standard deviation 10X smaller than the complex sine-waves amplitudes.\r\n')
    file_report.write('### Scenario parameters\r\n')
    file_report.write('* Amplitude of the 1st echo: '+str(amp_target1)+'\r\n')
    file_report.write('* Amplitude of the 2nd echo: '+str(amp_target2)+'\r\n')
    file_report.write('* Distance (m) between the 1st target and the radar:  '+str(dist_target1)+'\r\n')
    file_report.write('* Peak detection threshold on the amplitude of echoes: '+str(detection_level)+'\r\n')

    #Test parameters:
    file_report.write('## Test parameters\r\n')
    file_report.write('* PyBWE version: '+version('PyBWE')+'\r\n')
    file_report.write('* Tested function: PySSBWE.SSBWE\r\n')
    file_report.write('* Distances between targets **delta** (m): '+str(list_dist_targets)+'\r\n')
    file_report.write('* **SNR** levels (dB): '+str(list_snr_levels)+'\r\n')
    file_report.write('* Number of noise cases: '+str(nb_noise_case)+'\r\n')

    #Percentage of echoes detection:
    file_report.write('## Percentage of echoes detection\r\n')
    file_report.write('### Observability matrix method (1)\r\n')
    file_report.write(dataframe_ssbwe_method1_percentage_echoes_detection.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2)\r\n')
    file_report.write(dataframe_ssbwe_method2_percentage_echoes_detection.to_markdown()+'\r\n')

    #Error on the distance between targets:
    file_report.write('## Error on the distance between targets\r\n')
    file_report.write('### Observability matrix method (1) - Mean (m)\r\n')
    file_report.write(dataframe_ssbwe_method1_distance_mean_error.to_markdown()+'\r\n')
    file_report.write('### Observability matrix method (1) - STD (m)\r\n')
    file_report.write(dataframe_ssbwe_method1_distance_std_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - Mean (m)\r\n')
    file_report.write(dataframe_ssbwe_method2_distance_mean_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - STD (m)\r\n')
    file_report.write(dataframe_ssbwe_method2_distance_std_error.to_markdown()+'\r\n')

    #Error on the amplitude of the 1st echo:
    file_report.write('## Error on the amplitude of the 1st echo\r\n')
    file_report.write('### Observability matrix method (1) - Mean (%)\r\n')
    file_report.write(dataframe_ssbwe_method1_amplitude_1_mean_error.to_markdown()+'\r\n')
    file_report.write('### Observability matrix method (1) - STD (%)\r\n')
    file_report.write(dataframe_ssbwe_method1_amplitude_1_std_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - Mean (%)\r\n')
    file_report.write(dataframe_ssbwe_method2_amplitude_1_mean_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - STD (%)\r\n')
    file_report.write(dataframe_ssbwe_method2_amplitude_1_std_error.to_markdown()+'\r\n')

    #Error on the amplitude of the 2nd echo:
    file_report.write('## Error on the amplitude of the 2nd echo\r\n')
    file_report.write('### Observability matrix method (1) - Mean (%)\r\n')
    file_report.write(dataframe_ssbwe_method1_amplitude_2_mean_error.to_markdown()+'\r\n')
    file_report.write('### Observability matrix method (1) - STD (%)\r\n')
    file_report.write(dataframe_ssbwe_method1_amplitude_2_std_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - Mean (%)\r\n')
    file_report.write(dataframe_ssbwe_method2_amplitude_2_mean_error.to_markdown()+'\r\n')
    file_report.write('### Controllability matrix method (2) - STD (%)\r\n')
    file_report.write(dataframe_ssbwe_method2_amplitude_2_std_error.to_markdown()+'\r\n')