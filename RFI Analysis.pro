
;**************************************************************************************************
; Author:       R. Lord
; Organisation: SKA SA
; Ver:          1.0
;**************************************************************************************************


DataFilename             = 'tt382.tiq'
RBW                      = 100                          ;Resolution Bandwidth [Hz]
subtitle                 = 'Active Baseline Measurement #1'         ;Description of measurement
polarisation             = 0                           ;0=V polarisation, 1=H polarisation
chamber_type             = 1                           ;0=Anechoic Chamber, 1=Reverb Chamber


;---Calibration Data
AntGainFilenameNumber    = 1                           ;Refer to numbers below
LnaGainFilenameNumber    = 2                           ;Refer to numbers below
CableLossFilenameNumber  = 15                           ;Refer to numbers below


AntGainFilenameList      = ['Antenna_MESA_GLPDA'                        , $ ;#1
                            'Antenna_MESA_KLPDA1'                       , $ ;#2
                            'Antenna_MESA_KLPDA3'                       , $ ;#3
                            'Antenna_Houwteq_EMCO3115_FRH_Small'        , $ ;#4
                            'Antenna_Houwteq_EMCO3106B_FRH_Large'       , $ ;#5
                            'Antenna_Houwteq_EMCO3141_BiConiLog_FLPDA'  , $ ;#6
                            'No_Antenna']                                   ;#7 (Use this option for direct injection into RTSA)

LnaGainFilenameList      = ['LNA_1'                                     , $ ;#1
                            'LNA_2'                                     , $ ;#2
                            'LNA_1_plus_2'                              , $ ;#3
                            'No_LNA']                                       ;#4 (Use this option for passive antenna)

CableLossFilenameList    = ['CableLoss_Houwteq_1a'                      , $ ;#1
                            'CableLoss_Houwteq_1b'                      , $ ;#2
                            'CableLoss_Houwteq_1c'                      , $ ;#3
                            'CableLoss_Houwteq_1a_1b'                   , $ ;#4
                            'CableLoss_Houwteq_2a'                      , $ ;#5
                            'CableLoss_Houwteq_2b'                      , $ ;#6
                            'CableLoss_Houwteq_2c'                      , $ ;#7
                            'CableLoss_Houwteq_2a_2b_2c'                , $ ;#8
                            'CableLoss_Houwteq_3a'                      , $ ;#9
                            'CableLoss_Houwteq_3b'                      , $ ;#10
                            'CableLoss_Houwteq_3c'                      , $ ;#11
                            'CableLoss_Houwteq_3a_3b_3c'                , $ ;#12
                            'CableLoss_Houwteq_4a'                      , $ ;#13
                            'CableLoss_Houwteq_4b'                      , $ ;#14
                            'CableLoss_Pinelands_5a_5b'                 , $ ;#15
                            'CableLoss_Pinelands_6a_6b'                 , $ ;#16
                            'CableLoss_Pinelands_6a_6b_6c'              , $ ;#17
                            'CableLoss_Pinelands_7a'                    , $ ;#18
                            'No_Cable']                                     ;#19


;---Other
range                    = 1.0                         ;Distance from DUT to receiving antenna [m]
culprit_distance         = 1.0                         ;Distance between Antenna Focus and DUT [m] (E.g. 10m or 30m for AP, 1m for Receiver and Digitiser)
enclosure_attenuation    = 0.0                         ;Additional attenuation provided by enclosure [dB]
peak_hold                = 0                           ;0=no, 1=yes
remove_dc_offset         = 0                           ;0=no, 1=yes
write_spectrum_to_csv    = 0                           ;0=no, 1=yes (write spectrum to csv file)
transient_RFI_detection  = 0                           ;0=no, 1=yes (Find transient RFI by calculating normalised sigma for every frequency bin)
fraction_to_process      = 1.0                         ;Fraction of data to process (1.0 = 100%)
FFT_power_of_2           = 0                           ;0=no, 1=yes (if yes, FFT length will be rounded to nearest power of 2, but then RBW will be approximate)
num_peaks                = 20                          ;Max number of RFI peaks to identify
threshold_margin         = -3.0                        ;Margin within Spectral Line threshold level to identify RFI spike [dB]
input_dir                = 'data\'
output_dir               = 'plots\'
calibration_dir_ant      = 'Equipment_Database\Passive_Antennas\'
calibration_dir_lna      = 'Equipment_Database\LNAs\'
calibration_dir_cable    = 'Equipment_Database\Cables\'
plant_signals            = 0                           ;0=no, 1=yes
freq_signal_at_threshold = 856e6
freq_signal_at_Xsigma    = 856.4e6
Xsigma                   = 7.0
remove_datafile          = 1


;---RFI Flag List
generate_rfi_flaglist    = 0                           ;0=no, 1=yes
apply_rfi_flaglist       = 0                           ;0=no, 1=yes
cutoff_above_noisefloor  = 2.0                         ;everything above this level will be classified RFI [dB]
rfi_flaglist_filename    = 'test.csv'
rfi_flaglist_folder      = 'RFI_Flaglist\'


;---Reverberation Chamber
antenna_efficiency       = 0.8
CCF_filename             = 'Equipment_Database\Reverb_Chamber\ACF.csv'


;---Plot parameters
replot_old_data          = 0                           ;0=no, 1=yes
yaxis_bottom_add         = 0.0                         ;raises bottom of y-axis by this amount
yaxis_top_add            = 0.0                         ;drops top of y-axis by this amount
xaxis_centre             = 500e6                       ;zoom into spectrum at this centre frequency [Hz]
xaxis_span               = 1e10                        ;zoom into spectrum with this span (set to large value to see everything) [Hz]
show_threshold           = 1                           ;show threshold level (red dotted line)


;**************************************************************************************************
;**************************************************************************************************


;---Compile routines
Resolve_Routine, 'Subroutines'
headerlines = 2
show_plot = 1

if (remove_datafile eq 1) then begin
   answer = ''
   read, answer, prompt='Please press y to confirm that raw data should be deleted:'
   if (strlowcase(answer) ne 'y') then begin
      print, 'Script aborted'
      retall
   endif
endif

if (replot_old_data eq 0) then begin
   start_time = systime(1)
   prev_filename = DataFilename
   prev_RBW = RBW

   process_data, DataFilename, RBW, polarisation, AntGainFilenameNumber, AntGainFilenameList, $
                 LnaGainFilenameNumber, LnaGainFilenameList, CableLossFilenameNumber, $
                 CableLossFilenameList, AntGainFilename, LnaGainFilename, CableLossFilename, $
                 range, culprit_distance, enclosure_attenuation, peak_hold, remove_dc_offset, $
                 transient_RFI_detection, fraction_to_process, FFT_power_of_2, num_peaks, $
                 threshold_margin, headerlines, input_dir, output_dir, calibration_dir_ant, $
                 calibration_dir_lna, calibration_dir_cable, plant_signals, freq_signal_at_threshold, $
                 freq_signal_at_Xsigma, Xsigma, Spectrum, S, FreqAxis, noise_floor_E, $
                 SpectralE_interp, ContinuumE_interp, AntennaGain, Glna, CableLoss, SamplingFrequency, $
                 fft_length, num_ffts, Frequency, AcqBW, date, Preamp, rfi_SNR, RFAttenuation, $
                 rfi_frequency, rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, $
                 apply_rfi_flaglist, generate_rfi_flaglist, cutoff_above_noisefloor, $
                 rfi_flaglist_filename, rfi_flaglist_folder, chamber_type, antenna_efficiency, $
                 CCF_filename, CCF

endif else start_time = systime(1)

;---Test if replot_old_data is correctly set
test_replot_old_data, prev_filename, DataFilename, RBW, prev_RBW

save_result, Spectrum, S, FreqAxis, noise_floor_E, SpectralE_interp, ContinuumE_interp, AntennaGain, $
             Glna, CableLoss, DataFilename, AntGainFilename, LnaGainFilename, CableLossFilename, $
             subtitle, write_spectrum_to_csv, yaxis_bottom_add, yaxis_top_add, xaxis_centre, $
             xaxis_span, SamplingFrequency, fft_length, num_ffts, Frequency, AcqBW, output_dir, $
             date, transient_RFI_detection, polarisation,  Preamp, RFAttenuation, range, start_time, $
             fraction_to_process, num_peaks, plant_signals, freq_signal_at_threshold, $
             freq_signal_at_Xsigma, threshold_margin, rfi_frequency, rfi_plus_noise_power, $
             rfi_power, rfi_margin, rfi_3dB_bandwidth, rfi_SNR, show_plot, remove_dc_offset, $
             culprit_distance, enclosure_attenuation, apply_rfi_flaglist, generate_rfi_flaglist, $
             cutoff_above_noisefloor, rfi_flaglist_folder, rfi_flaglist_filename, chamber_type, $
             antenna_efficiency, CCF_filename, CCF, show_threshold

cleanup, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, AntennaGain, $
         Glna, CableLoss, S, CCF

if (remove_datafile eq 1) then begin
   print, 'Warning: Deleting raw data file ' + DataFilename
   File_Delete, input_dir + DataFilename
endif


End


;**************************************************************************************************
;**************************************************************************************************


