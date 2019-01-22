
;**************************************************************************************************
; Author:       R. Lord
; Organisation: SKA SA
; Ver:          1.0
;**************************************************************************************************


DataFilename_base        = 'tt'
DataFilename_start_num   = 381
DataFilename_end_num     = 381
RBW_list                 = [1,10,100]
subtitle                 = 'DUT Qualification'         ;Description of measurement
polarisation             = 0                           ;0=V polarisation, 1=H polarisation
chamber_type             = 1                           ;0=Anechoic Chamber, 1=Reverb Chamber
output_dir               = 'plots\'
calibration_dir_ant      = 'Equipment_Database\Passive_Antennas\'
AntGainFilenameNumber    = 1                                                             ;Refer to numbers below
AntGainFilenameList      = ['Antenna_MESA_GLPDA'                        , $ ;#1
                            'Antenna_MESA_KLPDA1'                       , $ ;#2
                            'Antenna_MESA_KLPDA3'                       , $ ;#3
                            'Antenna_Houwteq_EMCO3115_FRH_Small'        , $ ;#4
                            'Antenna_Houwteq_EMCO3106B_FRH_Large'       , $ ;#5
                            'Antenna_Houwteq_EMCO3141_BiConiLog_FLPDA'  , $ ;#6
                            'No_Antenna']                                   ;#7 (Use this option for direct injection into RTSA)
range                    = 1.0                         ;Distance from DUT to receiving antenna [m]
culprit_distance         = 1.0                         ;Distance between Antenna Focus and DUT [m] (E.g. 10m or 30m for AP, 1m for Receiver and Digitiser)
enclosure_attenuation    = 0.0                         ;Additional attenuation provided by enclosure [dB]
write_spectrum_to_csv    = 0                           ;0=no, 1=yes (write spectrum to csv file)
num_peaks                = 20                          ;Max number of RFI peaks to identify
threshold_margin         = -3.0                        ;Margin within Spectral Line threshold level to identify RFI spike [dB]
show_threshold           = 1                           ;show threshold level (red dotted line)


;**************************************************************************************************
;**************************************************************************************************


;---Compile routines
Resolve_Routine, 'Subroutines'
Resolve_Routine, 'strsplit', /is_function
Resolve_Routine, 'interpol', /is_function
headerlines = 2
show_plot = 0
maxbins = 100000000L
print, ''


;---Error checking
if (DataFilename_end_num lt DataFilename_start_num) then begin
   print, 'Error: DataFilename_end_num must be larger than DataFilename_start_num!'
   retall
endif

;---Get Antenna Gain filename
if ((AntGainFilenameNumber lt 1) OR (AntGainFilenameNumber gt n_elements(AntGainFilenameList))) then begin
   print, 'Error: Number supplied for Antenna Gain filename must be >=1 and <=' + strtrim(n_elements(AntGainFilenameList), 2)
   retall
endif
AntGainFilename = AntGainFilenameList[AntGainFilenameNumber - 1] + '.csv'

;---Start loop
for j=0, (n_elements(RBW_list) - 1) do begin
   start_time = systime(1)
   RBW = float(RBW_list[j])
   print, 'Stitching spectra with RBW = ' + strtrim(RBW_list[j], 2) + ' Hz'
   wait, 0.1

   Folder_list = DataFilename_base + strtrim(Indgen(DataFilename_end_num - DataFilename_start_num + 1) + DataFilename_start_num, 2)

   get_stitch_spectra_filenames, output_dir, Folder_list, RBW, stitch_filenames

   is_array = size(stitch_filenames)

   if (is_array[0] eq 1) then begin

      get_stitch_spectra_parameters, Folder_list, stitch_filenames, NumberSamples, CentreFrequency, SamplingFrequency, AcqBW

      get_stitch_spectra_Spectrum, Spectrum, FreqAxis, stitch_filenames, SamplingFrequency, NumberSamples, CentreFrequency, $
                                   AcqBW, shrinkfactor, maxbins

      read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation, AntennaGain

      get_stitch_spectra_values, FreqAxis, CentreFrequency, SamplingFrequency, NumberSamples, fft_length, AcqBW, $
                                 Frequency, shrinkfactor

      get_noise_floor, Spectrum, FreqAxis, SamplingFrequency, fft_length, noise_floor_E, AcqBW

      get_RFI_thresholds, FreqAxis, range, culprit_distance, enclosure_attenuation, ContinuumE_interp, SpectralE_interp

      SpectralE_interp = 10.0 * alog10((10.0 ^ (Temporary(SpectralE_interp) / 10.0)) + (10.0 ^ (noise_floor_E / 10.0)))

      get_stitch_spectra_plot, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, Folder_list, $
                               SamplingFrequency, fft_length, Frequency, AcqBW, subtitle, output_dir, DataFilename, $
                               AntennaGain, filename_base, show_plot, show_threshold

      find_rfi_spikes, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, AntennaGain, num_peaks, threshold_margin, $
                       SamplingFrequency, fft_length, 0, rfi_frequency, rfi_plus_noise_power, rfi_power, rfi_margin, $
                       rfi_3dB_bandwidth, rfi_SNR, Frequency, AcqBW

      stitch_spectra_summary, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, filename_base, Folder_list, $
                              AntGainFilename, AntennaGain, polarisation, Frequency, AcqBW, SamplingFrequency, $
                              fft_length, subtitle, range, culprit_distance, enclosure_attenuation, start_time, $
                              num_peaks, threshold_margin, rfi_frequency, rfi_plus_noise_power, rfi_power, $
                              rfi_margin, rfi_3dB_bandwidth, chamber_type

      save_spectrum, Spectrum, FreqAxis, write_spectrum_to_csv, filename_base, Frequency, SamplingFrequency, AcqBW

      cleanup, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, AntennaGain, Glna, CableLoss, S

   endif

endfor


End


