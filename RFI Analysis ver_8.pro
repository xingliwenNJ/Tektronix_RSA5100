
;**************************************************************************************************
; Author:       R. Lord
;		  :	    Braam Otto
; Organisation: SKA SA
; Ver:          8 (Centre Frequency List Mode)
;**************************************************************************************************


;---Input
; DataFilename_base        = '2019_01_22_Digitiser_UHF_026_Test 1_'
DataFilename_base = 'Sweep_Test_Ver111_'

;---Frequency List
;openr, lun, '.\freq_list.txt', /get_lun
openr, lun, 'C:\Users\RTA\Desktop\Ver_1_1_1\freq_list\freq_list.txt', /get_lun
n=0
line = ''
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  n=n+1
ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun

;openr, lun, '.\freq_list.txt', /get_lun
openr, lun, 'C:\Users\RTA\Desktop\Ver_1_1_1\freq_list\freq_list.txt', /get_lun
array_list = STRARR(n)
line = ''
n=0
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array_list[n] = DataFilename_base+line+'MHz.tiq'
  n=n+1
ENDWHILE
; Close the file and free the file unit
FREE_LUN, lun

print, array_list

DataFilename_list		 = array_list


DataFilename_start_num   = 1
DataFilename_end_num     = 1
Stitch_start_num         = 1
Stitch_end_num           = 1
RBW_list                 = [0.5, .1]
subtitle                 = 'Digitiser acceptance test 1 - BRP reverb chamber'                      ;Description of measurement
polarisation             = 0                                            ;0=V polarisation, 1=H polarisation
chamber_type             = 1                                            ;0=Anechoic Chamber, 1=Reverb Chamber, 2=Current Probe


;---Execution
perform_batch            = 1                                            ;0=no, 1=yes (process individual spectra)
perform_stitch           = 0                                            ;0=no, 1=yes (stitch spectra together)
view_results             = 0                                            ;0=no, 1=yes


;---Calibration Data
AntGainFilenameNumber    = 1                                            ;Refer to numbers below
LnaGainFilenameNumber    = 3                                            ;Refer to numbers below
CableLossFilenameNumber  = 17                                           ;Refer to numbers below
CCF_FilenameNumber       = 2                                            ;Refer to numbers below
CProbeFilenameNumber     = 2                                            ;Refer to numbers below
range                    = 1.0                                          ;Distance from DUT to receiving antenna


;---Threshold Level (red and blue dotted lines)
CulpritCategoryNumber    = 1                                            ;Refer to numbers below
culprit_distance         = 1.0                                          ;Distance between Antenna Focus and DUT [m]. Only for d <= 100m (Culprit Category #1 and #9)
enclosure_attenuation    = 0.0                                          ;Additional attenuation provided by enclosure [dB]
show_spectral_threshold  = 1                                            ;Show spectral line threshold level
show_continuum_threshold = 0                                            ;Show continuum threshold level


;---Noise Floor (green line)
show_noisefloor          = 0                                            ;Show noise floor level


AntGainFilenameList      = ['Antenna_MESA_GLPDA'                        , $ ;#1
                            'Antenna_MESA_KLPDA1'                       , $ ;#2
                            'Antenna_MESA_KLPDA3'                       , $ ;#3
                            'Antenna_Houwteq_EMCO3115_FRH_Small'        , $ ;#4
                            'Antenna_Houwteq_EMCO3106B_FRH_Large'       , $ ;#5
                            'Antenna_Houwteq_EMCO3141_BiConiLog_FLPDA'  , $ ;#6
                            'Antenna_RFI_Trailor_LPDA'                  , $ ;#7
                            'No_Antenna'                                , $ ;#8 (Use this option for direct injection into RTSA and for reverberation chamber)

                            'Antenna_Houwteq_BICON_HUF_Z2_16July2018'   , $ ;#9
                            'Antenna_Houwteq_HORN_EMCO3106_16July2018'  , $ ;#10
                            'Antenna_Houwteq_HORN_EMCO3115_16July2018']     ;#11

LnaGainFilenameList      = ['LNA_1'                                     , $ ;#1
                            'LNA_2'                                     , $ ;#2
                            'LNA_3'                                     , $ ;#3
                            'LNA_1_plus_2'                              , $ ;#4
                            'No_LNA'                                    , $ ;#5 (Use this option for passive antenna)

                            'Miteq_AFS4_Interp'                         , $ ;#6
                            'Miteq_JS4_Interp']                             ;#7

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
                            'CableLoss_RFI_Trailor'                     , $ ;#19
                            'No_Cable'                                  , $ ;#20

                            'CableLoss_Houwteq_1a_1b_16July2018']           ;#21

CCF_FilenameList         = ['ACF'                                       , $ ;#1
                            'ACF_03_06_15']                                 ;#2

CProbeFilenameList       = ['1 GHz Current Probe Transfer Impedance'    , $ ;#1
                            '3 GHz Current Probe Transfer Impedance']       ;#2

CulpritCategoryList      = ['d <= 100m for f > 500MHz'   , $ ;#1 d=culprit_distance for f>500 MHz, d=500m and h_tx=10m for f<500 MHz
                            'd = 1km'                    , $ ;#2 d=1km   and h_tx=1.5m for all f
                            'd = 10km'                   , $ ;#3 d=10km  and h_tx=1.5m for all f
                            'Core Weather Station'       , $ ;#4 d=100m  for f>350 MHz, d=2140m and h_tx=10m  for f<350 MHz
                            'Core Wind Station'          , $ ;#5 d=100m  for f>350 MHz, d=4760m and h_tx=10m  for f<350 MHz
                            'Meysdam Construction Camp'  , $ ;#6 d=1200m for f>350 MHz, d=3000m and h_tx=1.5m for f<350 MHz
                            'Losberg Farm'               , $ ;#7 d=1000m for f>350 MHz, d=4820m and h_tx=1.5m for f<350 MHz
                            'Losberg Site Complex'       , $ ;#8 d=2350m for f>350 MHz, d=3740m and h_tx=1.5m for f<350 MHz
                            'd <= 100m for all f']           ;#9 d=culprit_distance for all frequencies


;---Other
fraction_to_process      = 1.0                                          ;Fraction of data to process (1.0 = 100%)
num_peaks                = 100                                          ;Max number of RFI peaks to identify
threshold_margin         = -3.0                                         ;Margin within Spectral Line threshold level to identify RFI spike [dB]
plant_signals            = 0                                            ;0=no, 1=yes
freq_signal_at_threshold = 519.4e6
freq_signal_at_Xsigma    = 519.6e6
Xsigma                   = 7.0
remove_datafile          = 0                                            ;0=no, 1=yes (delete .tiq raw data files)
purge_spectrum_to_dat    = 0                                            ;0=no, 1=yes (delete binary spectra .dat files)
input_dir                = 'C:\RFI Raw Data\data\'
output_dir               = '.\'
bitmap_dir               = 'C:\RFI Archive\Equipment_Database\Bitmaps\'
calibration_dir_ant      = 'C:\RFI Archive\Equipment_Database\Passive_Antennas\'
calibration_dir_lna      = 'C:\RFI Archive\Equipment_Database\LNAs\'
calibration_dir_cable    = 'C:\RFI Archive\Equipment_Database\Cables\'
calibration_dir_ccf      = 'C:\RFI Archive\Equipment_Database\Reverb_Chamber\'
calibration_dir_cp       = 'C:\RFI Archive\Equipment_Database\Current_Probes\'


;---RFI Flag List
generate_rfi_flaglist    = 0                                            ;0=no, 1=yes (generate RFI flaglist for individual spectra)
apply_rfi_flaglist       = 0                                            ;0=no, 1=yes (apply RFI flaglist to individual spectra)
generate_flaglist_stitch = 0                                            ;0=no, 1=yes (generate RFI flaglist for stitched spectrum)
apply_flaglist_stitch    = 0                                            ;0=no, 1=yes (apply RFI flaglist to stitched spectrum)
cutoff_above_noisefloor  = 2.0                                          ;Everything above this level will be classified RFI [dB]
rfi_flaglist_filename    = 'rfi_flaglist.csv'
rfi_flaglist_folder      = 'RFI_Flaglist\'


;---Reverberation Chamber
antenna_efficiency       = 0.8


;---Plot parameters
yaxis_bottom_add         = 0.0                                          ;Raises bottom of y-axis by this amount
yaxis_top_add            = 0.0                                          ;Drops top of y-axis by this amount
xaxis_centre             = 500e6                                        ;Zoom into spectrum at this centre frequency [Hz]
xaxis_span               = 25e9                                         ;Zoom into spectrum with this span (set to large value to see everything) [Hz]
xsize                    = 1904                                         ;Horizontal size of plot window
ysize                    = 1000                                         ;Vertical size of plot window
log_xaxis                = 0                                            ;0=no, 1=yes, plot data with logarithmic x-axis


;**************************************************************************************************
;**************************************************************************************************


;---Compile routines
Resolve_Routine, 'rfi_library_ver_5'
Resolve_Routine, 'strsplit', /is_function
Resolve_Routine, 'interpol', /is_function
headerlines = 2
maxbins = 500000000L


if (purge_spectrum_to_dat eq 1) then begin
   purge_spectrum_dat, DataFilename_base, DataFilename_start_num, DataFilename_end_num, output_dir
   retall
endif


if (perform_batch eq 1) then begin

   progress_indicator, progress_widget, title='Batch Processing Individual Spectra'
   show_plot = 0

   if (remove_datafile eq 1) then begin
      result = Dialog_Message('Delete raw data?', /question, /default_no, title='Warning')
      if (strlowcase(result) eq 'no') then retall
   endif

   ;---Batch processing
   ; DataFilename_list = DataFilename_base + strtrim(Indgen(DataFilename_end_num - DataFilename_start_num + 1) + DataFilename_start_num, 2) + '.tiq'
   for i=0, (n_elements(DataFilename_list) - 1) do begin
      DataFilename = DataFilename_list[i]
      wait_for_file, input_dir + DataFilename
      if (File_Test(input_dir + DataFilename)) then begin   ;Check if file exists
         for j=0, (n_elements(RBW_list) - 1) do begin
            start_time = systime(1)
            start_mem = memory(/current)
            RBW = float(RBW_list[j])
            show_message = (1 - j) > 0
            first_run    = (1 - j) > 0
            progress_indicator, progress_widget, progress = $
                                round(((i*n_elements(RBW_list) + j) * 100.0) / (n_elements(DataFilename_list) * n_elements(RBW_list)))
            print, strtrim(round(((i*n_elements(RBW_list) + j) * 100.0) / (n_elements(DataFilename_list) * n_elements(RBW_list))), 2) + $
                   '% complete. Processing ' + DataFilename + ' with RBW = ' + strtrim(RBW_list[j], 2) + ' Hz'
            wait, 0.1
            process_data, DataFilename, RBW, polarisation, AntGainFilenameNumber, AntGainFilenameList, $
                          LnaGainFilenameNumber, LnaGainFilenameList, CableLossFilenameNumber, $
                          CableLossFilenameList, CCF_FilenameNumber, CCF_FilenameList, CProbeFilenameNumber, $
                          CProbeFilenameList, AntGainFilename, LnaGainFilename, CableLossFilename, $
                          CCF_filename, CProbeFilename, CulpritCategoryNumber, range, culprit_distance, $
                          enclosure_attenuation, fraction_to_process, num_peaks, threshold_margin, headerlines, $
                          input_dir, output_dir, calibration_dir_ant, calibration_dir_lna, calibration_dir_cable, $
                          calibration_dir_ccf, calibration_dir_cp, plant_signals, freq_signal_at_threshold, $
                          freq_signal_at_Xsigma, Xsigma, Signal, Spectrum, FreqAxis, noise_floor_E, $
                          SpectralE_spec, SpectralE_interp, ContinuumE_interp, AntennaGain, Glna, CableLoss, $
                          SamplingFrequency, fft_length, num_ffts, Frequency, AcqBW, date, Preamp, rfi_SNR, $
                          RFAttenuation, rfi_frequency, rfi_above_spec, rfi_plus_noise_power, rfi_power, $
                          rfi_margin, rfi_3dB_bandwidth, apply_rfi_flaglist, generate_rfi_flaglist, $
                          cutoff_above_noisefloor, rfi_flaglist_filename, rfi_flaglist_folder, chamber_type, $
                          antenna_efficiency, CCF, CurrentProbe, progress_widget, show_message

            save_result,  Signal, Spectrum, FreqAxis, noise_floor_E, SpectralE_spec, SpectralE_interp, $
                          ContinuumE_interp, AntennaGain, Glna, CableLoss, DataFilename, AntGainFilename, $
                          LnaGainFilename, CableLossFilename, subtitle, yaxis_bottom_add, $
                          yaxis_top_add, xaxis_centre, xaxis_span, SamplingFrequency, fft_length, num_ffts, $
                          Frequency, AcqBW, output_dir, date, polarisation,  Preamp, RFAttenuation, range, $
                          start_time, start_mem, fraction_to_process, num_peaks, plant_signals, $
                          freq_signal_at_threshold, freq_signal_at_Xsigma, threshold_margin, rfi_frequency, $
                          rfi_above_spec, rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, $
                          rfi_SNR, show_plot, culprit_distance, enclosure_attenuation, apply_rfi_flaglist, $
                          generate_rfi_flaglist, cutoff_above_noisefloor, rfi_flaglist_folder, $
                          rfi_flaglist_filename, chamber_type, antenna_efficiency, CCF_filename, CCF, $
                          CurrentProbe, CProbeFilename, show_spectral_threshold, show_continuum_threshold, $
                          show_noisefloor, xsize, ysize, log_xaxis=log_xaxis, first_run=first_run

            cleanup,      Signal, Spectrum, FreqAxis, SpectralE_spec, SpectralE_interp, ContinuumE_interp, $
                          noise_floor_E, AntennaGain, Glna, CableLoss, CCF
         endfor

         if (remove_datafile eq 1) then begin
            label = ['Warning:', $
                     'Deleting raw data file ' + DataFilename, $
                     '']
            info_widget, label, 2.0
            File_Delete, input_dir + DataFilename
         endif

      endif else begin
         label = ['Warning:', $
                  DataFilename + ' could not be opened!', $
                  '']
         info_widget, label, 2.0
      endelse
   endfor
   progress_indicator, progress_widget, /kill
   print, '100% complete with batch processing'
   print, ''

endif


if (perform_stitch eq 1) then begin

   progress_indicator, progress_widget, title='Stitch Processing'
   show_plot = 0

   ;---Error checking
   if (Stitch_end_num lt Stitch_start_num) then begin
      result = Dialog_Message(['Stitch_end_num must be >= than Stitch_start_num!', $
               '', 'Script will terminate.'], title='Error')
      retall
   endif

   ;---Get Filenames
   get_filenames, AntGainFilenameNumber, AntGainFilenameList, AntGainFilename, LnaGainFilenameNumber, $
                  LnaGainFilenameList, LnaGainFilename, CableLossFilenameNumber, CableLossFilenameList, $
                  CableLossFilename, CCF_FilenameNumber, CCF_FilenameList, CCF_filename, $
                  CProbeFilenameNumber, CProbeFilenameList, CProbeFilename

   ;---Start loop
   for j=0, (n_elements(RBW_list) - 1) do begin
      start_time = systime(1)
      start_mem = memory(/current)
      RBW = float(RBW_list[j])
      show_message = (1 - j) > 0
      progress_indicator, progress_widget, progress = round(((float(j) + 1.0) / n_elements(RBW_list)) * 100)
      print, 'Stitching spectra with RBW = ' + strtrim(RBW_list[j], 2) + ' Hz'
      wait, 0.1

      Folder_list = DataFilename_base + strtrim(Indgen(Stitch_end_num - Stitch_start_num + 1) + Stitch_start_num, 2)

      get_stitch_spectra_filenames, output_dir, Folder_list, RBW, stitch_filenames, signal_filenames

      is_array = size(stitch_filenames)

      if (is_array[0] eq 1) then begin

         get_stitch_spectra_parameters, Folder_list, stitch_filenames, NumberSamples, CentreFrequency, SamplingFrequency, AcqBW

         get_stitch_spectra_Spectrum, Spectrum, FreqAxis, stitch_filenames, SamplingFrequency, NumberSamples, CentreFrequency, $
                                      AcqBW, shrinkfactor, maxbins

         get_stitch_td_Signal, Signal, signal_filenames, time_window_size, time_resolution, time_windows, SamplingFrequency_TD

         get_stitch_spectra_values, FreqAxis, CentreFrequency, SamplingFrequency, NumberSamples, fft_length, AcqBW, $
                                    Frequency, shrinkfactor

         read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation, AntennaGain, $
                            chamber_type, show_message

         read_ccf, calibration_dir_ccf, FreqAxis, chamber_type, CCF_filename, 0, CCF

         read_current_probe, calibration_dir_cp, CProbeFilename, 10, FreqAxis, chamber_type, CurrentProbe

         get_noise_floor, Spectrum, FreqAxis, SamplingFrequency, fft_length, noise_floor_E, AcqBW

         get_rfi_flaglist, Spectrum, FreqAxis, noise_floor_E, SamplingFrequency, generate_flaglist_stitch, cutoff_above_noisefloor, $
                           rfi_flaglist_filename, rfi_flaglist_folder

         get_RFI_thresholds, FreqAxis, range, CulpritCategoryNumber, culprit_distance, enclosure_attenuation, ContinuumE_spec, $
                             SpectralE_spec, chamber_type

         SpectralE_interp = 10.0 * alog10((10.0 ^ (SpectralE_spec / 10.0)) + (10.0 ^ (noise_floor_E / 10.0)))
         ContinuumE_interp = 10.0 * alog10((10.0 ^ (ContinuumE_spec / 10.0)) + (10.0 ^ (noise_floor_E / 10.0)))

         remove_rfi_flags, Spectrum, FreqAxis, noise_floor_E, apply_flaglist_stitch, rfi_flaglist_filename, rfi_flaglist_folder

         get_stitch_spectra_plot, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, Folder_list, $
                                  SamplingFrequency, fft_length, Frequency, AcqBW, subtitle, output_dir, DataFilename, $
                                  DataFilename_base, Stitch_start_num, Stitch_end_num, AntennaGain, CurrentProbe, $
                                  filename_base, show_plot, show_spectral_threshold, show_continuum_threshold, $
                                  show_noisefloor, xsize, ysize, FolderName, chamber_type, antenna_efficiency, CCF, range, $
                                  log_xaxis=log_xaxis

         get_stitch_td_Signal_plot, Signal, SamplingFrequency_TD, NumberSamples, time_windows, subtitle, filename_base, $
                                    Folder_list, DataFilename_base, Stitch_start_num, Stitch_end_num, output_dir, $
                                    show_plot, xsize, ysize
         Signal = 0

         if (chamber_type ne 2) then begin
            find_rfi_spikes, Spectrum, FreqAxis, noise_floor_E, SpectralE_spec, SpectralE_interp, AntennaGain, num_peaks, $
                             threshold_margin, SamplingFrequency, fft_length, 0, rfi_frequency, rfi_above_spec, $
                             rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, rfi_SNR, Frequency, AcqBW
         endif

         stitch_spectra_summary, Spectrum, FreqAxis, noise_floor_E, SpectralE_spec, SpectralE_interp, filename_base, Folder_list, $
                                 AntGainFilename, AntennaGain, polarisation, Frequency, AcqBW, SamplingFrequency, $
                                 fft_length, subtitle, range, culprit_distance, enclosure_attenuation, start_time, $
                                 start_mem, num_peaks, threshold_margin, rfi_frequency, rfi_above_spec, rfi_plus_noise_power, $
                                 rfi_power, rfi_margin, rfi_3dB_bandwidth, chamber_type, antenna_efficiency, CCF, LnaGainFilename, $
                                 CableLossFilename, CCF_filename, CProbeFilename, CurrentProbe
         SpectralE_spec = 0
         SpectralE_interp = 0
         ContinuumE_interp = 0
         noise_floor_E = 0

         save_spectrum, Spectrum, FreqAxis, filename_base, Frequency, SamplingFrequency, AcqBW
         Spectrum = 0

         if (j eq 0) then begin
            save_stitch_cal_plots, FreqAxis, AntennaGain, AntGainFilename, polarisation, Glna, LnaGainFilename, CableLoss, $
                                   CableLossFilename, output_dir, chamber_type, CCF_filename, CCF, xsize, ysize, FolderName, $
                                   headerlines, calibration_dir_lna, calibration_dir_cable, calibration_dir_cp, CProbeFilename, $
                                   CurrentProbe
         endif
         FreqAxis = 0

      endif

   endfor
   progress_indicator, progress_widget, /kill
   print, '100% complete with stitch processing'
   print, ''

endif

if (view_results eq 1) then rfi_viewer, DataFilename_base, Stitch_start_num, Stitch_end_num, RBW_list, output_dir, bitmap_dir


End


