
;**************************************************************************************************
; Author:       R. Lord
; Organisation: SKA SA
; Ver:          1.0
;
; Good reference: http://www.bitweenie.com/listings/power-spectrum-matlab/
; Good reference: http://www.ee.umd.edu/~tretter/commlab/c6713slides/ch4.pdf (from Justin)
;
;**************************************************************************************************



;**************************************************************************************************
; Get Filenames
;**************************************************************************************************

Pro get_filenames, AntGainFilenameNumber, AntGainFilenameList, AntGainFilename, $
                   LnaGainFilenameNumber, LnaGainFilenameList, LnaGainFilename, $
                   CableLossFilenameNumber, CableLossFilenameList, CableLossFilename


if ((AntGainFilenameNumber lt 1) OR (AntGainFilenameNumber gt n_elements(AntGainFilenameList))) then begin
   print, 'Error: Number supplied for Antenna Gain filename must be >=1 and <=' + strtrim(n_elements(AntGainFilenameList), 2)
   retall
endif
if ((LnaGainFilenameNumber lt 1) OR (LnaGainFilenameNumber gt n_elements(LnaGainFilenameList))) then begin
   print, 'Error: Number supplied for LNA Gain filename must be >=1 and <=' + strtrim(n_elements(LnaGainFilenameList), 2)
   retall
endif
if ((CableLossFilenameNumber lt 1) OR (CableLossFilenameNumber gt n_elements(CableLossFilenameList))) then begin
   print, 'Error: Number supplied for Cable Loss filename must be >=1 and <=' + strtrim(n_elements(CableLossFilenameList), 2)
   retall
endif

AntGainFilename   = AntGainFilenameList[AntGainFilenameNumber - 1] + '.csv'
LnaGainFilename   = LnaGainFilenameList[LnaGainFilenameNumber - 1] + '.csv'
CableLossFilename = CableLossFilenameList[CableLossFilenameNumber - 1] + '.csv'


End


;**************************************************************************************************
; Check for errors
;**************************************************************************************************

Pro error_checks, fraction_to_process, input_dir, calibration_dir_ant, calibration_dir_lna, $
                  calibration_dir_cable, DataFilename, AntGainFilename, LnaGainFilename, $
                  CableLossFilename, culprit_distance, enclosure_attenuation, chamber_type, CCF_filename


if (NOT File_Test(input_dir + DataFilename)) then begin
   print, 'Error: Cannot open data file ' + input_dir + DataFilename
   retall
endif
if (NOT File_Test(calibration_dir_ant + AntGainFilename)) then begin
   print, 'Error: Cannot open Antenna Gain file ' + calibration_dir_ant + AntGainFilename
   retall
endif
if (NOT File_Test(calibration_dir_lna + LnaGainFilename)) then begin
   print, 'Error: Cannot open LNA Gain file ' + calibration_dir_lna + LnaGainFilename
   retall
endif
if (NOT File_Test(calibration_dir_cable + CableLossFilename)) then begin
   print, 'Error: Cannot open Cable Loss file ' + calibration_dir_cable + CableLossFilename
   retall
endif
if ((chamber_type eq 1) AND (NOT File_Test(CCF_filename))) then begin
   print, 'Error: Cannot open CCF file ' + CCF_filename
   retall
endif
if ((fraction_to_process gt 1.0) OR (fraction_to_process le 0.0)) then begin
   print, 'Error: fraction_to_process must be between 0.0 and 1.0 !'
   retall
endif
if ((culprit_distance le 0.0) OR (culprit_distance gt 100.0)) then begin
   print, 'Error: culprit_distance (' + strtrim(culprit_distance, 2) + ') must be >0m and <= 100m'
   retall
endif
if (enclosure_attenuation lt 0.0) then begin
   print, 'Error: enclosure_attenuation (' + strtrim(enclosure_attenuation, 2) + ') must be >=0 dB'
   retall
endif

End


;**************************************************************************************************
; Read parameters
;**************************************************************************************************

Pro read_params, params_filename, DataFilename_base, DataFilename_start_num, DataFilename_end_num, $
                 RBW_list, subtitle, polarisation, chamber_type, AntGainFilenameNumber, $
                 LnaGainFilenameNumber, CableLossFilenameNumber, range, culprit_distance, $
                 enclosure_attenuation, antenna_efficiency, CCF_filename


if (NOT File_Test(params_filename)) then begin
   print, 'Error: Cannot open parameter file ' + params_filename
   retall
endif

openr, lun, params_filename, /get_lun
temp = ''
readf, lun, temp

while (strpos(temp, 'DataFilename_base') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
DataFilename_base = temp[1]
if (strpos(DataFilename_base, "'") ge 0) then begin
   DataFilename_base = strsplit(DataFilename_base, "'", /extract)
   DataFilename_base = DataFilename_base[1]
endif
if (strpos(DataFilename_base, '"') ge 0) then begin
   DataFilename_base = strsplit(DataFilename_base, '"', /extract)
   DataFilename_base = DataFilename_base[1]
endif

temp = ''
while (strpos(temp, 'DataFilename_start_num') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
DataFilename_start_num = fix(temp[1])

temp = ''
while (strpos(temp, 'DataFilename_end_num') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
DataFilename_end_num = fix(temp[1])

temp = ''
while (strpos(temp, 'RBW_list') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
RBW_list = temp[1]
RBW_list = float(strsplit(strmid(RBW_list, strpos(RBW_list, '[')+1, strpos(RBW_list, ']')-1), ',', /extract))

temp = ''
while (strpos(temp, 'subtitle') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
subtitle = temp[1]
if (strpos(subtitle, "'") ge 0) then begin
   subtitle = strsplit(subtitle, "'", /extract)
   subtitle = subtitle[1]
endif
if (strpos(subtitle, '"') ge 0) then begin
   subtitle = strsplit(subtitle, '"', /extract)
   subtitle = subtitle[1]
endif

temp = ''
while (strpos(temp, 'polarisation') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
polarisation = fix(temp[1])

temp = ''
while (strpos(temp, 'chamber_type') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
chamber_type = fix(temp[1])

temp = ''
while (strpos(temp, 'AntGainFilenameNumber') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
AntGainFilenameNumber = fix(temp[1])

temp = ''
while (strpos(temp, 'LnaGainFilenameNumber') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
LnaGainFilenameNumber = fix(temp[1])

temp = ''
while (strpos(temp, 'CableLossFilenameNumber') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
CableLossFilenameNumber = fix(temp[1])

temp = ''
while (strpos(temp, 'range') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
range = float(temp[1])

temp = ''
while (strpos(temp, 'culprit_distance') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
culprit_distance = float(temp[1])

temp = ''
while (strpos(temp, 'enclosure_attenuation') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
enclosure_attenuation = float(temp[1])

temp = ''
while (strpos(temp, 'antenna_efficiency') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
antenna_efficiency = float(temp[1])

temp = ''
while (strpos(temp, 'CCF_filename') eq -1) do readf, lun, temp
temp = strsplit(temp, '=', /extract)
CCF_filename = temp[1]
if (strpos(CCF_filename, "'") ge 0) then begin
   CCF_filename = strsplit(CCF_filename, "'", /extract)
   CCF_filename = CCF_filename[1]
endif
if (strpos(CCF_filename, '"') ge 0) then begin
   CCF_filename = strsplit(CCF_filename, '"', /extract)
   CCF_filename = CCF_filename[1]
endif

free_lun, lun

End


;**************************************************************************************************
; Read header
;**************************************************************************************************

Pro read_header, lun, offset, input_dir, DataFilename, SamplingFrequency, NumberSamples, Scaling, $
                 date, Frequency, AcqBW, RFAttenuation, Preamp


openr, lun, input_dir + DataFilename, /get_lun

temp = ''
readf, lun, temp
offset = long(strmid(temp, 18, 9))

temp = ''
while (strpos(temp, 'SamplingFrequency') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
SamplingFrequency = float(temp[0])

temp = ''
while (strpos(temp, 'NumberSamples') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
NumberSamples = long(temp[0])

temp = ''
while (strpos(temp, 'Scaling') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
Scaling = float(temp[0])

temp = ''
while (strpos(temp, 'DateTime') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
date = strmid(temp[0], 0, 10)

temp = ''
while (strpos(temp, 'Frequency') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
Frequency = double(temp[0])

temp = ''
while (strpos(temp, 'AcquisitionBandwidth') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
AcqBW = float(temp[0])

temp = ''
while (strpos(temp, 'RFAttenuation') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
RFAttenuation = float(temp[0])

temp = ''
while (strpos(temp, 'Preamp') eq -1) do readf, lun, temp
temp = strsplit(temp, '>', /extract)
temp = strsplit(temp[1], '<', /extract)
Preamp = temp[0]

free_lun, lun

End


;**************************************************************************************************
; Get FFT Length
;**************************************************************************************************

Pro get_fft_length, FFT_power_of_2, SamplingFrequency, RBW, NumberSamples, fraction_to_process, $
                    fft_length, num_ffts


if (FFT_power_of_2 eq 1) then begin
   fft_length = 2L ^ round(alog(SamplingFrequency / RBW) / alog(2))  ;power of 2
endif else begin
   fft_length = round((SamplingFrequency / RBW) / 2.0) * 2L          ;even number
endelse

if (fft_length lt 16) then begin
   RBW = SamplingFrequency / 16.0
   print, 'Warning: RBW too large! New RBW = ' + strtrim(RBW, 2) + ' Hz'
   fft_length = round((SamplingFrequency / RBW) / 2.0) * 2L          ;even number
endif

num_ffts = floor((NumberSamples * fraction_to_process) / fft_length)
if (num_ffts lt 1) then begin
   print, 'Warning: Too little data to process with desired FFT length of ' + strtrim(fft_length, 2)
   fft_length = floor(NumberSamples * fraction_to_process)
   num_ffts = 1
   print, 'New FFT length = ' + strtrim(fft_length, 2)
   print, 'Number of FFTs = ' + strtrim(num_ffts, 2)
endif

End


;**************************************************************************************************
; Get Power Sepctrum in dBm
;**************************************************************************************************

Pro get_power_spectrum_dBm, Spectrum, S, input_dir, DataFilename, offset, fft_length, num_ffts, $
                            Scaling, transient_RFI_detection, peak_hold


Spectrum = fltarr(fft_length)
temp = lonarr(2, fft_length)
if ((transient_RFI_detection eq 1) AND (num_ffts ge 10)) then begin
   M = fltarr(fft_length)
   S = fltarr(fft_length)
endif

openr, lun, input_dir + DataFilename, /get_lun
point_lun, lun, offset

if (peak_hold eq 1) then begin
   for i=0L, (num_ffts - 1L) do begin
      readu, lun, temp
      Spectrum_temp = abs(FFT(complex(temp[0, *], temp[1, *]) * Scaling)) ^ 2.0
      index = where(Spectrum_temp gt Spectrum, count)
      if (count gt 0) then Spectrum[index] = Spectrum_temp[index]
   endfor
   index = 0
   Spectrum_temp = 0
   scale_factor = (1000.0 / 50.0)
endif else begin
   for i=0L, (num_ffts - 1L) do begin
      readu, lun, temp
      if ((transient_RFI_detection eq 1) AND (num_ffts ge 10)) then begin
         power_FFT = abs(FFT(complex(temp[0, *], temp[1, *]) * Scaling)) ^ 2.0
         Spectrum = Temporary(Spectrum) + power_FFT
         power_FFT = Temporary(power_FFT) * 1e12
         Mprev = M
         M = M + ((power_FFT - M) / (i + 1.0))
         S = Temporary(S) + ((power_FFT - Mprev) * (power_FFT - M))
      endif else begin
         Spectrum = Temporary(Spectrum) + (abs(FFT(complex(temp[0, *], temp[1, *]) * Scaling)) ^ 2.0)
      endelse
   endfor
   if ((transient_RFI_detection eq 1) AND (num_ffts ge 10)) then begin
      S = 10.0 * alog10(shift(sqrt(Temporary(S) / (num_ffts - 1L)) / M, fft_length / 2L))
      power_FFT = 0
      Mprev = 0
      M = 0
   endif
   scale_factor = (1000.0 / 50.0) / num_ffts
endelse

Spectrum = 10.0 * alog10(shift(Temporary(Spectrum) * scale_factor, fft_length / 2L))   ;[dBm]
free_lun, lun
temp = 0

End


;**************************************************************************************************
; Convert Power to Electric Field Strength
;**************************************************************************************************

Pro get_efield_strength, Spectrum, FreqAxis, AntennaGain, Glna, CableLoss, CCF, chamber_type, $
                         antenna_efficiency, range


;---Constants
Z     = 119.9169832 * !pi
light = 299792458.0


if (chamber_type eq 0) then begin ;Anechoic Chamber

   Spectrum = Temporary(Spectrum) + (20.0 * alog10(FreqAxis)) - AntennaGain - Glna + CableLoss + $
              (10.0 * alog10(4.0 * !pi * Z / (light ^ 2.0))) + 90.0                            ;[dBuV/m]

endif else begin                  ;Reverberation Chamber

   Spectrum = Temporary(Spectrum) - Glna + CableLoss - (10.0 * alog10(antenna_efficiency)) - $
              CCF + (10.0 * alog10(Z / (4.0 * !pi * (range ^ 2.0)))) + 90.0   ;[dBuV/m]

endelse


End


;**************************************************************************************************
; Find Noise Floor
;**************************************************************************************************

Pro get_noise_floor, Spectrum, FreqAxis, SamplingFrequency, fft_length, noise_floor_E, AcqBW


fbinsize      = SamplingFrequency / fft_length
numpoints_s1  = 10 < (round(((SamplingFrequency - AcqBW) / 2.0) / fbinsize) / 5)
numpoints_s2  = round(AcqBW / 1e6) > 5
numpoints_s3  = 10 < (round(((SamplingFrequency - AcqBW) / 2.0) / fbinsize) / 5)

if (numpoints_s1 ge 2) then begin

   FreqAxis_temp = dblarr(numpoints_s1 + numpoints_s2 + numpoints_s3)
   Spectrum_temp = fltarr(numpoints_s1 + numpoints_s2 + numpoints_s3)

   ;---Section 1
   m_length = round((((SamplingFrequency - AcqBW) / 2.0) / fbinsize) / numpoints_s1)
   for i=0, (numpoints_s1 - 1) do begin
      FreqAxis_temp[i] = FreqAxis[((i * m_length) + (m_length / 2)) < (fft_length - 1)]
      Spectrum_temp[i] = median(Spectrum[i * m_length : (((i + 1) * m_length) - 1) < (fft_length - 1)])
   endfor

   ;---Section 2
   m_length = round((AcqBW / fbinsize) / numpoints_s2)
   index_start = round(((SamplingFrequency - AcqBW) / 2.0) / fbinsize)
   for i=0, (numpoints_s2 - 1) do begin
      index_nan = where(~Finite(Spectrum[index_start + (i * m_length) : (index_start + (((i + 1) * m_length) - 1)) < (fft_length - 1)]), count)
      if (count eq 0) then begin
         FreqAxis_temp[numpoints_s1 + i] = FreqAxis[(index_start + (i * m_length) + (m_length / 2)) < (fft_length - 1)]
         Spectrum_temp[numpoints_s1 + i] = median(Spectrum[index_start + (i * m_length) : (index_start + (((i + 1) * m_length) - 1)) < (fft_length - 1)])
      endif
   endfor

   ;---Section 3
   m_length = round((((SamplingFrequency - AcqBW) / 2.0) / fbinsize) / numpoints_s3)
   index_start = round((((SamplingFrequency - AcqBW) / 2.0) + AcqBW) / fbinsize)
   for i=0, (numpoints_s3 - 1) do begin
      FreqAxis_temp[numpoints_s1 + numpoints_s2 + i] = FreqAxis[(index_start + (i * m_length) + (m_length / 2)) < (fft_length - 1)]
      Spectrum_temp[numpoints_s1 + numpoints_s2 + i] = median(Spectrum[index_start + (i * m_length) : (index_start + (((i + 1) * m_length) - 1)) < (fft_length - 1)])
   endfor

   index = where(FreqAxis_temp eq 0.0, count, complement=nonzero_index)
   if (count gt 0) then begin
      FreqAxis_temp = FreqAxis_temp[nonzero_index]
      Spectrum_temp = Spectrum_temp[nonzero_index]
   endif

   noise_floor_E = interpol(Spectrum_temp, FreqAxis_temp, FreqAxis, /spline)           ;[dBuV/m]

endif else begin
   noise_floor_E = fltarr(fft_length)
   noise_floor_E[*] = median(Spectrum)
endelse


End


;**************************************************************************************************
; Get RFI Flag List
;**************************************************************************************************

Pro get_rfi_flaglist, Spectrum, FreqAxis, noise_floor_E, SamplingFrequency, generate_rfi_flaglist, $
                      cutoff_above_noisefloor, rfi_flaglist_filename, rfi_flaglist_folder

if (generate_rfi_flaglist eq 1) then begin
   index = where(Spectrum ge (noise_floor_E + cutoff_above_noisefloor), count)
   if (count gt 0) then begin
      fbinsize = SamplingFrequency / n_elements(Spectrum)
      start_f = [double(FreqAxis[index[0]])]
      end_f = [0d]
      for i=1L, (count - 1L) do begin
         if (index[i] gt (index[i-1] + 1)) then begin   ;new spike
            start_f = [start_f, FreqAxis[index[i]]]
            end_f = [end_f, FreqAxis[index[i-1]]]
         endif
      endfor
      end_f = [end_f, FreqAxis[index[count - 1]]]
      end_f = end_f[1 : n_elements(end_f) - 1]

      file_mkdir, rfi_flaglist_folder
      openw, lun, rfi_flaglist_folder + rfi_flaglist_filename, /get_lun
      printf, lun, 'Centre Frequency [Hz], Bandwidth [Hz]'
      if (n_elements(start_f) eq 1) then begin
         printf, lun, strtrim(string((start_f + end_f) / 2d, format='(F15.1)'), 2) + ', ' + $
                      strtrim(string((end_f - start_f) + fbinsize, format='(F15.1)'), 2)
      endif else begin
         printf, lun, transpose(strtrim(string((start_f + end_f) / 2d, format='(F15.1)'), 2) + ', ' + $
                      strtrim(string((end_f - start_f) + fbinsize, format='(F15.1)'), 2))
      endelse
      free_lun, lun
   endif else begin
      print, 'Warning: No RFI spikes detected to write to RFI Flag List!'
   endelse
endif

End


;**************************************************************************************************
; Remove DC-Offset Spike
;**************************************************************************************************

Pro remove_DC, Spectrum, noise_floor_E, remove_dc_offset, fft_length

if (remove_dc_offset eq 1) then begin
   Spectrum[fft_length / 2L] = noise_floor_E[fft_length / 2L]
endif

End


;**************************************************************************************************
; Apply RFI Flag List
;**************************************************************************************************

Pro remove_rfi_flags, Spectrum, FreqAxis, noise_floor_E, apply_rfi_flaglist, rfi_flaglist_filename, $
                      rfi_flaglist_folder


if (apply_rfi_flaglist eq 1) then begin

   if (NOT File_Test(rfi_flaglist_folder + rfi_flaglist_filename)) then begin
      print, 'Warning: Could not open RFI Flaglist file ' + rfi_flaglist_folder + rfi_flaglist_filename
   endif else begin

      headerlines = 1
      numsamples = File_Lines(rfi_flaglist_folder + rfi_flaglist_filename) - headerlines
      centre_f = fltarr(numsamples)
      bandwidth = fltarr(numsamples)

      ;---Read header
      openr, lun, rfi_flaglist_folder + rfi_flaglist_filename, /get_lun
      if (headerlines gt 0) then begin
         header = strarr(headerlines)
         readf, lun, header
      endif

      ;---Read data
      row = ''
      for i=0L, (numsamples - 1L) do begin
         readf, lun, row
         temp = strsplit(row, ',', /extract)
         if (strtrim(temp[0], 2) ne '') then begin
            centre_f = double(temp[0])
            bandwidth = double(temp[1])
            start_f = centre_f - (bandwidth / 2.0)
            end_f = centre_f + (bandwidth / 2.0)

            index_start = where(FreqAxis le start_f, count)
            if (count gt 0) then begin
               index_start = index_start[count - 1]
            endif else begin
               index_start = 0
            endelse
            index_stop = where(FreqAxis ge end_f, count)
            if (count gt 0) then begin
               index_stop = index_stop[0]
            endif else begin
               index_stop = n_elements(FreqAxis) - 1
            endelse
            Spectrum[index_start : index_stop] = noise_floor_E[index_start : index_stop]
         endif
      endfor
      free_lun, lun

   endelse
endif

End


;**************************************************************************************************
; Find RFI Spikes
;**************************************************************************************************

Pro find_rfi_spikes, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, AntennaGain, $
                     num_peaks, threshold_margin, SamplingFrequency, fft_length, num_ffts, $
                     rfi_frequency, rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, $
                     rfi_SNR, Frequency, AcqBW


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


fbinsize             = SamplingFrequency / fft_length
rfi_frequency        = fltarr(num_peaks)
rfi_plus_noise_power = fltarr(num_peaks)
rfi_power            = fltarr(num_peaks)
rfi_margin           = fltarr(num_peaks)
rfi_3dB_bandwidth    = fltarr(num_peaks)
rfi_SNR              = fltarr(num_peaks)
rfi_level            = fltarr(n_elements(Spectrum), /nozero)
rfi_level[*]         = threshold_margin - 1.0
good                 = where(Finite(Spectrum))
rfi_level[good]      = Spectrum[good] - SpectralE_interp[good]


;---Only consider data in the valid Acquisition Bandwidth region
index = where((FreqAxis lt (Frequency - (AcqBW / 2.0))) OR (FreqAxis gt (Frequency + (AcqBW / 2.0))), count)
if (count gt 0L) then begin
   rfi_level[index] = threshold_margin - 1.0
endif
index = 0


for i=0, (num_peaks - 1) do begin
   peak_index = where(rfi_level eq max(rfi_level))
   index_centre = peak_index[0]
   if (rfi_level[index_centre] ge threshold_margin) then begin

      rfi_frequency[i] = FreqAxis[index_centre]                                            ;[Hz]
      rfi_plus_noise_power[i] = Spectrum[index_centre]                                     ;[dBuV/m]
      rfi_power[i] = 10.0 * alog10((10.0 ^ (rfi_plus_noise_power[i] / 10.0)) - $
                                   (10.0 ^ (noise_floor_E[index_centre] / 10.0)))          ;[dBuV/m]
      rfi_margin[i] = rfi_level[index_centre]                                              ;[dB]

      ;---Find SNR (RFI-to-sigma ratio)
      if (num_ffts gt 0) then begin
         noise_floor_dBm = noise_floor_E[index_centre] - (20.0 * alog10(FreqAxis[index_centre])) + $
                           AntennaGain[index_centre] - (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
         Tsys = (10.0 ^ (noise_floor_dBm / 10.0)) / (k * fbinsize * 1000.0)
         sigma_dBm = 10.0 * alog10(k * Tsys * SamplingFrequency * 1000.0 / (fft_length * sqrt(num_ffts)))
         rfi_SNR[i] = rfi_power[i] - (20.0 * alog10(FreqAxis[index_centre])) + AntennaGain[index_centre] - $
                      (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0 - sigma_dBm
      endif

      ;---Find 3dB bandwidth
      deltaP = 0.0
      left_count = 0L
      new_P = rfi_plus_noise_power[i]
      prev_P = new_P
      while ((deltaP lt 3.0) AND (new_P le prev_P)) do begin
         left_count += 1L
         if ((index_centre - left_count) gt 0L) then begin
            prev_P = new_P
            new_P = Spectrum[index_centre - left_count]
            deltaP = rfi_plus_noise_power[i] - new_P
         endif else break
      endwhile
      deltaP = 0.0
      right_count = 0L
      new_P = rfi_plus_noise_power[i]
      prev_P = new_P
      while ((deltaP lt 3.0) AND (new_P le prev_P)) do begin
         right_count += 1L
         if ((index_centre + right_count) lt n_elements(Spectrum)) then begin
            prev_P = new_P
            new_P = Spectrum[index_centre + right_count]
            deltaP = rfi_plus_noise_power[i] - new_P
         endif else break
      endwhile
      rfi_3dB_bandwidth[i] = (left_count + right_count) * fbinsize                         ;[Hz]

      ;---Mask out current peak
      rfi_level[(index_centre - left_count)  > 0 : $
                (index_centre + right_count) < (n_elements(Spectrum) - 1)] = threshold_margin - 1.0
   endif else break  ;all peaks have been found
endfor

End


;**************************************************************************************************
; Plant Signals
;**************************************************************************************************

Pro plant_test_signals, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, plant_signals, $
                        freq_signal_at_threshold, freq_signal_at_Xsigma, Xsigma, Frequency, AcqBW

if (plant_signals eq 1) then begin
   index1 = where(FreqAxis ge freq_signal_at_threshold, count1)
   index2 = where(FreqAxis ge freq_signal_at_Xsigma, count2)
   if (count1 gt 0) then begin
      Spectrum[index1[0]] = SpectralE_interp[index1[0]]
   endif else begin
      print, 'Error: Frequency of planted signal is out of range!'
   endelse
   if (count2 gt 0) then begin
      index = where(FreqAxis gt ((FreqAxis[index2[0]] - 0.5e6) > (Frequency - (AcqBW / 2.0))) AND $
                    FreqAxis lt ((FreqAxis[index2[0]] + 0.5e6) < (Frequency + (AcqBW / 2.0))), count)
      noise_sigma_E = (10.0 * alog10(stddev(10.0 ^ (Spectrum[index] / 10.0))))            ;[dBuV/m]
      Spectrum[index2[0]] = 10.0 * alog10((10.0 ^ (noise_floor_E[index2[0]] / 10.0)) + $
                            ((10.0 ^ (noise_sigma_E / 10.0)) * Xsigma))
   endif else begin
      print, 'Error: Frequency of planted signal is out of range!'
   endelse
endif

End


;**************************************************************************************************
; Read LNA Gain
;**************************************************************************************************

Pro read_lna_gain, calibration_dir_lna, LnaGainFilename, headerlines, FreqAxis, Glna


if (NOT File_Test(calibration_dir_lna + LnaGainFilename)) then begin
   print, 'Error: Cannot open LNA gain file ' + calibration_dir_lna + LnaGainFilename
   retall
endif

numsamples = File_Lines(calibration_dir_lna + LnaGainFilename) - headerlines
frequency = fltarr(numsamples)
Glna = fltarr(numsamples)
freq_scale = 1.0

;---Read header
openr, lun, calibration_dir_lna + LnaGainFilename, /get_lun
if (headerlines gt 0) then begin
   header = strarr(headerlines)
   readf, lun, header

   ;---Detect frequency units
   if (strpos(header[headerlines - 1], 'MHz') ge 0) then begin
      freq_scale = 1e6
   endif
endif

;---Read data
row = ''
for i=0L, (numsamples - 1L) do begin
   readf, lun, row
   temp = strsplit(row, ',', /extract)
   if (strtrim(temp[0], 2) ne '') then begin
      frequency[i] = float(temp[0])
      Glna[i] = float(temp[1])
   endif
endfor
free_lun, lun

;---Discard empty rows
valid_data = where(frequency ne 0.0)
frequency = frequency(valid_data)
Glna = Glna(valid_data)

;---Warn if frequency is out of range
frequency = frequency * freq_scale
if ((min(FreqAxis) lt min(frequency)) OR (max(FreqAxis) gt max(frequency))) then begin
   print, 'Warning: Requested frequency for LNA gain is out of range!'
endif

;---Calculate interpolated LNA gain
Glna = interpol(Glna, frequency, FreqAxis, /spline)

;---Clip extrapolated values to last valid value
index = where(FreqAxis lt min(frequency), count)
if (count gt 0) then begin
   Glna[index] = Glna[(index[count - 1] + 1) < (n_elements(Glna) - 1)]
endif
index = where(FreqAxis gt max(frequency), count)
if (count gt 0) then begin
   Glna[index] = Glna[(index[0] - 1) > 0]
endif

End


;**************************************************************************************************
; Read Cable Loss
;**************************************************************************************************

Pro read_cable_loss, calibration_dir_cable, CableLossFilename, headerlines, FreqAxis, CableLoss


if (NOT File_Test(calibration_dir_cable + CableLossFilename)) then begin
   print, 'Error: Cannot open Cable Loss file ' + calibration_dir_cable + CableLossFilename
   retall
endif

numsamples = File_Lines(calibration_dir_cable + CableLossFilename) - headerlines
frequency = fltarr(numsamples)
CableLoss = fltarr(numsamples)
freq_scale = 1.0

;---Read header
openr, lun, calibration_dir_cable + CableLossFilename, /get_lun
if (headerlines gt 0) then begin
   header = strarr(headerlines)
   readf, lun, header

   ;---Detect frequency units
   if (strpos(header[headerlines - 1], 'MHz') ge 0) then begin
      freq_scale = 1e6
   endif
endif

;---Read data
row = ''
for i=0L, (numsamples - 1L) do begin
   readf, lun, row
   temp = strsplit(row, ',', /extract)
   if (strtrim(temp[0], 2) ne '') then begin
      frequency[i] = float(temp[0])
      CableLoss[i] = float(temp[1])
   endif
endfor
free_lun, lun

;---Discard empty rows
valid_data = where(frequency ne 0.0)
frequency = frequency(valid_data)
CableLoss = CableLoss(valid_data)

;---Warn if frequency is out of range
frequency = frequency * freq_scale
if ((min(FreqAxis) lt min(frequency)) OR (max(FreqAxis) gt max(frequency))) then begin
   print, 'Warning: Requested frequency for Cable Loss is out of range!'
endif

;---Calculate interpolated Cable Loss
CableLoss = interpol(CableLoss, frequency, FreqAxis, /spline)

;---Clip extrapolated values to last valid value
index = where(FreqAxis lt min(frequency), count)
if (count gt 0) then begin
   CableLoss[index] = CableLoss[(index[count - 1] + 1) < (n_elements(CableLoss) - 1)]
endif
index = where(FreqAxis gt max(frequency), count)
if (count gt 0) then begin
   CableLoss[index] = CableLoss[(index[0] - 1) > 0]
endif


End


;**************************************************************************************************
; Read Antenna Gain
;**************************************************************************************************

Pro read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation, AntennaGain


if (NOT File_Test(calibration_dir_ant + AntGainFilename)) then begin
   print, 'Error: Cannot open file containing antenna gain pattern ' + calibration_dir_ant + AntGainFilename
   retall
endif

numsamples = File_Lines(calibration_dir_ant + AntGainFilename) - headerlines
frequency = dblarr(numsamples)
AntennaGain = fltarr(2, numsamples)
freq_scale = 1d

;---Read header
openr, lun, calibration_dir_ant + AntGainFilename, /get_lun
if (headerlines gt 0) then begin
   header = strarr(headerlines)
   readf, lun, header

   ;---Detect frequency units
   if (strpos(header[headerlines - 1], 'MHz') ge 0) then begin
      freq_scale = 1d6
   endif
endif

;---Read data
row = ''
for i=0L, (numsamples - 1L) do begin
   readf, lun, row
   temp = strsplit(row, ',', /extract)
   if ((n_elements(temp) lt 3) AND (polarisation eq 1)) then begin
      print, 'Error: Antenna H-Pol calibration data not available!'
      retall
   endif
   if (strtrim(temp[0], 2) ne '') then begin
      frequency[i] = double(temp[0])
      AntennaGain[0, i] = float(temp[1])
      if (polarisation eq 1) then AntennaGain[1, i] = float(temp[2])
   endif
endfor
free_lun, lun

;---Discard empty rows
valid_data = where(frequency ne 0.0)
frequency = frequency(valid_data)
AntennaGain = AntennaGain(*, valid_data)

;---Warn if frequency is out of range
frequency = frequency * freq_scale
if ((min(FreqAxis) lt min(frequency)) OR (max(FreqAxis) gt max(frequency))) then begin
   print, 'Warning: Requested frequency for Antenna Gain Pattern is out of range!'
endif

;---Return interpolated antenna gain
if (polarisation eq 0) then begin
   index = where(AntennaGain[0, *] ne 0.0, count)
   if (count eq 0) then begin
      print, 'Warning: Antenna V-Pol calibration data not available! Using zeroes!'
   endif
   AntennaGain = interpol(reform(AntennaGain[0, *]), frequency, FreqAxis, /spline)
endif else begin
   index = where(AntennaGain[1, *] ne 0.0, count)
   if (count eq 0) then begin
      print, 'Warning: Antenna H-Pol calibration data not available! Using zeroes!'
   endif
   AntennaGain = interpol(reform(AntennaGain[1, *]), frequency, FreqAxis, /spline)
endelse

;---Clip extrapolated values to last valid value
index = where(FreqAxis lt min(frequency), count)
if (count gt 0) then begin
   AntennaGain[index] = AntennaGain[(index[count - 1] + 1) < (n_elements(AntennaGain) - 1)]
endif
index = where(FreqAxis gt max(frequency), count)
if (count gt 0) then begin
   AntennaGain[index] = AntennaGain[(index[0] - 1) > 0]
endif

End


;**************************************************************************************************
; Get CCF
;**************************************************************************************************

Pro read_ccf, FreqAxis, chamber_type, CCF_filename, headerlines, CCF


if (chamber_type eq 1) then begin

   if (NOT File_Test(CCF_filename)) then begin
      print, 'Error: Cannot open CCF file ' + CCF_filename
      retall
   endif

   numsamples = File_Lines(CCF_filename) - headerlines
   frequency = fltarr(numsamples)
   CCF = fltarr(numsamples)

   ;---Read header
   openr, lun, CCF_filename, /get_lun
   if (headerlines gt 0) then begin
      header = strarr(headerlines)
      readf, lun, header
   endif

   ;---Read data
   row = ''
   for i=0L, (numsamples - 1L) do begin
      readf, lun, row
      temp = strsplit(row, ',', /extract)
      if (strtrim(temp[0], 2) ne '') then begin
         frequency[i] = float(temp[0])
         CCF[i] = float(temp[1])
      endif
   endfor
   free_lun, lun

   ;---Discard empty rows
   valid_data = where(frequency ne 0.0)
   frequency = frequency(valid_data)
   CCF = CCF(valid_data)

   ;---Warn if frequency is out of range
   if ((min(FreqAxis) lt min(frequency)) OR (max(FreqAxis) gt max(frequency))) then begin
      print, 'Warning: Requested frequency for CCF is out of range!'
   endif

   ;---Calculate interpolated Cable Loss
   CCF = interpol(CCF, frequency, FreqAxis, /spline)

   ;---Clip extrapolated values to last valid value
   index = where(FreqAxis lt min(frequency), count)
   if (count gt 0) then begin
      CCF[index] = CCF[(index[count - 1] + 1) < (n_elements(CCF) - 1)]
   endif
   index = where(FreqAxis gt max(frequency), count)
   if (count gt 0) then begin
      CCF[index] = CCF[(index[0] - 1) > 0]
   endif

   ;---Convert do dB
   CCF = 10.0 * alog10(CCF)

endif

End


;**************************************************************************************************
; Get RFI Thresholds
;**************************************************************************************************

Pro get_RFI_thresholds, FreqAxis, range, culprit_distance, enclosure_attenuation, ContinuumE_interp, $
                        SpectralE_interp


;---Constants
Z = 119.9169832 * !pi

;---Threshold requirements for a device located 10m from the nearest antenna focus
Freq       = [  70.000, 499.9999,   500.000, 2000.000, 25500.000] * 1d6                     ;Hz
ContinuumP = [ -89.937, -81.22459, -125.268, -117.604,   -84.512] + enclosure_attenuation   ;dBm
SpectralP  = [-104.937, -96.22459, -140.268, -132.604,   -99.512] + enclosure_attenuation   ;dBm

;---Apply adjustment above 500 MHz
adjustment = 10.0 * alog10((culprit_distance / 10.0) ^ 2.0)
index = where(Freq ge 500d6)
ContinuumP[index] = ContinuumP[index] + adjustment
SpectralP[index]  = SpectralP[index]  + adjustment

;---Convert to electric field strength
ContinuumE = ContinuumP + (10.0 * alog10(Z / (4.0 * !pi * (range ^ 2.0)))) + 90.0    ;dBuV/m
SpectralE  = SpectralP  + (10.0 * alog10(Z / (4.0 * !pi * (range ^ 2.0)))) + 90.0    ;dBuV/m

;---Interpolate
N = n_elements(Freq)
;ContinuumE_interp = fltarr(n_elements(FreqAxis))
SpectralE_interp = fltarr(n_elements(FreqAxis))
for i=0L, (N - 2L) do begin
   index = where((FreqAxis ge Freq[i]) AND (FreqAxis lt Freq[i+1]), count)
   if (count gt 0L) then begin
      ;ContinuumE_interp[index] = ContinuumE[i] + (ContinuumE[i+1] - ContinuumE[i]) * alog10(FreqAxis[index] / Freq[i]) / alog10(Freq[i+1] / Freq[i])
      SpectralE_interp[index] = SpectralE[i] + (SpectralE[i+1] - SpectralE[i]) * alog10(FreqAxis[index] / Freq[i]) / alog10(Freq[i+1] / Freq[i])
   endif
endfor
index = where(FreqAxis lt Freq[0], count)
if (count gt 0L) then begin
   ;ContinuumE_interp[index] = ContinuumE[0] + (ContinuumE[1] - ContinuumE[0]) * alog10(FreqAxis[index] / Freq[0]) / alog10(Freq[1] / Freq[0])
   SpectralE_interp[index] = SpectralE[0] + (SpectralE[1] - SpectralE[0]) * alog10(FreqAxis[index] / Freq[0]) / alog10(Freq[1] / Freq[0])
endif
index = where(FreqAxis ge Freq[N - 1L], count)
if (count gt 0L) then begin
   ;ContinuumE_interp[index] = ContinuumE[N-2] + (ContinuumE[N-1] - ContinuumE[N-2]) * alog10(FreqAxis[index] / Freq[N-2]) / alog10(Freq[N-1] / Freq[N-2])
   SpectralE_interp[index] = SpectralE[N-2] + (SpectralE[N-1] - SpectralE[N-2]) * alog10(FreqAxis[index] / Freq[N-2]) / alog10(Freq[N-1] / Freq[N-2])
endif

End


;**************************************************************************************************
; Test if replot_old_data is incorrectly set
;**************************************************************************************************

Pro test_replot_old_data, prev_filename, DataFilename, RBW, prev_RBW

if (n_elements(prev_filename) eq 0) then begin
   print, 'Error: First set replot_old_data=0.'
   retall
endif
if (DataFilename ne prev_filename) then begin
   print, 'Error: Filename has changed, but replot_old_data=1. First set replot_old_data=0.
   retall
endif
if (RBW ne prev_RBW) then begin
   print, 'Error: RBW has changed, but replot_old_data=1. First set replot_old_data=0.
   retall
endif

End


;**************************************************************************************************
; Plot Spectrum
;**************************************************************************************************

Pro plot_spectrum, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, SamplingFrequency, $
                   fft_length, num_ffts, Frequency, AcqBW, subtitle, output_dir, filename_base, DataFilename, $
                   date, xaxis_centre, xaxis_span, yaxis_bottom_add, yaxis_top_add, show_plot, AntennaGain, $
                   generate_rfi_flaglist, cutoff_above_noisefloor, show_threshold


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Find peak
fbinsize = SamplingFrequency / fft_length
tint = (fft_length * num_ffts) / SamplingFrequency
peak_index = where(Spectrum eq max(Spectrum))
peak_freq = (peak_index[0] * fbinsize) + FreqAxis[0]               ;[Hz]
peak = Spectrum[peak_index[0]]                                     ;[dBuV/m]


;---Find Tsys
index_fc = fft_length / 2
noise_floor_dBm = noise_floor_E[index_fc] - (20.0 * alog10(FreqAxis[index_fc])) + AntennaGain[index_fc] - $
                  (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
Tsys = (10.0 ^ (noise_floor_dBm / 10.0)) / (k * fbinsize * 1000.0)


;---Prepare plot annotations
ytitle = 'Electric Field Strength [dBuV/m]'
noise_unit = ' dBµV/m'
if (subtitle ne '') then begin
   title = 'Spectrum of ' + DataFilename + ' on ' + date + '!C!C' + subtitle
   ymargin = [8.5, 6]
endif else begin
   title = 'Spectrum of ' + DataFilename + ' on ' + date
   ymargin = [8.5, 4]
endelse
if (FreqAxis[0] ge 1e9) then begin
   freq_scale = 1e9
   freq_unit = '[GHz]'
endif else begin
   freq_scale = 1e6
   freq_unit = '[MHz]'
endelse
if (AcqBW ge 1e6) then begin
   acqbw_scale = 1e6
   acqbw_unit = ' [MHz]'
endif else begin
   acqbw_scale = 1e3
   acqbw_unit = ' [kHz]'
endelse
if (SamplingFrequency ge 1e6) then begin
   fadc_scale = 1e6
   fadc_unit = ' [MHz]'
endif else begin
   fadc_scale = 1e3
   fadc_unit = ' [kHz]'
endelse
if (tint ge 1.0) then begin
   tint_scale = 1.0
   tint_unit = ' s'
endif else begin
   tint_scale = 1e3
   tint_unit = ' ms'
endelse
if (fbinsize ge 1e3) then begin
   rbw_scale = 1e3
   rbw_unit = ' kHz'
endif else begin
   rbw_scale = 1.0
   rbw_unit = ' Hz'
endelse
if (show_threshold eq 1) then begin
   minval = min([Spectrum, SpectralE_interp])
   maxval = max([Spectrum, SpectralE_interp])
endif else begin
   minval = min(Spectrum)
   maxval = max(Spectrum)
endelse
xleft  = ((Frequency - (SamplingFrequency / 2.0)) > (xaxis_centre - (xaxis_span / 2.0))) / freq_scale
xright = ((Frequency + (SamplingFrequency / 2.0)) < (xaxis_centre + (xaxis_span / 2.0))) / freq_scale
subtitle1 = '!CPeak of ' + strtrim(string(peak, format='(F15.3)'), 2) + ' dBuV/m at ' + $
            strtrim(string(peak_freq / 1e6, format='(F15.3)'), 2) + ' MHz, Tsys = ' + $
            strtrim(round(Tsys), 2) + ' K'
subtitle2 = '!C!Cf!Dc!N = ' + strtrim(string(Frequency / 1e6, format='(F15.3)'), 2) + ' MHz' + $
            ', AcqBW = ' + strtrim(string(AcqBW / acqbw_scale, format='(F15.3)'), 2) + acqbw_unit + $
            ', f!Dadc!N = ' + strtrim(string(SamplingFrequency / fadc_scale, format='(F15.3)'), 2) + fadc_unit + $
            ', t!Dint!N = ' + strtrim(string(tint / tint_scale, format='(F15.3)'), 2) + tint_unit + $
            ', RBW = ' + strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2) + rbw_unit


;---Plot unfiltered spectrum
window, /free, xsize=1000, ysize=600, title=DataFilename
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
charsize  = 1.08
plot, FreqAxis / freq_scale, Spectrum, $
      title  = title, $
      subtitle = subtitle1 + subtitle2, $
      xtitle = 'Frequency ' + freq_unit, $
      ytitle = ytitle, $
      xstyle = 1, ystyle = 1, $
      xmargin = [10, 3], ymargin = ymargin, $
      charsize = charsize, $
      xrange = [xleft, xright], $
      yrange = [minval - 0.0 + yaxis_bottom_add, maxval + 1.0 + yaxis_top_add]


;---Over-Plot noise floor
oplot, FreqAxis / freq_scale, noise_floor_E, color='18E218'xL


;---Over-Plot threshold levels
if (show_threshold eq 1) then begin
   oplot, FreqAxis / freq_scale, SpectralE_interp, linestyle = 2, color='0000FF'xL
   ;oplot, FreqAxis / freq_scale, ContinuumE_interp, linestyle = 4, color='0000FF'xL
endif


;---Over-Plot RFI Flag List threshold level
if (generate_rfi_flaglist eq 1) then begin
   oplot, FreqAxis / freq_scale, noise_floor_E + cutoff_above_noisefloor, color='FF0000'xL
endif


;---Save plot
file_mkdir, output_dir + strmid(DataFilename, 0, strlen(DataFilename) - 4)
filename_base = output_dir + strmid(DataFilename, 0, strlen(DataFilename) - 4) + '/' + $
                strmid(DataFilename, 0, strlen(DataFilename) - 4) + ', RBW = ' + $
                strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2) + rbw_unit
write_png, filename_base + ', Spectrum.png', tvrd(true=1)
if (show_plot eq 0) then wdelete


End


;**************************************************************************************************
; Plot Transients
;**************************************************************************************************

Pro plot_transients, S, FreqAxis, transient_RFI_detection, SamplingFrequency, num_ffts, fft_length, $
                     Frequency, filename_base, show_plot

if (transient_RFI_detection eq 1) then begin
   if (num_ffts ge 10) then begin
      sigma_peak_index = where(S eq max(S))
      sigma_peak_freq = (sigma_peak_index[0] * (SamplingFrequency / fft_length)) + Frequency - (SamplingFrequency / 2.0)
      if (FreqAxis[0] ge 1e9) then begin
         freq_scale = 1e9
         freq_unit = '[GHz]'
      endif else begin
         freq_scale = 1e6
         freq_unit = '[MHz]'
      endelse
      window, /free, xsize=1000, ysize=600
      plot, FreqAxis / freq_scale, S, $
            title = 'Normalised Standard Deviation', $
            subtitle = '!CPeak at ' + strtrim(string(sigma_peak_freq / 1e6, format='(F15.3)'), 2) + ' MHz', $
            xtitle = 'Frequency ' + freq_unit, $
            ytitle = 'Normalised StdDev [dB]', $
            xstyle = 1, ystyle = 1, $
            xmargin = [9, 5], ymargin = [6.5, 4], $
            charsize = 1.18, $
            yrange = [min(S) - 2, max(S) + 2], $
            xrange = [min(FreqAxis), max(FreqAxis)] / freq_scale
      write_png, filename_base + ', Normalised StdDev.png', tvrd(true=1)
      if (show_plot eq 0) then wdelete
   endif else begin
      print, 'Warning: Transient RFI detection was not performed, because number of FFTs was less than 10!'
   endelse
endif

End


;**************************************************************************************************
; Save Calibration Plots
;**************************************************************************************************

Pro save_calibration_plots, FreqAxis, AntennaGain, AntGainFilename, polarisation, Glna, $
                            LnaGainFilename, CableLoss, CableLossFilename, output_dir, DataFilename, $
                            chamber_type, CCF_filename, CCF


if (polarisation eq 0) then begin
   polarisation_str = 'V-pol'
endif else begin
   polarisation_str = 'H-pol'
endelse
if (FreqAxis[0] ge 1e9) then begin
   freq_scale = 1e9
   freq_unit = '[GHz]'
endif else begin
   freq_scale = 1e6
   freq_unit = '[MHz]'
endelse
filename_base = output_dir + strmid(DataFilename, 0, strlen(DataFilename) - 4) + '\' + $
                strmid(DataFilename, 0, strlen(DataFilename) - 4)


;---Save antenna gain pattern plot
window, /free, xsize=1000, ysize=600
plot, FreqAxis / freq_scale, AntennaGain, $
      title = 'Antenna Gain Pattern (' + AntGainFilename + ') ' + polarisation_str, $
      xtitle = 'Frequency ' + freq_unit, $
      ytitle = 'Antenna Gain [dBi]', $
      xstyle = 1, ystyle = 1, $
      xmargin = [9, 5], ymargin = [5.5, 6], $
      charsize = 1.18, $
      yrange = [min(AntennaGain) - 2, max(AntennaGain) + 2], $
      xrange = [min(FreqAxis), max(FreqAxis)] / freq_scale
write_png, filename_base + ', Grx, ' + polarisation_str + '.png', tvrd(true=1)
wdelete


;---Save LNA gain pattern plot
window, /free, xsize=1000, ysize=600
plot, FreqAxis / freq_scale, Glna, $
      title = 'LNA Gain Pattern (' + LnaGainFilename + ')', $
      xtitle = 'Frequency ' + freq_unit, $
      ytitle = 'LNA Gain [dB]', $
      xstyle = 1, ystyle = 1, $
      xmargin = [9, 5], ymargin = [5.5, 6], $
      charsize = 1.18, $
      yrange = [min(Glna) - 2, max(Glna) + 2], $
      xrange = [min(FreqAxis), max(FreqAxis)] / freq_scale
write_png, filename_base + ', LNA Gain.png', tvrd(true=1)
wdelete


;---Save cable loss pattern plot
window, /free, xsize=1000, ysize=600
plot, FreqAxis / freq_scale, CableLoss, $
      title = 'Cable Loss Pattern (' + CableLossFilename + ')', $
      xtitle = 'Frequency ' + freq_unit, $
      ytitle = 'Cable Loss [dB]', $
      xstyle = 1, ystyle = 1, $
      xmargin = [9, 5], ymargin = [5.5, 6], $
      charsize = 1.18, $
      yrange = [min(CableLoss) - 2, max(CableLoss) + 2], $
      xrange = [min(FreqAxis), max(FreqAxis)] / freq_scale
write_png, filename_base + ', Cable Loss.png', tvrd(true=1)
wdelete


;---Save CCF plot
if (chamber_type eq 1) then begin
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / freq_scale, CCF, $
         title = 'CCF (' + CCF_filename + ')', $
         xtitle = 'Frequency ' + freq_unit, $
         ytitle = 'CCF [dB]', $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min(CCF) - 2, max(CCF) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / freq_scale
   write_png, filename_base + ', CCF.png', tvrd(true=1)
   wdelete
endif


End


;**************************************************************************************************
; Save Summary
;**************************************************************************************************

Pro save_summary, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, filename_base, date, DataFilename, $
                  AntGainFilename, LnaGainFilename, CableLossFilename, polarisation, Preamp, RFAttenuation, $
                  AntennaGain, Glna, CableLoss, Frequency, AcqBW, SamplingFrequency, num_ffts, fft_length, $
                  range, start_time, fraction_to_process, num_peaks, plant_signals, freq_signal_at_threshold, $
                  freq_signal_at_Xsigma, threshold_margin, rfi_frequency, rfi_plus_noise_power, rfi_power, $
                  rfi_margin, rfi_3dB_bandwidth, rfi_SNR, subtitle, remove_dc_offset, culprit_distance, $
                  enclosure_attenuation, apply_rfi_flaglist, rfi_flaglist_folder, rfi_flaglist_filename, $
                  chamber_type, antenna_efficiency, CCF_filename


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Statistics at Centre Frequency
index_fc = fft_length / 2
fbinsize = SamplingFrequency / fft_length
noise_floor_dBm = noise_floor_E[index_fc] - (20.0 * alog10(FreqAxis[index_fc])) + AntennaGain[index_fc] - $
                  (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
noise_floor_dBmHz = noise_floor_dBm - (10.0 * alog10(fbinsize))
Tsys = (10.0 ^ (noise_floor_dBm / 10.0)) / (k * fbinsize * 1000.0)
sigma_dBm = 10.0 * alog10(k * Tsys * SamplingFrequency * 1000.0 / (fft_length * sqrt(num_ffts)))
threshold_level_dBuVm = 10.0 * alog10((10.0 ^ (SpectralE_interp[index_fc] / 10.0)) -  (10.0 ^ (noise_floor_E[index_fc] / 10.0))) ;[dBuV/m]
threshold_level_dBm = threshold_level_dBuVm - (20.0 * alog10(FreqAxis[index_fc])) + AntennaGain[index_fc] - $
                      (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
SNR = threshold_level_dBm - sigma_dBm
tint = (fft_length * num_ffts) / SamplingFrequency

if (polarisation eq 0) then polarisation_str = 'V-pol' else polarisation_str = 'H-pol'
if (Preamp eq 'true') then Preamp_str = 'ON' else Preamp_str = 'OFF'
if (remove_dc_offset eq 1) then dc_offset_str = 'yes' else dc_offset_str = 'no'
if (chamber_type eq 0) then begin
   chamber_str = 'Anechoic Chamber'
endif else begin
   chamber_str = 'Reverberation Chamber'
endelse
if (FreqAxis[index_fc] ge 1e9) then begin
   freq_scale = 1e9
   freq_unit = '[GHz]'
endif else begin
   freq_scale = 1e6
   freq_unit = '[MHz]'
endelse
if (AcqBW ge 1e6) then begin
   acqbw_scale = 1e6
   acqbw_unit = '[MHz]'
endif else begin
   acqbw_scale = 1e3
   acqbw_unit = '[kHz]'
endelse
if (SamplingFrequency ge 1e6) then begin
   fadc_scale = 1e6
   fadc_unit = '[MHz]'
endif else begin
   fadc_scale = 1e3
   fadc_unit = '[kHz]'
endelse
if (tint ge 1.0) then begin
   tint_scale = 1.0
   tint_unit = '[s]'
endif else begin
   tint_scale = 1e3
   tint_unit = '[ms]'
endelse
if (fbinsize ge 1e3) then begin
   rbw_scale = 1e3
   rbw_unit = '[kHz]'
endif else begin
   rbw_scale = 1.0
   rbw_unit = '[Hz]'
endelse


openw, lun, filename_base + ', Spreadsheet.csv', /get_lun, error=err
if (err ne 0) then begin
   openw, lun, filename_base + ', Spreadsheet (1).csv', /get_lun
endif
printf, lun, ''
printf, lun, 'Data Filename'         + ',' + DataFilename
printf, lun, 'Acquisition Date'      + ',' + date
printf, lun, 'Description'           + ',' + subtitle
printf, lun, 'Chamber Type'          + ',' + chamber_str
printf, lun, 'Antenna Gain Filename' + ',' + AntGainFilename
printf, lun, 'LNA Gain Filename'     + ',' + LnaGainFilename
printf, lun, 'Cable Loss Filename'   + ',' + CableLossFilename
if (chamber_type eq 1) then begin
   printf, lun, 'CCF Filename'       + ',' + CCF_filename
   printf, lun, 'Antenna efficiency' + ',' + strtrim(string(antenna_efficiency, format='(F15.3)'), 2)
   printf, lun, 'Polarisation'       + ',' + 'Not Applicable'
endif else begin
   printf, lun, 'Polarisation'       + ',' + polarisation_str
endelse
printf, lun, 'Pre-amplifier'         + ',' + Preamp_str
printf, lun, 'RF attenuation'        + ',' + strtrim(string(RFAttenuation, format='(F15.1)'), 2)                   + ',' + '[dB]'
printf, lun, 'f_c'                   + ',' + strtrim(string(FreqAxis[index_fc] / freq_scale, format='(F18.6)'), 2) + ',' + freq_unit
printf, lun, 'B_t'                   + ',' + strtrim(string(AcqBW / acqbw_scale, format='(F18.6)'), 2)             + ',' + acqbw_unit
printf, lun, 'f_adc'                 + ',' + strtrim(string(SamplingFrequency / fadc_scale, format='(F18.6)'), 2)  + ',' + fadc_unit
printf, lun, 'Number of FFTs (N)'    + ',' + strtrim(num_ffts, 2)
printf, lun, 'FFT length (L)'        + ',' + strtrim(fft_length, 2)
printf, lun, 'B_ch'                  + ',' + strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2)            + ',' + rbw_unit
printf, lun, 't_i'                   + ',' + strtrim(string(tint / tint_scale, format='(F15.3)'), 2)               + ',' + tint_unit
printf, lun, 'range'                 + ',' + strtrim(string(range, format='(F15.1)'), 2)                           + ',' + '[m]'
printf, lun, 'Culprit distance'      + ',' + strtrim(string(culprit_distance, format='(F15.1)'), 2)                + ',' + '[m]'
printf, lun, 'Enclosure attenuation' + ',' + strtrim(string(enclosure_attenuation, format='(F15.1)'), 2)           + ',' + '[dB]'
printf, lun, 'Remove DC offset'      + ',' + strtrim(dc_offset_str, 2)
if (apply_rfi_flaglist eq 1) then begin
   if (NOT File_Test(rfi_flaglist_folder + rfi_flaglist_filename)) then begin
      printf, lun, 'Apply RFI Flaglist' + ',' + 'RFI Flaglist File could not be opened'
   endif else begin
      printf, lun, 'Apply RFI Flaglist' + ',' + 'yes'
   endelse
   printf, lun, 'RFI Flaglist File'  + ',' + rfi_flaglist_folder + rfi_flaglist_filename
endif else begin
   printf, lun, 'Apply RFI Flaglist' + ',' + 'no'
endelse
printf, lun, 'Processing time'       + ',' + strtrim(string(systime(1) - start_time, format='(F15.1)'), 2)         + ',' + '[s]'
if (fraction_to_process lt 1.0) then begin
   printf, lun, 'Warning: Only a fraction (' + strtrim(fraction_to_process, 2) + ') of the entire data was processed!'
endif
printf, lun, ''
printf, lun, 'At f_c (' + strtrim(FreqAxis[index_fc] / freq_scale, 2) + ' ' + freq_unit + '):'
printf, lun, 'G_rx'                  + ',' + strtrim(string(AntennaGain[index_fc], format='(F15.1)'), 2)           + ',' + '[dB]'
printf, lun, 'G_lna'                 + ',' + strtrim(string(Glna[index_fc], format='(F15.1)'), 2)                  + ',' + '[dB]'
printf, lun, 'L_cable'               + ',' + strtrim(string(CableLoss[index_fc], format='(F15.1)'), 2)             + ',' + '[dB]'
printf, lun, 'NF_E'                  + ',' + strtrim(string(noise_floor_E[index_fc], format='(F15.1)'), 2)         + ',' + '[dBµV/m]'
printf, lun, 'NF_dBm'                + ',' + strtrim(string(noise_floor_dBm, format='(F15.1)'), 2)                 + ',' + '[dBm]'    + ',' + 'referenced at LNA input'
printf, lun, 'NF_dBm/Hz'             + ',' + strtrim(string(noise_floor_dBmHz, format='(F15.1)'), 2)               + ',' + '[dBm/Hz]' + ',' + 'referenced at LNA input'
printf, lun, 'Tsys'                  + ',' + strtrim(string(Tsys, format='(F15.1)'), 2)                            + ',' + '[K]'      + ',' + 'referenced at LNA input'
printf, lun, 'sigma_dBm'             + ',' + strtrim(string(sigma_dBm, format='(F15.1)'), 2)                       + ',' + '[dBm]'    + ',' + 'referenced at LNA input'
printf, lun, 'SNR'                   + ',' + strtrim(string(SNR, format='(F15.1)'), 2)                             + ',' + '[dB]'     + ',' + 'Threshold Level to Sigma Ratio'
printf, lun, ''

rfi_comment = strarr(n_elements(rfi_frequency))
if (plant_signals eq 1) then begin
   index = where(rfi_frequency eq freq_signal_at_threshold, count)
   if (count gt 0) then rfi_comment[index] = 'Planted signal'
   index = where(rfi_frequency eq freq_signal_at_Xsigma, count)
   if (count gt 0) then rfi_comment[index] = 'Planted signal'
endif
valid_index = where(rfi_frequency ne 0.0, count)
if (count gt 0) then begin
   rfi_frequency        = rfi_frequency[valid_index]
   rfi_plus_noise_power = rfi_plus_noise_power[valid_index]
   rfi_power            = rfi_power[valid_index]
   rfi_3dB_bandwidth    = rfi_3dB_bandwidth[valid_index]
   rfi_SNR              = rfi_SNR[valid_index]
   rfi_comment          = rfi_comment[valid_index]
   sort_index           = sort(rfi_frequency)
   printf, lun, 'Top ' + strtrim(num_peaks, 2) + ' RFI spikes within valid AcqBW region:' + ',' + $
                strtrim(string((FreqAxis[index_fc] - (AcqBW / 2.0)) / freq_scale, format='(F18.6)'), 2) + ' ' + freq_unit + ' < f < ' + $
                strtrim(string((FreqAxis[index_fc] + (AcqBW / 2.0)) / freq_scale, format='(F18.6)'), 2) + ' ' + freq_unit
   printf, lun, ''
   printf, lun, '#' + ',' + 'Frequency' + ',' + 'Above Threshold' + ',' + 'RFI + Noise' + ',' + 'RFI' + ',' + '3dB BW' + ',' + 'SNR' + ',' + 'Comment'
   printf, lun, '' + ',' + '[MHz]' + ',' + '[dB]' + ',' + '[dBµV/m]' + ',' + '[dBµV/m]' + ',' + '[Hz]' + ',' + '[dB]'
   for i=0, (n_elements(valid_index) - 1) do begin
      index = sort_index[i]
      printf, lun, strtrim(i+1, 2) + ',' + $
                   strtrim(string(rfi_frequency[index] / 1e6,  format='(F15.6)'), 2) + ',' + $
                   strtrim(string(rfi_margin[index],           format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_plus_noise_power[index], format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_power[index],            format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_3dB_bandwidth[index],    format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_SNR[index],              format='(F15.1)'), 2) + ',' + $
                   rfi_comment[index]
   endfor
endif else begin
   printf, lun, 'No RFI spikes within ' + strtrim(string(threshold_margin, format='(F15.1)'), 2) + ' dB of threshold'
endelse

free_lun, lun

End


;**************************************************************************************************
; Save Spectrum
;**************************************************************************************************

Pro save_spectrum, Spectrum, FreqAxis, write_spectrum_to_csv, filename_base, Frequency, $
                   SamplingFrequency, AcqBW

if (write_spectrum_to_csv eq 1) then begin
   openw, lun, filename_base + ', Spectrum.csv', /get_lun
   printf, lun, 'NumberSamples, ' + strtrim(n_elements(FreqAxis), 2)
   printf, lun, 'CentreFrequency, ' + string(Frequency, format='(F12.1)') + ', [Hz]'
   printf, lun, 'SamplingFrequency, ' + string(SamplingFrequency, format='(F12.1)') + ', [Hz]'
   printf, lun, 'AcquisitionBandwidth, ' + string(AcqBW, format='(F12.1)') + ', [Hz]'
   printf, lun, transpose(strtrim(string(FreqAxis, format='(F15.1)'), 2) + ', ' + $
                strtrim(string(Spectrum, format='(F15.1)'), 2))
   free_lun, lun
endif

End


;**************************************************************************************************
; Get Stitch Spectra Filenames
;**************************************************************************************************

Pro get_stitch_spectra_filenames, output_dir, Folder_list, RBW, stitch_filenames

stitch_filenames = strarr(n_elements(Folder_list))

for i=0, (n_elements(Folder_list) - 1) do begin
   filenames = File_Search(output_dir + Folder_list[i] + '\*Spectrum.csv', count=count)
   if (count gt 0) then begin
      for j=0, (n_elements(filenames) - 1) do begin
         filename_rbw = strsplit(filenames[j], '=', /extract)
         if (strpos(filename_rbw[1], 'kHz') ge 0) then begin
            filename_rbw = strsplit(filename_rbw[1], 'kHz', /extract)
            filename_rbw = float(filename_rbw[0]) * 1e3
         endif else begin
            filename_rbw = strsplit(filename_rbw[1], 'Hz', /extract)
            filename_rbw = float(filename_rbw[0])
         endelse
         if ((filename_rbw gt (RBW - (RBW / 10.0))) AND (filename_rbw lt (RBW + (RBW / 10.0)))) then begin
            stitch_filenames[i] = filenames[j]
         endif
      endfor
      if (stitch_filenames[i] eq '') then begin
         print, 'Warning: ' + output_dir + Folder_list[i] + $
                '\ does not contain a Spectrum.csv file with RBW = ' + strtrim(RBW, 2) + ' Hz'
      endif
   endif else begin
      print, 'Warning: ' + output_dir + Folder_list[i] + $
             '\ does not exist or does not contain any valid Spectrum.csv files'
   endelse
endfor

valid_index = where(stitch_filenames ne '', count)
if (count gt 0) then begin
   stitch_filenames = stitch_filenames[valid_index]
   Folder_list = Folder_list[valid_index]
endif else begin
   print, 'Error: No valid Spectrum.csv files were found to stitch together for RBW = ' + $
          strtrim(RBW, 2) + ' Hz'
   stitch_filenames = 0
endelse

End


;**************************************************************************************************
; Get Stitch Spectra Parameters
;**************************************************************************************************

Pro get_stitch_spectra_parameters, Folder_list, stitch_filenames, NumberSamples, CentreFrequency, $
                                   SamplingFrequency, AcqBW

num_spectra       = n_elements(stitch_filenames)
NumberSamples     = lonarr(num_spectra)
CentreFrequency   = dblarr(num_spectra)
SamplingFrequency = fltarr(num_spectra)
AcqBW             = fltarr(num_spectra)

for i=0, (num_spectra - 1) do begin
   openr, lun, stitch_filenames[i], /get_lun
   temp = ''
   readf, lun, temp
   temp = strsplit(temp, ',', /extract)
   NumberSamples[i] = long(temp[1])
   temp = ''
   readf, lun, temp
   temp = strsplit(temp, ',', /extract)
   CentreFrequency[i] = double(temp[1])
   temp = ''
   readf, lun, temp
   temp = strsplit(temp, ',', /extract)
   SamplingFrequency[i] = float(temp[1])
   temp = ''
   readf, lun, temp
   temp = strsplit(temp, ',', /extract)
   AcqBW[i] = float(temp[1])
   free_lun, lun
endfor

;---Verify that all spectra have same fbinsize
fbinsize = SamplingFrequency / NumberSamples
index = where (fbinsize ne fbinsize[0], count)
if (count gt 0) then begin
   print, 'Error: Different spectra have different frequency bin sizes!'
   retall
endif

;---Sort in ascending order
index = sort(CentreFrequency)
Folder_list       = Folder_list[index]
stitch_filenames  = stitch_filenames[index]
NumberSamples     = NumberSamples[index]
CentreFrequency   = CentreFrequency[index]
SamplingFrequency = SamplingFrequency[index]
AcqBW             = AcqBW[index]

End


;**************************************************************************************************
; Get Stitch Spectra Spectrum
;**************************************************************************************************

Pro get_stitch_spectra_Spectrum, Spectrum, FreqAxis, stitch_filenames, SamplingFrequency, $
                                 NumberSamples, CentreFrequency, AcqBW, shrinkfactor, maxbins

num_spectra = n_elements(CentreFrequency)
fbinsize = SamplingFrequency[0] / NumberSamples[0]
start_freq = CentreFrequency[0] - (SamplingFrequency[0] / 2.0)
stop_freq  = CentreFrequency[num_spectra - 1] + (SamplingFrequency[num_spectra - 1] / 2.0)
numbins = long((stop_freq - start_freq) / fbinsize)
shrinkfactor = 1

if (numbins gt maxbins) then begin
   shrinkfactor = ceil(float(numbins) / maxbins)
   print, 'Warning: Composite array is too large. Summing ' + strtrim(shrinkfactor, 2) + ' bins to reduce size.'
   numbins = numbins / shrinkfactor
endif


if (num_spectra gt 1) then begin

   FreqAxis    = (Dindgen(numbins) * fbinsize * shrinkfactor) + start_freq
   Spectrum    = fltarr(numbins)
   Spectrum[*] = !values.f_nan

   for i=0, (num_spectra - 1) do begin
      print, '        Stitching ' + stitch_filenames[i]
      header = strarr(4)
      data = dblarr(2, NumberSamples[i])
      openr, lun, stitch_filenames[i], /get_lun
      readf, lun, header
      readf, lun, data
      free_lun, lun

      if (i eq 0) then begin                         ;Keep the first bit of the first spectrum
         start_index = where(FreqAxis ge data[0, 0])
         start_index = start_index[0]
         if (start_index gt 0) then begin            ;Cut off first few empty bins
            FreqAxis = FreqAxis[start_index : *]
            Spectrum = Spectrum[start_index : *]
            start_index = 0
         endif
         num_data_bins = (floor(NumberSamples[i] / float(shrinkfactor)) < n_elements(FreqAxis)) * shrinkfactor
         end_data_index = num_data_bins - 1L
         Spectrum[0 : (num_data_bins / shrinkfactor) - 1] = 10.0 * alog10(Total(Reform(10.0 ^ $
                  (data[1, 0 : end_data_index] / 10.0), shrinkfactor, (num_data_bins / shrinkfactor)), 1))
      endif else begin
         start_freq = ((((CentreFrequency[i-1] + (AcqBW[i-1] / 2.0)) + (CentreFrequency[i] - (AcqBW[i] / 2.0))) / 2.0) > $
                      (CentreFrequency[i-1] + (AcqBW[i-1] / 2.0))) > data[0, 0]
         start_data_index = where(data[0, *] ge start_freq)
         start_data_index = start_data_index[0]
         start_index = where(FreqAxis ge start_freq)
         start_index = start_index[0]
         bins_remaining = floor((NumberSamples[i] - start_data_index) / float(shrinkfactor)) * shrinkfactor
         end_index = (start_index + (bins_remaining / shrinkfactor) - 1) < (n_elements(FreqAxis) - 1L)
         bins_remaining = (end_index - start_index + 1) * shrinkfactor
         end_data_index = start_data_index + bins_remaining - 1
         Spectrum[start_index : end_index] = 10.0 * alog10(Total(Reform(10.0 ^ $
                  (data[1, start_data_index : end_data_index] / 10.0), shrinkfactor, (bins_remaining / shrinkfactor)), 1))
      endelse
   endfor

   if (end_index lt (n_elements(FreqAxis) - 1L)) then begin   ;Cut off last few empty bins
      FreqAxis = FreqAxis[0 : end_index]
      Spectrum = Spectrum[0 : end_index]
   endif

endif else begin
   header = strarr(4)
   data = dblarr(2, NumberSamples[0])
   openr, lun, stitch_filenames[0], /get_lun
   readf, lun, header
   readf, lun, data
   free_lun, lun
   FreqAxis = data[0, *]
   Spectrum = data[1, *]
endelse


End


;**************************************************************************************************
; Get Stitch Spectra Values
;**************************************************************************************************

Pro get_stitch_spectra_values, FreqAxis, CentreFrequency, SamplingFrequency, NumberSamples, $
                               fft_length, AcqBW, Frequency, shrinkfactor

num_spectra = n_elements(CentreFrequency)
fft_length = n_elements(FreqAxis)
fbinsize = (SamplingFrequency[0] / NumberSamples[0]) * shrinkfactor
SamplingFrequency = fft_length * fbinsize
AcqBW = (CentreFrequency[num_spectra - 1] + (AcqBW[num_spectra - 1] / 2.0)) - $
        (CentreFrequency[0] - (AcqBW[0] / 2.0))
Frequency = FreqAxis[fft_length / 2]

End


;**************************************************************************************************
; Get Stitch Spectra Plot
;**************************************************************************************************

Pro get_stitch_spectra_plot, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, $
                             Folder_list, SamplingFrequency, fft_length, Frequency, AcqBW, subtitle, $
                             output_dir, DataFilename, AntennaGain, filename_base, show_plot, $
                             show_threshold


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Find peak
fbinsize = SamplingFrequency / fft_length
good = where(Finite(Spectrum))
peak_index = good[where(Spectrum[good] eq max(Spectrum[good]))]
peak_freq = (peak_index[0] * fbinsize) + FreqAxis[0]               ;[Hz]
peak = Spectrum[peak_index[0]]                                     ;[dBuV/m]


;---Find Tsys
index_fc = fft_length / 2
noise_floor_dBm = noise_floor_E[index_fc] - (20.0 * alog10(FreqAxis[index_fc])) + $
                  AntennaGain[index_fc] - (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
Tsys = (10.0 ^ (noise_floor_dBm / 10.0)) / (k * fbinsize * 1000.0)


;---Prepare plot annotations
FolderName = Folder_list[0] + ' - ' + Folder_list[n_elements(Folder_list) - 1]
ytitle = 'Electric Field Strength [dBuV/m]'
noise_unit = ' dBµV/m'
if (subtitle ne '') then begin
   title = 'Combined Spectrum of ' + Folder_list[0] + ' through ' + $
            Folder_list[n_elements(Folder_list) - 1] + '!C!C' + subtitle
   ymargin = [8.5, 6]
endif else begin
   title = 'Combined Spectrum of ' + Folder_list[0] + ' through ' + $
            Folder_list[n_elements(Folder_list) - 1]
   ymargin = [8.5, 4]
endelse
if (FreqAxis[0] ge 1e9) then begin
   freq_scale = 1e9
   freq_unit = '[GHz]'
endif else begin
   freq_scale = 1e6
   freq_unit = '[MHz]'
endelse
if (AcqBW ge 1e6) then begin
   acqbw_scale = 1e6
   acqbw_unit = ' [MHz]'
endif else begin
   acqbw_scale = 1e3
   acqbw_unit = ' [kHz]'
endelse
if (SamplingFrequency ge 1e6) then begin
   fadc_scale = 1e6
   fadc_unit = ' [MHz]'
endif else begin
   fadc_scale = 1e3
   fadc_unit = ' [kHz]'
endelse
if (fbinsize ge 1e3) then begin
   rbw_scale = 1e3
   rbw_unit = ' kHz'
endif else begin
   rbw_scale = 1.0
   rbw_unit = ' Hz'
endelse
if (show_threshold eq 1) then begin
   minval = min([Spectrum[good], SpectralE_interp])
   maxval = max([Spectrum[good], SpectralE_interp])
endif else begin
   minval = min(Spectrum[good])
   maxval = max(Spectrum[good])
endelse
subtitle1 = '!CPeak of ' + strtrim(string(peak, format='(F15.3)'), 2) + ' dBuV/m at ' + $
            strtrim(string(peak_freq / 1e6, format='(F15.3)'), 2) + ' MHz, Tsys = ' + $
            strtrim(round(Tsys), 2) + ' K'
subtitle2 = '!C!Cf!Dc!N = ' + strtrim(string(Frequency / 1e6, format='(F15.3)'), 2) + ' MHz' + $
            ', AcqBW = ' + strtrim(string(AcqBW / acqbw_scale, format='(F15.3)'), 2) + acqbw_unit + $
            ', f!Dadc!N = ' + strtrim(string(SamplingFrequency / fadc_scale, format='(F15.3)'), 2) + fadc_unit + $
            ', RBW = ' + strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2) + rbw_unit


;---Plot unfiltered spectrum
window, /free, xsize=1000, ysize=600, title=FolderName
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
charsize  = 1.08
plot, FreqAxis / freq_scale, Spectrum, $
      title  = title, $
      subtitle = subtitle1 + subtitle2, $
      xtitle = 'Frequency ' + freq_unit, $
      ytitle = ytitle, $
      xstyle = 1, ystyle = 1, $
      xmargin = [10, 3], ymargin = ymargin, $
      charsize = charsize, $
      xrange = [(Frequency - (SamplingFrequency / 2.0)) / freq_scale, $
                (Frequency + (SamplingFrequency / 2.0)) / freq_scale], $
      yrange = [minval - 0.0, maxval + 1.0]


;---Over-Plot noise floor
oplot, FreqAxis / freq_scale, noise_floor_E, color='18E218'xL


;---Over-Plot threshold levels
if (show_threshold eq 1) then begin
   oplot, FreqAxis / freq_scale, SpectralE_interp, linestyle = 2, color='0000FF'xL
   ;oplot, FreqAxis / freq_scale, ContinuumE_interp, linestyle = 4, color='0000FF'xL
endif


;---Save plot
file_mkdir, output_dir + FolderName
filename_base = output_dir + FolderName + '\' + FolderName + ', RBW = ' + $
                strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2) + rbw_unit
write_png, filename_base + ', Spectrum.png', tvrd(true=1)
if (show_plot eq 0) then wdelete

End


;**************************************************************************************************
; Save Stitch Spectra Summary
;**************************************************************************************************

Pro stitch_spectra_summary, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, filename_base, $
                            Folder_list, AntGainFilename, AntennaGain, polarisation, Frequency, $
                            AcqBW, SamplingFrequency, fft_length, subtitle, range, culprit_distance, $
                            enclosure_attenuation, start_time, num_peaks, threshold_margin, $
                            rfi_frequency, rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, $
                            chamber_type


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Statistics at Centre Frequency
index_fc = fft_length / 2
fbinsize = SamplingFrequency / fft_length
noise_floor_dBm = noise_floor_E[index_fc] - (20.0 * alog10(FreqAxis[index_fc])) + AntennaGain[index_fc] - $
                  (10.0 * alog10((4.0 * !pi * Z) / (light ^ 2.0))) - 90.0
noise_floor_dBmHz = noise_floor_dBm - (10.0 * alog10(fbinsize))
Tsys = (10.0 ^ (noise_floor_dBm / 10.0)) / (k * fbinsize * 1000.0)

if (chamber_type eq 0) then begin
   chamber_str = 'Anechoic Chamber'
endif else begin
   chamber_str = 'Reverberation Chamber'
endelse
if (polarisation eq 0) then polarisation_str = 'V-pol' else polarisation_str = 'H-pol'
if (FreqAxis[index_fc] ge 1e9) then begin
   freq_scale = 1e9
   freq_unit = '[GHz]'
endif else begin
   freq_scale = 1e6
   freq_unit = '[MHz]'
endelse
if (AcqBW ge 1e6) then begin
   acqbw_scale = 1e6
   acqbw_unit = '[MHz]'
endif else begin
   acqbw_scale = 1e3
   acqbw_unit = '[kHz]'
endelse
if (SamplingFrequency ge 1e6) then begin
   fadc_scale = 1e6
   fadc_unit = '[MHz]'
endif else begin
   fadc_scale = 1e3
   fadc_unit = '[kHz]'
endelse
if (fbinsize ge 1e3) then begin
   rbw_scale = 1e3
   rbw_unit = '[kHz]'
endif else begin
   rbw_scale = 1.0
   rbw_unit = '[Hz]'
endelse

openw, lun, filename_base + ', Spreadsheet.csv', /get_lun, error=err
if (err ne 0) then begin
   openw, lun, filename_base + ', Spreadsheet (1).csv', /get_lun
endif
printf, lun, ''
printf, lun, 'Filename Range'        + ',' + Folder_list[0] + ' to ' + Folder_list[n_elements(Folder_list) - 1]
printf, lun, 'Description'           + ',' + subtitle
printf, lun, 'Chamber Type'          + ',' + chamber_str
printf, lun, 'Antenna Gain Filename' + ',' + AntGainFilename
if (chamber_type eq 1) then begin
   printf, lun, 'Polarisation'       + ',' + 'Not Applicable'
endif else begin
   printf, lun, 'Polarisation'       + ',' + polarisation_str
endelse
printf, lun, 'f_c'                   + ',' + strtrim(string(FreqAxis[index_fc] / freq_scale, format='(F18.6)'), 2) + ',' + freq_unit
printf, lun, 'B_t'                   + ',' + strtrim(string(AcqBW / acqbw_scale, format='(F18.6)'), 2)             + ',' + acqbw_unit
printf, lun, 'f_adc'                 + ',' + strtrim(string(SamplingFrequency / fadc_scale, format='(F18.6)'), 2)  + ',' + fadc_unit
printf, lun, 'FFT length (L)'        + ',' + strtrim(fft_length, 2)
printf, lun, 'B_ch'                  + ',' + strtrim(string(fbinsize / rbw_scale, format='(F15.3)'), 2)            + ',' + rbw_unit
printf, lun, 'range'                 + ',' + strtrim(string(range, format='(F15.1)'), 2)                           + ',' + '[m]'
printf, lun, 'Culprit distance'      + ',' + strtrim(string(culprit_distance, format='(F15.1)'), 2)                + ',' + '[m]'
printf, lun, 'Enclosure attenuation' + ',' + strtrim(string(enclosure_attenuation, format='(F15.1)'), 2)           + ',' + '[dB]'
printf, lun, 'Processing time'       + ',' + strtrim(systime(1) - start_time, 2)                                   + ',' + '[s]'
printf, lun, ''
printf, lun, 'At f_c (' + strtrim(FreqAxis[index_fc] / freq_scale, 2) + ' ' + freq_unit + '):'
printf, lun, 'G_rx'                  + ',' + strtrim(AntennaGain[index_fc], 2)                                     + ',' + '[dB]'
printf, lun, 'NF_E'                  + ',' + strtrim(noise_floor_E[index_fc], 2)                                   + ',' + '[dBµV/m]'
printf, lun, 'NF_dBm'                + ',' + strtrim(noise_floor_dBm, 2)                                           + ',' + '[dBm]'    + ',' + 'referenced at LNA input'
printf, lun, 'NF_dBm/Hz'             + ',' + strtrim(noise_floor_dBmHz, 2)                                         + ',' + '[dBm/Hz]' + ',' + 'referenced at LNA input'
printf, lun, 'Tsys'                  + ',' + strtrim(Tsys, 2)                                                      + ',' + '[K]'      + ',' + 'referenced at LNA input'
printf, lun, ''

valid_index = where(rfi_frequency ne 0.0, count)
if (count gt 0) then begin
   rfi_frequency        = rfi_frequency[valid_index]
   rfi_plus_noise_power = rfi_plus_noise_power[valid_index]
   rfi_power            = rfi_power[valid_index]
   rfi_3dB_bandwidth    = rfi_3dB_bandwidth[valid_index]
   sort_index           = sort(rfi_frequency)
   printf, lun, 'Top ' + strtrim(num_peaks, 2) + ' RFI spikes within valid AcqBW region:' + ',' + $
                strtrim(string((FreqAxis[index_fc] - (AcqBW / 2.0)) / freq_scale, format='(F18.6)'), 2) + ' ' + freq_unit + ' < f < ' + $
                strtrim(string((FreqAxis[index_fc] + (AcqBW / 2.0)) / freq_scale, format='(F18.6)'), 2) + ' ' + freq_unit
   printf, lun, ''
   printf, lun, '#' + ',' + 'Frequency' + ',' + 'Above Threshold' + ',' + 'RFI + Noise' + ',' + 'RFI' + ',' + '3dB BW'
   printf, lun, '' + ',' + '[MHz]' + ',' + '[dB]' + ',' + '[dBµV/m]' + ',' + '[dBµV/m]' + ',' + '[Hz]'
   for i=0, (n_elements(valid_index) - 1) do begin
      index = sort_index[i]
      printf, lun, strtrim(i+1, 2) + ',' + $
                   strtrim(string(rfi_frequency[index] / 1e6,  format='(F15.6)'), 2) + ',' + $
                   strtrim(string(rfi_margin[index],           format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_plus_noise_power[index], format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_power[index],            format='(F15.1)'), 2) + ',' + $
                   strtrim(string(rfi_3dB_bandwidth[index],    format='(F15.1)'), 2)
   endfor
endif else begin
   printf, lun, 'No RFI spikes within ' + strtrim(threshold_margin, 2) + ' dB of threshold'
endelse

free_lun, lun

End


;**************************************************************************************************
; Cleanup
;**************************************************************************************************

Pro cleanup, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, AntennaGain, $
             Glna, CableLoss, S, CCF

Spectrum = 0
FreqAxis = 0
SpectralE_interp = 0
ContinuumE_interp = 0
noise_floor_E = 0
AntennaGain = 0
Glna = 0
CableLoss = 0
S = 0
CCF = 0

End


;**************************************************************************************************
; Process Data
;**************************************************************************************************

Pro process_data, DataFilename, RBW, polarisation, AntGainFilenameNumber, AntGainFilenameList, $
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


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Get Filenames
get_filenames, AntGainFilenameNumber, AntGainFilenameList, AntGainFilename, $
               LnaGainFilenameNumber, LnaGainFilenameList, LnaGainFilename, $
               CableLossFilenameNumber, CableLossFilenameList, CableLossFilename


;---Error checks
error_checks, fraction_to_process, input_dir, calibration_dir_ant, calibration_dir_lna, $
              calibration_dir_cable, DataFilename, AntGainFilename, LnaGainFilename, CableLossFilename, $
              culprit_distance, enclosure_attenuation, chamber_type, CCF_filename


;---Read header information
read_header, lun, offset, input_dir, DataFilename, SamplingFrequency, NumberSamples, Scaling, $
             date, Frequency, AcqBW, RFAttenuation, Preamp


;---Calculate FFT length and number of FFTs to be averaged
get_fft_length, FFT_power_of_2, SamplingFrequency, RBW, NumberSamples, fraction_to_process, $
                fft_length, num_ffts


;---Get frequency axis
FreqAxis = (Dindgen(fft_length) * (SamplingFrequency / fft_length)) + (Frequency - (SamplingFrequency / 2.0))


;---Get calibration data as a function of frequency
read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation, AntennaGain
read_cable_loss, calibration_dir_cable, CableLossFilename, headerlines, FreqAxis, CableLoss
read_lna_gain, calibration_dir_lna, LnaGainFilename, headerlines, FreqAxis, Glna
read_ccf, FreqAxis, chamber_type, CCF_filename, 0, CCF


;---Calculate averaged magnitude FFT [dBm]
get_power_spectrum_dBm, Spectrum, S, input_dir, DataFilename, offset, fft_length, num_ffts, Scaling, $
                        transient_RFI_detection, peak_hold


;---Convert to electric field strength [dBuV/m]
get_efield_strength, Spectrum, FreqAxis, AntennaGain, Glna, CableLoss, CCF, chamber_type, $
                     antenna_efficiency, range


;---Find median-filtered noise floor level [dBuV/m]
get_noise_floor, Spectrum, FreqAxis, SamplingFrequency, fft_length, noise_floor_E, AcqBW


;---Generate RFI Flag List
get_rfi_flaglist, Spectrum, FreqAxis, noise_floor_E, SamplingFrequency, generate_rfi_flaglist, $
                  cutoff_above_noisefloor, rfi_flaglist_filename, rfi_flaglist_folder


;---Get RFI threshold levels
get_RFI_thresholds, FreqAxis, range, culprit_distance, enclosure_attenuation, ContinuumE_interp, SpectralE_interp


;---Add median-filtered noise floor level to RFI threshold levels
;ContinuumE_interp = 10.0 * alog10((10.0 ^ (Temporary(ContinuumE_interp) / 10.0)) + (10.0 ^ (noise_floor_E / 10.0)))  ;[dBuV/m]
SpectralE_interp = 10.0 * alog10((10.0 ^ (Temporary(SpectralE_interp) / 10.0)) + (10.0 ^ (noise_floor_E / 10.0)))    ;[dBuV/m]


;---Remove DC offset spike
remove_DC, Spectrum, noise_floor_E, remove_dc_offset, fft_length


;---Apply RFI Flag List
remove_rfi_flags, Spectrum, FreqAxis, noise_floor_E, apply_rfi_flaglist, rfi_flaglist_filename, rfi_flaglist_folder


;---Plant signals
plant_test_signals, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, plant_signals, freq_signal_at_threshold, $
                    freq_signal_at_Xsigma, Xsigma, Frequency, AcqBW


;---Find RFI spikes
find_rfi_spikes, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, AntennaGain, $
                 num_peaks, threshold_margin, SamplingFrequency, fft_length, num_ffts, $
                 rfi_frequency, rfi_plus_noise_power, rfi_power, rfi_margin, rfi_3dB_bandwidth, $
                 rfi_SNR, Frequency, AcqBW

End


;**************************************************************************************************
; Save Result
;**************************************************************************************************

Pro save_result, Spectrum, S, FreqAxis, noise_floor_E, SpectralE_interp, ContinuumE_interp, AntennaGain, $
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


;---Constants
Z     = 119.9169832 * !pi
k     = 1.3806503e-23
light = 299792458.0


;---Plot spectrum
plot_spectrum, Spectrum, FreqAxis, SpectralE_interp, ContinuumE_interp, noise_floor_E, SamplingFrequency, $
               fft_length, num_ffts, Frequency, AcqBW, subtitle, output_dir, filename_base, DataFilename, $
               date, xaxis_centre, xaxis_span, yaxis_bottom_add, yaxis_top_add, show_plot, AntennaGain, $
               generate_rfi_flaglist, cutoff_above_noisefloor, show_threshold


;---Plot transients
plot_transients, S, FreqAxis, transient_RFI_detection, SamplingFrequency, num_ffts, fft_length, $
                 Frequency, filename_base, show_plot


;---Save calibration plots
save_calibration_plots, FreqAxis, AntennaGain, AntGainFilename, polarisation, Glna, LnaGainFilename, $
                        CableLoss, CableLossFilename, output_dir, DataFilename, chamber_type, $
                        CCF_filename, CCF


;---Save useful information to file
save_summary, Spectrum, FreqAxis, noise_floor_E, SpectralE_interp, filename_base, date, DataFilename, $
              AntGainFilename, LnaGainFilename, CableLossFilename, polarisation, Preamp, RFAttenuation, $
              AntennaGain, Glna, CableLoss, Frequency, AcqBW, SamplingFrequency, num_ffts, fft_length, $
              range, start_time, fraction_to_process, num_peaks, plant_signals, freq_signal_at_threshold, $
              freq_signal_at_Xsigma, threshold_margin, rfi_frequency, rfi_plus_noise_power, rfi_power, $
              rfi_margin, rfi_3dB_bandwidth, rfi_SNR, subtitle, remove_dc_offset, culprit_distance, $
              enclosure_attenuation, apply_rfi_flaglist, rfi_flaglist_folder, rfi_flaglist_filename, $
              chamber_type, antenna_efficiency, CCF_filename


;---Save spectrum to csv file
save_spectrum, Spectrum, FreqAxis, write_spectrum_to_csv, filename_base, Frequency, SamplingFrequency, AcqBW


End


;**************************************************************************************************
; Dummy procedure to ensure that all above routines are compiled
;**************************************************************************************************

Pro Subroutines


End


;**************************************************************************************************
; This routine (subroutines.pro) may also be run as a main procedure
;**************************************************************************************************


test_antenna_gain        = 0
test_lna_gain            = 1
test_cable_loss          = 0
test_rfi_thresholds      = 0
test_ccf                 = 0


AntGainFilenameNumber    = 1                           ;Refer to numbers below
LnaGainFilenameNumber    = 4                           ;Refer to numbers below
CableLossFilenameNumber  = 1                           ;Refer to numbers below

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
                            'No_Cable']                                     ;#16

output_dir               = 'plots\'
calibration_dir_ant      = 'Equipment_Database\Passive_Antennas\'
calibration_dir_lna      = 'Equipment_Database\LNAs\'
calibration_dir_cable    = 'Equipment_Database\Cables\'
CCF_filename             = 'Equipment_Database\Reverb_Chamber\ACF.csv'


if (test_antenna_gain eq 1) then begin
   AntGainFilename  = AntGainFilenameList[AntGainFilenameNumber - 1] + '.csv'
   fmin             = 300.0
   fmax             = 400.0
   polarisation1    = 0                          ;0 = V polarisation, 1 = H polarisation
   polarisation2    = 1                          ;0 = V polarisation, 1 = H polarisation
   headerlines      = 2
   ;numpoints = round(fmax - fmin + 1)
   numpoints = 100e6 / 10000
   FreqAxis = (((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin) * 1e6
   read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation1, AntennaGain_V
   read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis, polarisation2, AntennaGain_H
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / 1e6, AntennaGain_V, $
         title = 'Antenna Gain Pattern (' + AntGainFilename + ')', $
         xtitle = 'Frequency [MHz]', $
         ytitle = 'Antenna Gain [dBi]', $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min([AntennaGain_V, AntennaGain_H]) - 2, max([AntennaGain_V, AntennaGain_H]) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / 1e6
   oplot, FreqAxis / 1e6, AntennaGain_H, linestyle=2
   write_png, output_dir + 'Antenna_Gain.png', tvrd(true=1)
endif

if (test_lna_gain eq 1) then begin
   LnaGainFilename   = LnaGainFilenameList[LnaGainFilenameNumber - 1] + '.csv'
   fmin              = 50.0
   fmax              = 3500.0
   headerlines       = 2
   numpoints = round(fmax - fmin + 1)
   FreqAxis = (((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin) * 1e6
   read_lna_gain, calibration_dir_lna, LnaGainFilename, headerlines, FreqAxis, Glna
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / 1e6, Glna, $
         title = 'LNA Gain (' + LnaGainFilename + ')', $
         xtitle = 'Frequency [MHz]', $
         ytitle = 'LNA Gain [dB]', $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min(Glna) - 2, max(Glna) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / 1e6
   write_png, output_dir + 'LNA_Gain.png', tvrd(true=1)
endif

if (test_cable_loss eq 1) then begin
   CableLossFilename = CableLossFilenameList[CableLossFilenameNumber - 1] + '.csv'
   fmin              = 50.0
   fmax              = 6000.0
   headerlines       = 2
   numpoints = round(fmax - fmin + 1)
   FreqAxis = (((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin) * 1e6
   read_cable_loss, calibration_dir_cable, CableLossFilename, headerlines, FreqAxis, CableLoss
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / 1e6, CableLoss, $
         title = 'Cable Loss (' + CableLossFilename + ')', $
         xtitle = 'Frequency [MHz]', $
         ytitle = 'Cable Loss [dB]', $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min(CableLoss) - 2, max(CableLoss) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / 1e6
   write_png, output_dir + 'Cable_Loss.png', tvrd(true=1)
endif

if (test_rfi_thresholds eq 1) then begin
   fmin = 70.0
   fmax = 3000.0
   range = 1.0
   culprit_distance = 10.0
   enclosure_attenuation = 0.0
   numpoints = round(fmax - fmin + 1)
   FreqAxis = ((((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin)) * 1d6
   get_RFI_thresholds, FreqAxis, range, culprit_distance, enclosure_attenuation, ContinuumE_interp, $
                       SpectralE_interp
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / 1d6, SpectralE_interp, $
         title = 'Spectral Line RFI Threshold', $
         xtitle = 'Frequency [MHz]', $
         ytitle = 'RFI Threshold Level [dBuV/m]', $
         /xlog, $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min(SpectralE_interp) - 2, max(SpectralE_interp) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / 1d6
   write_png, output_dir + 'RFI_spectral_threshold.png', tvrd(true=1)
endif

if (test_ccf eq 1) then begin
   chamber_type      = 1
   fmin              = 70.0
   fmax              = 8000.0
   headerlines       = 0
   numpoints = round(fmax - fmin + 1)
   FreqAxis = (((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin) * 1e6
   read_ccf, FreqAxis, chamber_type, CCF_filename, headerlines, CCF
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL
   window, /free, xsize=1000, ysize=600
   plot, FreqAxis / 1e6, CCF, $
         title = 'CCF (' + CCF_filename + ')', $
         xtitle = 'Frequency [MHz]', $
         ytitle = 'CCF [dB]', $
         xstyle = 1, ystyle = 1, $
         xmargin = [9, 5], ymargin = [5.5, 6], $
         charsize = 1.18, $
         yrange = [min(CCF) - 2, max(CCF) + 2], $
         xrange = [min(FreqAxis), max(FreqAxis)] / 1e6
   write_png, output_dir + 'CCF.png', tvrd(true=1)
endif


End