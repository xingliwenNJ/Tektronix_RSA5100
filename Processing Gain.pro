
;**************************************************************************************************
; Author:       R. Lord
; Organisation: SKA SA
; Ver:          1.0
;**************************************************************************************************


tint                     = [7.05, 76.3, 343.52]       ;[s]
SamplingFrequency        = [150e6, 12.5e6, 1.5625e6]  ;[Hz]
Bt                       = [85e6, 10e6, 1.25e6]       ;[Hz]
RBW                      = 10.0                       ;Desired approximate RBW (FFT length will be rounded to nearest power of 2)
Tsys                     = 870.0                      ;[K]
range                    = 1.0                        ;Distance from DUT to receiving antenna
culprit_distance         = 10.0                       ;Distance between Antenna Focus and DUT [m] (E.g. 10m or 30m for AP, 1m for Receiver and Digitiser)
enclosure_attenuation    = 0.0                        ;Additional attenuation provided by enclosure [dB]
FFT_power_of_2           = 0                          ;0=no, 1=yes (if yes, FFT length will be rounded to nearest power of 2, but then RBW will be approximate)
polarisation             = 0                          ;0=V polarisation, 1=H polarisation



;---Constants
k = 1.3806503e-23
Z = 119.9169832 * !pi
light = 299792458.0
fmin = 70.0
fmax = 3000.0


;---Compile routines
Resolve_Routine, 'Subroutines'


;---Get Antenna Positioner threshold levels in dBuV/m
numpoints = round(fmax - fmin + 1)
FreqAxis = ((Findgen(numpoints) / (numpoints - 1.0)) * (fmax - fmin)) + fmin
get_RFI_thresholds, FreqAxis * 1d6, range, culprit_distance, enclosure_attenuation, ContinuumE_interp, SpectralE_interp


;---Get calibration data as a function of frequency
headerlines              = 2
calibration_dir_ant      = 'Equipment_Database\Passive_Antennas\'
AntGainFilenameNumber    = 1                             ;Refer to numbers below
AntGainFilenameList      = ['Antenna_MESA_GLPDA'                        , $ ;#1
                            'Antenna_MESA_KLPDA1'                       , $ ;#2
                            'Antenna_MESA_KLPDA3'                       , $ ;#3
                            'Antenna_Houwteq_EMCO3115_FRH_Small'        , $ ;#4
                            'Antenna_Houwteq_EMCO3106B_FRH_Large'       , $ ;#5
                            'Antenna_Houwteq_EMCO3141_BiConiLog_FLPDA'  , $ ;#6
                            'No_Antenna']                                   ;#7 (Use this option for direct injection into RTSA)

;---Get Antenna Gain filename
if ((AntGainFilenameNumber lt 1) OR (AntGainFilenameNumber gt n_elements(AntGainFilenameList))) then begin
   print, 'Error: Number supplied for Antenna Gain filename must be >=1 and <=' + strtrim(n_elements(AntGainFilenameList), 2)
   retall
endif
AntGainFilename = AntGainFilenameList[AntGainFilenameNumber - 1] + '.csv'
read_antenna_gain, calibration_dir_ant, AntGainFilename, headerlines, FreqAxis * 1e6, polarisation, AntennaGain


;---Calculate Antenna Positioner threshold levels in dBm, as measured at the Spectrum Analyzer
Pr = SpectralE_interp - 90.0 + AntennaGain - (20.0 * alog10(FreqAxis * 1e6)) + (10.0 * alog10((light ^ 2.0) / (4.0 * !pi * Z)))


for i=0, (n_elements(tint) - 1) do begin

   ;---Calculate FFT length and number of FFTs to be averaged
   NumberSamples = floor(tint[i] * SamplingFrequency[i])

   ;---Calculate FFT length and number of FFTs to be averaged
   if (FFT_power_of_2 eq 1) then begin
      fft_length = 2L ^ round(alog(SamplingFrequency[i] / RBW) / alog(2))  ;power of 2
   endif else begin
      fft_length = round((SamplingFrequency[i] / RBW) / 2.0) * 2L          ;even number
   endelse
   num_ffts = floor(NumberSamples / fft_length)
   if (num_ffts lt 1) then begin
      fft_length = NumberSamples
      num_ffts = 1
   endif
   fbinsize = SamplingFrequency[i] / fft_length
   tint_actual = (fft_length * num_ffts) / SamplingFrequency[i]

   print, 'RBW = ' + strtrim(fbinsize, 2)
   print, 'FFT length = ' + strtrim(fft_length, 2)
   print, 'Number of FFTs averaged = ' + strtrim(num_ffts, 2)


   ;---Calculate sigma (stddev) of noise floor
   sigma = 10.0 * alog10((k * Tsys * SamplingFrequency[i] * 1000.0) / (fft_length * sqrt(num_ffts)))


   ;---Calculate RFI-to-Sigma ratio in dB
   SNR = Pr - sigma


   ;---Calculate RFI-to-Noisefloor ratio in dB
   NF = 10.0 * alog10(k * Tsys * (SamplingFrequency[i] / fft_length) * 1000.0)
   RFI_to_noisefloor = 10.0 * alog10((10.0 ^ (Pr / 10.0)) + (10.0 ^ (NF / 10.0))) - NF


   ;---Prepare plot annotations
   Tsys_str = strtrim(round(Tsys * 10.0) / 10.0, 2)
   Tsys_str = strmid(Tsys_str, 0, strpos(Tsys_str, '.') + 2)
   Bt_str = strtrim(round(Bt[i] / 1e6 * 1000.0) / 1000.0, 2)
   Bt_str = strmid(Bt_str, 0, strpos(Bt_str, '.') + 4)
   if (tint_actual ge 1.0) then begin
      tint_str = strtrim(round(tint_actual * 1000.0) / 1000.0, 2)
      tint_str = strmid(tint_str, 0, strpos(tint_str, '.') + 4) + ' s'
   endif else begin
      tint_str = strtrim(round(tint_actual * 1000.0 * 10.0) / 10.0, 2)
      tint_str = strmid(tint_str, 0, strpos(tint_str, '.') + 2) + ' ms'
   endelse
   if (fbinsize ge 1000.0) then begin
      fbinsize_str = strtrim(round(fbinsize / 1000.0 * 1000.0) / 1000.0, 2)
      fbinsize_str = strmid(fbinsize_str, 0, strpos(fbinsize_str, '.') + 4) + ' kHz'
   endif else begin
      fbinsize_str = strtrim(round(fbinsize * 1000.0) / 1000.0, 2)
      fbinsize_str = strmid(fbinsize_str, 0, strpos(fbinsize_str, '.') + 4) + ' Hz'
   endelse


   ;---Plot results
   output_dir = 'plots/'
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL

   xstart = 600
   xstop = 790
   y1 = 23.2
   y2 = 25.2
   y3 = 27.2
   if (i eq 0) then begin
      window, 0, xsize=1000, ysize=600
      plot, FreqAxis, SNR, $
            /xlog, $
            title = 'RSA RFI-to-Sigma Ratio', $
            xtitle = 'Frequency [MHz]', $
            ytitle = 'RFI-to-Sigma [dB]', $
            xstyle = 1, ystyle = 1, $
            xmargin = [9, 5], ymargin = [5.5, 4], $
            charsize = 1.18, $
            xrange = [500, fmax], $
            yrange = [0, 30]

      plots, [xstart, xstop], [y1, y1]
      xyouts, 820, y1-0.2, 'B!Dt!N = ' + Bt_str + ' MHz', charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif

   if (i eq 1) then begin
      oplot, FreqAxis, SNR, $
             linestyle=4

      plots, [xstart, xstop], [y2, y2], linestyle=4
      xyouts, 820, y2-0.2, 'B!Dt!N = ' + Bt_str + ' MHz', charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif

   if (i eq 2) then begin
      oplot, FreqAxis, SNR, $
             linestyle=2

      oplot, [fmin, fmax], [8.45, 8.45], $
             linestyle=1

      plots, [xstart, xstop], [y3, y3], linestyle=2
      xyouts, 820, y3-0.2, 'B!Dt!N = ' + Bt_str + ' MHz', charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif


   index = where(FreqAxis ge 520.0)
   print, 'Expected RFI-to-Sigma Ratio at 520 MHz = ' + strtrim(SNR[index[0]], 2) + ' [dB]'
;   print, 'Expected RFI-to-Noisefloor Ratio at 520 MHz = ' + strtrim(RFI_to_noisefloor[index[0]], 2) + ' [dB]'
   print, ''


endfor

write_png, output_dir + 'RSA RFI-to-Sigma Ratio.png', tvrd(true=1)







;---Now do the same for Ratty

tint              = [100, 1000, 5000.0]       ;[s]
SamplingFrequency = [856e6, 856e6, 856e6]    ;[Hz]
Bt                = [770e6, 770e6, 770e6]       ;[Hz]
RBW               = 10.0                 ;Desired approximate RBW (FFT length will be rounded to nearest power of 2)
Tsys              = 500.0                ;[K]
range             = 1.0                  ;Distance from DUT to receiving antenna
relaxation        = 0.0                  ;Relaxation of threshold level


for i=0, (n_elements(tint) - 1) do begin


   ;---Calculate FFT length and number of FFTs to be averaged
   NumberSamples = floor(tint[i] * SamplingFrequency[i], /L64)
   fft_length = 1048576L
   num_ffts = floor(NumberSamples / fft_length)
   if (num_ffts lt 1) then begin
      fft_length = NumberSamples
      num_ffts = 1
   endif
   fbinsize = SamplingFrequency[i] / fft_length
   tint_actual = (fft_length * num_ffts) / SamplingFrequency[i]

   print, 'RBW = ' + strtrim(fbinsize, 2)
   print, 'FFT length = ' + strtrim(fft_length, 2)
   print, 'Number of FFTs averaged = ' + strtrim(num_ffts, 2)


   ;---Calculate sigma (stddev) of noise floor
   sigma = 10.0 * alog10((k * Tsys * SamplingFrequency[i] * 1000.0) / (fft_length * sqrt(num_ffts)))


   ;---Calculate RFI-to-Sigma ratio in dB
   SNR = Pr - sigma


   ;---Calculate RFI-to-Noisefloor ratio in dB
   NF = 10.0 * alog10(k * Tsys * (SamplingFrequency[i] / fft_length) * 1000.0)
   RFI_to_noisefloor = 10.0 * alog10((10.0 ^ (Pr / 10.0)) + (10.0 ^ (NF / 10.0))) - NF


   ;---Prepare plot annotations
   Tsys_str = strtrim(round(Tsys * 10.0) / 10.0, 2)
   Tsys_str = strmid(Tsys_str, 0, strpos(Tsys_str, '.') + 2)
   Bt_str = strtrim(round(Bt[i] / 1e6 * 1000.0) / 1000.0, 2)
   Bt_str = strmid(Bt_str, 0, strpos(Bt_str, '.') + 4)
   tint_str = strtrim(round(tint_actual), 2) + ' s'
   if (fbinsize ge 1000.0) then begin
      fbinsize_str = strtrim(round(fbinsize / 1000.0 * 1000.0) / 1000.0, 2)
      fbinsize_str = strmid(fbinsize_str, 0, strpos(fbinsize_str, '.') + 4) + ' kHz'
   endif else begin
      fbinsize_str = strtrim(round(fbinsize * 1000.0) / 1000.0, 2)
      fbinsize_str = strmid(fbinsize_str, 0, strpos(fbinsize_str, '.') + 4) + ' Hz'
   endelse


   ;---Plot results
   output_dir = 'plots/'
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL

   xstart = 600
   xstop = 790
   y1 = 23.2
   y2 = 25.2
   y3 = 27.2
   if (i eq 0) then begin
      window, 1, xsize=1000, ysize=600
      plot, FreqAxis, SNR, $
            /xlog, $
            title = 'Ratty RFI-to-Sigma Ratio', $
            xtitle = 'Frequency [MHz]', $
            ytitle = 'RFI-to-Sigma [dB]', $
            xstyle = 1, ystyle = 1, $
            xmargin = [9, 5], ymargin = [5.5, 4], $
            charsize = 1.18, $
            xrange = [500, fmax], $
            yrange = [0, 30]

      plots, [xstart, xstop], [y1, y1]
      xyouts, 820, y1-0.2, 't!Di!N = ' + tint_str, charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif

   if (i eq 1) then begin
      oplot, FreqAxis, SNR, $
             linestyle=4

      plots, [xstart, xstop], [y2, y2], linestyle=4
      xyouts, 820, y2-0.2, 't!Di!N = ' + tint_str, charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif

   if (i eq 2) then begin
      oplot, FreqAxis, SNR, $
             linestyle=2

      oplot, [fmin, fmax], [8.45, 8.45], $
             linestyle=1

      plots, [xstart, xstop], [y3, y3], linestyle=2
      xyouts, 820, y3-0.2, 't!Di!N = ' + tint_str, charsize=1.4 ;, t!Di!N = ' + tint_str + ', B!Dch!N = ' + fbinsize_str
   endif


   index = where(FreqAxis ge 520.0)
   print, 'Expected RFI-to-Sigma Ratio at 520 MHz = ' + strtrim(SNR[index[0]], 2) + ' [dB]'
;   print, 'Expected RFI-to-Noisefloor Ratio at 520 MHz = ' + strtrim(RFI_to_noisefloor[index[0]], 2) + ' [dB]'
   print, ''


endfor

write_png, output_dir + 'Ratty RFI-to-Sigma Ratio.png', tvrd(true=1)





End

