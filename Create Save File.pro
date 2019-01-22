
;**************************************************************************************************
; Author:       R. Lord
; Organisation: SKA SA
; Ver:          1.0
;**************************************************************************************************


Resolve_Routine, 'python_batch_processing'
Resolve_Routine, 'Subroutines'
Resolve_All
save, /routines, filename='python_batch_processing.sav'

Resolve_Routine, 'python_stitch_spectra'
Resolve_Routine, 'Subroutines'
Resolve_All
save, /routines, filename='python_stitch_spectra.sav'


End


