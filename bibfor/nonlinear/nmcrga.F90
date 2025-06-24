! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine nmcrga(sderro)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24), intent(in) :: sderro
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD GESTION ALGO)
!
! CREATION DE LA SD
!
! --------------------------------------------------------------------------------------------------
!
! NB: LA SD S'APPELLE SDERRO
!
! IN  SDERRO : SD ERREUR
!
! --------------------------------------------------------------------------------------------------
!
! - Name of events
    character(len=16), parameter :: eventName(ZEVEN) = (/'ERRE_INTE', 'INTE_NPHY', 'DIVE_DEBO', &
                                                         'INTE_BORN', 'ERRE_NPHY', &
                                                         'ERRE_PILO', 'CONV_PILO', 'ERRE_FACS', &
                                                         'ERRE_FACT', 'ERRE_CTD1', 'ERRE_CTD2', &
                                                         'ERRE_TIMN', 'ERRE_TIMP', 'ERRE_EXCP', &
                                                         'ITER_MAXI', &
                                                         'DIVE_RESI', 'RESI_MAXR', 'RESI_MAXN', &
                                                         'CRIT_STAB', 'DIVE_FIXG', 'RESI_MAXI', &
                                                         'DIVE_FIXF', 'DIVE_FIXC', 'ERRE_CTCG', &
                                                         'ERRE_CTCF', 'ERRE_CTCC', 'DIVE_FROT', &
                                                         'DIVE_GEOM', 'DIVE_RELA', 'DIVE_MAXI', &
                                                         'DIVE_REFE', 'DIVE_COMP', 'DIVE_CTCC', &
                                                         'SOLV_ITMX', 'DIVE_HROM', 'DIVE_PENE', &
                                                         'ERRE_APPA'/)
! - Return code (name)
    character(len=8), parameter :: eventReturnCode(ZEVEN) = (/'LDC', 'LDC', 'LDC', &
                                                              'LDC', 'XXX', &
                                                              'PIL', 'PIL', 'FAC', &
                                                              'FAC', 'CTC', 'CTC', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'XXX', 'XXX', 'XXX', &
                                                              'RES', 'XXX', 'XXX', &
                                                              'XXX'/)
! - Return code (value)
    integer(kind=8), parameter :: eventReturnValue(ZEVEN) = (/1, 2, 3, &
                                                              4, 99, &
                                                              1, 2, 1, &
                                                              2, 1, 2, &
                                                              99, 99, 99, &
                                                              99, &
                                                              99, 99, 99, &
                                                              99, 99, 99, &
                                                              99, 99, 99, &
                                                              99, 99, 99, &
                                                              99, 99, 99, &
                                                              99, 99, 99, &
                                                              1, 99, 99, &
                                                              99/)
!
! --- TYPE ET NIVEAU DE DECLENCHEMENT POSSIBLES DE L'EVENEMENT
! TROIS TYPES
! EVEN  : EVENEMENT A CARACTERE PUREMENT INFORMATIF
!          -> PEUT ETRE TRAITE SI UTILISATEUR LE DEMANDE DANS
!             DEFI_LIST_INST
! ERRI_ : EVENEMENT A TRAITER IMMEDIATEMENT SI ON VEUT CONTINUER
! ERRC_ : EVENEMENT A TRAITER A CONVERGENCE
! CONV_ : EVENEMENT A TRAITER POUR DETERMINER LA CONVERGENCE
!
    character(len=16), parameter :: eventLevel(ZEVEN) = (/'ERRI_NEWT', 'ERRC_NEWT', 'CONV_NEWT', &
                                                          'EVEN     ', 'ERRI_NEWT', &
                                                          'ERRI_NEWT', 'CONV_CALC', 'ERRI_NEWT', &
                                                          'ERRI_NEWT', 'ERRI_NEWT', 'ERRI_NEWT', &
                                                          'ERRI_CALC', 'ERRI_CALC', 'ERRI_CALC', &
                                                          'ERRI_NEWT', &
                                                          'EVEN     ', 'EVEN     ', 'EVEN     ', &
                                                          'EVEN     ', 'CONV_FIXE', 'EVEN     ', &
                                                          'CONV_FIXE', 'CONV_FIXE', 'ERRI_FIXE', &
                                                          'ERRI_FIXE', 'ERRI_FIXE', 'CONV_RESI', &
                                                          'CONV_NEWT', 'CONV_RESI', 'CONV_RESI', &
                                                          'CONV_RESI', 'CONV_RESI', 'CONV_NEWT', &
                                                          'ERRI_NEWT', 'CONV_FIXE', 'CONV_RESI', &
                                                          'ERRI_NEWT'/)
!
! --- FONCTIONNALITE ACTIVE SI NECESSAIRE POUR CONVERGENCE
!
    character(len=24), parameter :: eventActiFunc(ZEVEN) = &
                                    (/'         ', '         ', '         ', &
                                      '         ', '         ', &
                                      '         ', 'PILOTAGE ', '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', '         ', '         ', &
                                      '         ', 'RESI_RELA', 'RESI_MAXI', &
                                      'RESI_REFE', 'RESI_COMP', '         ', &
                                      'LDLT_SP  ', '         ', '         ', &
                                      '         '/)
!
! --- CODE DU MESSAGE A AFFICHER
!
    character(len=24), parameter :: eventMesg(ZEVEN) = (/ &
                                    'MECANONLINE10_1 ', 'MECANONLINE10_13', '                ', &
                                    'MECANONLINE10_25', 'MECANONLINE10_13', &
                                    'MECANONLINE10_2 ', '                ', 'MECANONLINE10_6 ', &
                                    'MECANONLINE10_6 ', 'MECANONLINE10_4 ', 'MECANONLINE10_4 ', &
                                    'MECANONLINE10_7 ', 'MECANONLINE10_5 ', 'MECANONLINE10_8 ', &
                                    'MECANONLINE10_3 ', &
                                    '                ', '                ', '                ', &
                                    'MECANONLINE10_20', '                ', 'MECANONLINE10_26', &
                                    '                ', '                ', 'MECANONLINE10_9 ', &
                                    'MECANONLINE10_10', 'MECANONLINE10_11', '                ', &
                                    '                ', '                ', '                ', &
                                    '                ', '                ', '                ', &
                                    'MECANONLINE10_12', '                ', '                ', &
                                    'MECANONLINE10_14'/)
!
    integer(kind=8) :: iEvent
    character(len=24) :: eventECONJv, eventECOVJv, eventENIVJv, eventEFCTJv, eventEMSGJv
    character(len=24) :: eventCONVJv, eventEEVTJv, eventENOMJv, eventEACTJv
    integer(kind=8), pointer :: eventEACT(:) => null(), eventECOV(:) => null()
    integer(kind=8), pointer :: eventEEVT(:) => null(), eventCONV(:) => null()
    character(len=16), pointer :: eventENOM(:) => null()
    character(len=16), pointer :: eventENIV(:) => null()
    character(len=24), pointer :: eventEMSG(:) => null(), eventEFCT(:) => null()
    character(len=8), pointer :: eventECON(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - OBJETS
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventECOVJv = sderro(1:19)//'.ECOV'
    eventECONJv = sderro(1:19)//'.ECON'
    eventENIVJv = sderro(1:19)//'.ENIV'
    eventEFCTJv = sderro(1:19)//'.EFCT'
    eventEACTJv = sderro(1:19)//'.EACT'
    eventCONVJv = sderro(1:19)//'.CONV'
    eventEEVTJv = sderro(1:19)//'.EEVT'
    eventEMSGJv = sderro(1:19)//'.EMSG'
    call wkvect(eventENOMJv, 'V V K16', ZEVEN, vk16=eventENOM)
    call wkvect(eventECOVJv, 'V V I', ZEVEN, vi=eventECOV)
    call wkvect(eventECONJv, 'V V K8', ZEVEN, vk8=eventECON)
    call wkvect(eventENIVJv, 'V V K16', ZEVEN, vk16=eventENIV)
    call wkvect(eventEFCTJv, 'V V K24', ZEVEN, vk24=eventEFCT)
    call wkvect(eventEACTJv, 'V V I', ZEVEN, vi=eventEACT)
    call wkvect(eventCONVJv, 'V V I', NB_LOOP+1, vi=eventCONV)
    call wkvect(eventEEVTJv, 'V V I', 2, vi=eventEEVT)
    call wkvect(eventEMSGJv, 'V V K24', ZEVEN, vk24=eventEMSG)
!
    do iEvent = 1, ZEVEN
        eventENOM(iEvent) = eventName(iEvent)
        eventEACT(iEvent) = EVENT_IS_INACTIVE
        eventECON(iEvent) = eventReturnCode(iEvent)
        eventECOV(iEvent) = eventReturnValue(iEvent)
        eventENIV(iEvent) = eventLevel(iEvent)
        eventEFCT(iEvent) = eventActiFunc(iEvent)
        eventEMSG(iEvent) = eventMesg(iEvent)
    end do
!
    call jedema()
end subroutine
