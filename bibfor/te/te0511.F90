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

subroutine te0511(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/evala1.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
    character(len=16) :: option, nomte
! =====================================================================
!    - FONCTION REALISEE:  CALCUL DU MODULE DE RIGIDITE
!                          DE MICRO-DILTATION CONTENU
!                          AUX POINTS DE GAUSS
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! =====================================================================
! =====================================================================
    aster_logical :: logthm
    integer(kind=8) :: imate, ivarip, icontp, ilocal, ic, iv
    integer(kind=8) :: nbvari, kpg, nbsig, imat, ibid
    integer(kind=8) :: icode, iret, tabthm(3), dimmax, npgu
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano
    real(kind=8) :: module, sig(6)
    real(kind=8) :: vin(100)
    character(len=8) :: mod, alias8
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relcom
    integer(kind=8), parameter :: rindic = 1
! =====================================================================
! --- RINDIC EST LE NOMBRE DE PARAMETRE DE LOCALISATION DEFINIT -------
! --- SOUS LE MOT-CLE INDL_R DANS GRANDEUR_SIMPLE.CATA --------------
! =====================================================================

! =====================================================================
    call teattr('S', 'ALIAS8', alias8, ibid)
    if (option .eq. 'PDIL_ELGA') then
! =====================================================================
! --- VERIFICATION DE COHERENCE ---------------------------------------
! --- LE TENSEUR ACOUSTIQUE EST DEVELOPPE EN 2D UNIQUEMENT ------------
! =====================================================================
! --- CAS D'UN POST-TRAITEMENT EN MECANIQUE DRAINE --------------------
! =====================================================================
        logthm = .false.
        if ((alias8(3:5) .eq. 'DPL') .or. (alias8(3:5) .eq. 'DPS')) then
            mod(1:6) = 'D_PLAN'
            nbsig = nbsigm()
        else if (alias8(3:5) .eq. 'CPL') then
            mod(1:6) = 'C_PLAN'
            nbsig = nbsigm()
        else if (alias8(3:5) .eq. 'AX_') then
            mod(1:4) = 'AXIS'
            nbsig = nbsigm()
        else
! =====================================================================
! --- CAS D'UN POST-TRAITEMENT EN MECANIQUE THM -----------------------
! =====================================================================
            logthm = .true.
            if (alias8(3:5) .eq. 'AH2') then
                mod(1:4) = 'AXIS'
            else if ((alias8(3:5) .eq. 'DH2') .or. (alias8(3:5) .eq. 'DR1') .or. &
                     (alias8(3:5) .eq. 'DM1')) then
                mod(1:6) = 'D_PLAN'
            else
! =====================================================================
! --- CAS NON TRAITE --------------------------------------------------
! =====================================================================
                call utmess('F', 'ELEMENTS_11', sk=nomte)
            end if
        end if
! =====================================================================
! --- RECUPERATION DU ELREFE ------------------------------------------
! =====================================================================
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
! =====================================================================
! --- PARAMETRES EN ENTREE --------------------------------------------
! =====================================================================
        call jevech('PMATERC', 'L', imate)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PVARIPR', 'L', ivarip)
        if (logthm) then
! =====================================================================
! --- DANS LE CADRE THM ON FAIT UN TECACH PLUTOT QU'UN JEVECH POUR ----
! --- RECUPERER EGALEMENT LA DIMENSION DU VECTEUR QUI DIFFERE SUIVANT -
! --- LA MODELISATION THM ---------------------------------------------
! =====================================================================
            call tecach('OOO', 'PCONTPR', 'L', iret, nval=3, &
                        itab=tabthm)
            icontp = tabthm(1)
            dimmax = tabthm(2)
            npgu = tabthm(3)
! =====================================================================
! --- on teste la coherence des recuperations elrefe_info et tecach sur ----
! --- LE NOMBRE DE POINTS DE GAUSS ------------------------------------
! =====================================================================
            ASSERT(npgu .eq. npg)
            nbsig = dimmax/npg
! =====================================================================
! --- DANS LE CADRE DE LA THM ON RECUPERE DIRECTEMENT LA RELATION -----
! --- DE COMPORTEMENT DE TYPE MECANIQUE -------------------------------
! =====================================================================
            relcom = compor(MECA_NAME)
        else
            call jevech('PCONTPR', 'L', icontp)
            relcom = compor(RELA_NAME)
        end if
! =====================================================================
! --- NOMBRE DE VARIABLES INTERNES ASSOCIE A LA LOI DE COMPORTEMENT ---
! =====================================================================
        read (compor(NVAR), '(I16)') nbvari
! =====================================================================
! --- PARAMETRES EN SORTIE --------------------------------------------
! =====================================================================
        call jevech('PPDIL', 'E', ilocal)
! =====================================================================
! --- BOUCLE SUR LES POINTS DE GAUSS ----------------------------------
! =====================================================================
        do kpg = 1, npg
! =====================================================================
! --- INITIALISATIONS -------------------------------------------------
! =====================================================================
            module = 0.0d0
! =====================================================================
! --- CALCUL DU MODULE DE RIGIDITE DE MICRO-DILTATION -----------------
! =====================================================================
            do ic = 1, nbsig
                sig(ic) = zr(icontp-1+(kpg-1)*nbsig+ic)
            end do
            do iv = 1, nbvari
                vin(iv) = zr(ivarip-1+(kpg-1)*nbvari+iv)
            end do
            imat = zi(imate)
            call evala1('RIGI', kpg, 1, mod, relcom, &
                        sig, vin, imat, module, icode)
! =====================================================================
! --- SURCHARGE DE L'INDICATEUR DE LOCALISATION -----------------------
! =====================================================================
            zr(ilocal-1+1+(kpg-1)*rindic) = module
        end do
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    end if
! =====================================================================
end subroutine
