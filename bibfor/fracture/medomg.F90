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

subroutine medomg(result, numord, modele, mate, mateco, lischa)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lisccm.h"
#include "asterfort/liscnv.h"
#include "asterfort/liscom.h"
#include "asterfort/lislec.h"
#include "asterfort/lisnnb.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rslesd.h"
!
!
    integer(kind=8) :: numord
    character(len=8) :: modele, result
    character(len=24) :: mate, mateco
    character(len=19) :: lischa
!
! ----------------------------------------------------------------------
!
!  OPERATEUR CALC_G
!
!  SAISIE ET VERIFICATION DE LA COHERENCE DES DONNEES MECANIQUES
!  DU PROBLEME
!
! ----------------------------------------------------------------------
!
! IN  RESULT : NOM DE LA SD RESULTAT
! IN  NUMORD : NUMERO D'ORDRE DANS SD RESULTAT
! OUT MODELE : NOM DU MODELE
! OUT MATE   : MATERIAU CODE
! OUT LISCHA : LISTE DES CHARGES
!
! ----------------------------------------------------------------------
!
    character(len=8) :: materi
    character(len=16) :: phenom, motfac, nomcmd
    character(len=19) :: lisold
    integer(kind=8) :: iexcit, nbchar
    character(len=1) :: base, codarr
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    materi = ' '
    modele = ' '
    mateco = ' '
    nomcmd = 'CALC_G'
    phenom = 'MECANIQUE'
    motfac = 'EXCIT'
    base = 'V'
    codarr = 'F'
!
! - SUPPRESSION ANCIENNES LISTE_CHARGES
!
    call detrsd(' ', lischa)
!
! - RECUPERATION MODELE, MATERIAU, CARA_ELEM ET LISCHA DANS SD RESU
!
    call rslesd(result, numord, &
                model_=modele, materi_=materi, &
                list_load_=lisold, iexcit_=iexcit)
!
! - CODAGE DU MATERIAU
!
    if (materi .ne. ' ') call rcmfmc(materi, mateco, l_ther_=ASTER_FALSE)
    mate = materi
!
! - ON PREND LE CHARGEMENT DANS LA SD
!
    if (iexcit .eq. 0) then
        call liscnv(phenom, base, lisold, lischa)
    end if
!
! - ON PREND LE CHARGEMENT DONNE PAR L'UTILISATEUR
!
    if (iexcit .eq. 1) then
        call lislec(motfac, phenom, base, lischa)
    end if
!
    call lisnnb(lischa, nbchar)
!
! - VERIFICATION DE LA COHERENCE DES MODELES
!
    if (nbchar .ne. 0) call liscom(modele, codarr, lischa)
!
! - VERIFICATION CHARGE PRISE EN COMPTE DANS CALC_G ?
!
    if (nbchar .ne. 0) call lisccm(nomcmd, 'F', lischa)
!
    call jedema()
end subroutine
