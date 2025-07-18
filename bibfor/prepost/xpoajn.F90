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
subroutine xpoajn(maxfem, ino, lsn, jdirno, prefno, &
                  nfiss, he, nnn, inn, inntot, &
                  nbnoc, nbnofi, inofi, co, iacoo2)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    character(len=2) :: prefno(4)
    character(len=8) :: maxfem
    integer(kind=8) :: jdirno, nnn, inn, inntot, nbnoc, ino
    integer(kind=8) :: nbnofi, inofi, iacoo2, nfiss, he(nfiss)
    real(kind=8) :: lsn(nfiss), co(3)
!
!            ON AJOUTE UN NOUVEAU NOEUD AU NOUVEAU MAILLAGE X-FEM
!
!   IN
!     INO   : NUMÉRO DU NOEUD OU DU POINT D'INTERSECTION
!     LSN    : LEVEL SETS NORMALES EN INO
!     JDIRNO : ADRESSE DU TABLEAU DIRNO LOCAL
!     PREFNO : PREFERENCES POUR LE NOMAGE DES NOUVELLES ENTITES
!     NFISS  : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT PARENT
!     HE     : VALEURS DE(S) FONCTION(S) HEAVISIDE SUR LE SOUS ÉLÉMENT
!     NNN    : NOMBRE DE NOUVEAU NOEUDS A CREER SUR LA MAILLE PARENT
!     INN    : COMPTEUR LOCAL DU NOMBRE DE NOUVEAUX NOEUDS CREES
!     INNTOT : COMPTEUR TOTAL DU NOMBRE DE NOUVEAUX NOEUDS CREES
!     NBNOC  : NOMBRE DE NOEUDS CLASSIQUES DU MAILLAGE FISSURE
!     NBNOFI : NOMBRE DE NOEUDS SITUES SUR LA FISSURE
!     INOFI  : LISTE DES NOEUDS SITUES SUR LA FISSURE
!     NBNOLA : NOMBRE DE NOEUDS SITUES SUR LA FISSURE AVEC DES LAGS
!     INOLA  : LISTE DES NOEUDS SITUES SUR LA FISSURE AVEC DES LAGS
!     CO     : COORDONNEES  INITIALES DE INO
!     IACOO2 : ADRESSE DES COORDONNES DES NOEUDS DU MAILLAGE FISSURE
!     DDLC :  NOMBRE DE DDLS DE CONTACT DE L'ÉLÉMENT PARENT
!   OUT
!     MAXFEM : NOM DU MAILLAGE FISSURE
!     INN    : COMPTEUR LOCAL DU NOMBRE DE NOUVEAU NOEUDS CREES
!     INNTOT : COMPTEUR TOTAL DU NOMBRE DE NOUVEAU NOEUDS CREES
!     NBNOFI : NOMBRE DE NOEUDS SITUES SUR LA FISSURE
!     INOFI  : LISTE DES NOEUDS SITUES SUR LA FISSURE
!     NBNOLA : NOMBRE DE NOEUDS SITUES SUR LA FISSURE AVEC DES LAGS
!     INOLA  : LISTE DES NOEUDS SITUES SUR LA FISSURE AVEC DES LAGS
!     IACOO2 :  ADRESSE DES COORDONNES DES NOEUDS DU MAILLAGE FISSURE
!
!
    real(kind=8) :: crilsn, minlsn
    integer(kind=8) :: j, ifiss, fiss
    character(len=2) :: nm
    character(len=6) :: chn
    character(len=8) :: valk(2)
    parameter(crilsn=1.d-4)
    aster_logical :: lpint
    data valk/'NOEUDS', 'XPOAJN'/
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- LPINT EST VRAI SI LE NOEUD DU MAILLAGE X-FEM EST SUR LA FISSURE.
! --- ON ATTACHERA DANS CE CAS CE NOEUDS AU GROUPE NFISSU
    if (ino .lt. 1000) then
        lpint = .false.
        do ifiss = 1, nfiss
            if (lsn(ifiss) .eq. 0.d0) lpint = .true.
        end do
    else if (ino .gt. 1000 .and. ino .lt. 2000) then
        lpint = .true.
    else if (ino .gt. 2000) then
        lpint = .false.
        do ifiss = 1, nfiss
            if (abs(lsn(ifiss)) .lt. crilsn) lpint = .true.
        end do
    end if
!
    if (lpint) then
        minlsn = r8maem()
        do ifiss = 1, nfiss
!     ON DETECTE LA FISSURE CORESPONDANTE AU POINT D'INTERSECTION
!     ATTENTION, IL PEUT Y AVOIR PLUSIEURS CANDIDAT AU NIV DE L'INTER
            if (abs(lsn(ifiss)) .lt. minlsn .and. he(ifiss) .ne. 0) then
                minlsn = abs(lsn(ifiss))
                fiss = ifiss
            end if
        end do
        if (he(fiss) .eq. -1) then
            nm = prefno(2)
        else
            nm = prefno(3)
        end if
    else
        nm = prefno(1)
    end if
!
!     COMPTEUR DES NOMS DES NOEUDS
    if (inntot .ge. 1291467968) then
        call utmess('F', 'XFEM_8', sk=valk(1))
    end if
    inn = inn+1
    inntot = inntot+1
    ASSERT(inn .le. nnn)
!
    zi(jdirno-1+(2+nfiss)*(inn-1)+1) = ino
    zi(jdirno-1+(2+nfiss)*(inn-1)+2) = nbnoc+inntot
    do ifiss = 1, nfiss
        zi(jdirno-1+(2+nfiss)*(inn-1)+2+ifiss) = he(ifiss)
    end do
    call codlet(inntot, 'G', chn)
!
    call jecroc(jexnom(maxfem//'.NOMNOE', nm//chn))
    do j = 1, 3
        zr(iacoo2-1+3*(nbnoc+inntot-1)+j) = co(j)
    end do
!       LISTE DES NOEUDS SUR LA FISSURE
    if (lpint) then
        nbnofi = nbnofi+1
        zi(inofi-1+nbnofi) = nbnoc+inntot
!        IF (HE(FISS).EQ.-1.AND.DDLC.GT.0) THEN
!        ATTENTION, IL FAUDRA RÉCUPÉRER DDLC AVEC XPOCMP DS XPOMAX
!       LISTE DES NOEUDS PORTANT DDLS DE CONTACT (COTÉ ESCLAVE)
!          NBNOLA=NBNOLA+1
!          ZI(INOLA-1+NBNOLA)=NBNOC+INNTOT
!        ENDIF
    end if
!
    call jedema()
end subroutine
