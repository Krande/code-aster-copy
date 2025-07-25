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
subroutine lisccc(nomcmd, motclc, nbauth, nbnaut, mclaut)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/iscode.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lisdef.h"
    character(len=16) :: nomcmd
    integer(kind=8) :: motclc(2)
    integer(kind=8) :: nbnaut, nbauth, mclaut(2)
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (LISTE_CHARGES)
!
! VERIFICATION COMPATIBILITE CHARGE/COMMANDE - SOUS
!
! ----------------------------------------------------------------------
!
!
! IN  NOMCMD : NOM DE LA COMMANDE
! IN  MOTCLC : CODE (ENTIER CODE) CONTENANT LES MOTS-CLEFS
! OUT NBAUTH : NOMBRE DE MOTS-CLEFS AUTORISES DANS CETTE COMMANDE
! OUT NBNAUT : NOMBRE DE MOTS-CLEFS NON AUTORISES DANS CETTE COMMANDE
! OUT MCLAUT : CODE (ENTIER CODE) CONTENANT LES MOTS-CLEFS AUTORISES
!              DANS CETTE COMMANDE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: zbgdlh, zbgccg
    parameter(zbgdlh=17, zbgccg=9)
    character(len=16) :: autdlh(zbgdlh), autccg(zbgccg)
!
    integer(kind=8) :: tabcox(60), tabaut(60)
    integer(kind=8) :: nbtota, nbgcmd
    integer(kind=8) :: ipose, iposit(2), iauth, ibid
    character(len=8) :: k8bid
    character(len=16) :: motclf
    aster_logical :: lfind
!
! --- DYNA_VIBRA//HARM/GENE
!
    data autdlh/&
     &     'DIRI_DUAL', 'FORCE_NODALE', 'EPSI_INIT',&
     &     'PRES_REP', 'FLUX_THM_REP', 'PESANTEUR',&
     &     'ROTATION', 'FORCE_CONTOUR', 'FORCE_INTERNE#3D',&
     &     'FORCE_INTERNE#2D', 'FORCE_ARETE', 'FORCE_FACE',&
     &     'FORCE_POUTRE', 'FORCE_COQUE#3D', 'FORCE_COQUE#2D',&
     &     'VECT_ASSE', 'VECT_ASSE_GENE'/
!
! --- CALC_G
!
    data autccg/&
     &     'DIRI_DUAL', 'EPSI_INIT', 'PRES_REP',&
     &     'PESANTEUR', 'ROTATION', 'FORCE_CONTOUR',&
     &     'FORCE_INTERNE#3D', 'FORCE_INTERNE#2D', 'FORCE_FACE'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    nbnaut = 0
    nbauth = 0
    nbtota = 0
    do ipose = 1, 60
        tabaut(ipose) = 0
    end do
!
! --- DECODAGE DU CHARGEMENT
!
    call isdeco(motclc, tabcox, 60)
!
! --- SELECTION DATA DE LA COMMANDE ACTIVE
!
    if (nomcmd .eq. 'DYNA_VIBRA') then
        nbgcmd = zbgdlh
    else if (nomcmd .eq. 'CALC_G') then
        nbgcmd = zbgccg
    else
        ASSERT(.false.)
    end if
!
! --- BOUCLE SUR LES MOTS-CLEFS ACTIFS DANS LA CHARGE
!
    do ipose = 1, 60
        if (tabcox(ipose) .eq. 1) then
            lfind = .false.
            nbtota = nbtota+1
!
! --------- BOUCLE SUR LES MOTS-CLEFS AUTORISES
!
            do iauth = 1, nbgcmd
                if (nomcmd .eq. 'DYNA_VIBRA') then
                    motclf = autdlh(iauth)
                    call lisdef('POSM', motclf, ibid, k8bid, iposit)
                else if (nomcmd .eq. 'CALC_G') then
                    motclf = autccg(iauth)
                    call lisdef('POSM', motclf, ibid, k8bid, iposit)
                else
                    ASSERT(.false.)
                end if
                ASSERT(iposit(1) .ne. 0)
                if (iposit(1) .eq. ipose) lfind = .true.
            end do
            if (lfind) then
                nbauth = nbauth+1
                tabaut(ipose) = 1
            end if
        else
            tabaut(ipose) = 0
        end if
    end do
!
! --- CODAGE MOT-CLEFS AUTORISES
!
    call iscode(tabaut, mclaut, 60)
!
    call jedema()
end subroutine
