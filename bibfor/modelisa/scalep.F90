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
subroutine scalep(spectr, noma, base, nuor, nbm, &
                  imodi, nbmr, nbexcp, ltable, iaxe, &
                  scal)
    implicit none
!     PROJECTION D'UN SPECTRE D'EXCITATION TURBULENTE LOCALISEE SUR UNE
!     BASE MODALE PERTURBEE PAR COUPLAGE FLUIDE-STRUCTURE
!     CALCUL DES PRODUITS SCALAIRES PHII(XK).NK ET PHII'(XM).XM
!     APPELANT : SPECEP
!-----------------------------------------------------------------------
! IN  : SPECTR : NOM DU CONCEPT SPECTRE
! IN  : NOMA   : NOM DU CONCEPT MAILLAGE
! IN  : BASE   : NOM DU CONCEPT MELASFLU
! IN  : NUOR   : NUMEROS D'ORDRE DES MODES DU CONCEPT MELASFLU
! IN  : NBM    : NOMBRE DE MODES DU CONCEPT MELASFLU
! IN  : IMODI  : INDICE DU PREMIER MODE PRIS EN COMPTE
! IN  : NBMR   : NOMBRE DE MODES PRIS EN COMPTE
! IN  : NBEXCP : NOMBRE D'EXCITATIONS PONCTUELLES
! IN  : LTABLE : BOOLEEN, DONNE LE CAS DE CALCUL
!       LTABLE = .TRUE.  => SPECTRE DEFINI PAR L'UTILISATEUR
!       LTABLE = .FALSE. => RECOURS A UN SPECTRE GRAPPE2 PREDEFINI
! IN  : IAXE   : ENTIER DONNANT L'AXE DIRECTEUR DE LA POUTRE
! OUT : SCAL   : TABLEAU DES PRODUITS SCALAIRES (NBEXCP,NBMR)
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nbm, nuor(nbm), imodi, nbmr, nbexcp, iaxe
    aster_logical :: ltable
    character(len=8) :: noma
    character(len=19) :: spectr, base
    real(kind=8) :: scal(nbexcp, nbmr)
!
    real(kind=8) :: dgrd
    character(len=24) :: spnnoe, spvare, spvate, chvale
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: idec, idefm, idir1, idir2, iex, imod, imodf
    integer(kind=8) :: imr, inat, inatu, inuno, irot1, irot2, ispno
    integer(kind=8) :: ispre, ispte, iteta, iv, ivale, nuno
    real(kind=8) :: theta
!-----------------------------------------------------------------------
    call jemarq()
!
! --- 1.INITIALISATIONS
!
    dgrd = r8dgrd()
    imodf = imodi+nbmr-1
    iv = 1
    if (iaxe .eq. 1) then
        idir1 = 2
        idir2 = 3
        irot1 = 6
        irot2 = 5
    else if (iaxe .eq. 2) then
        idir1 = 3
        idir2 = 1
        irot1 = 4
        irot2 = 6
    else
        idir1 = 1
        idir2 = 2
        irot1 = 5
        irot2 = 4
    end if
!
! --- 2.RECUPERATIONS SIMULTANEES
!       - DES NUMEROS DES NOEUDS D'APPLICATION DES EXCITATIONS
!       - DES DIRECTIONS D'APPLICATION DE L'EXCITATION EN CHAQUE NOEUD
!       - DE LA NATURE DE L'EXCITATION EN CHAQUE NOEUD
!
    spnnoe = spectr//'.NNOE'
    call jeveuo(spnnoe, 'L', ispno)
!
    call wkvect('&&SCALEP.TEMP.NUNO', 'V V I', nbexcp, inuno)
    call wkvect('&&SCALEP.TEMP.TETA', 'V V R', 2*nbexcp, iteta)
    call wkvect('&&SCALEP.TEMP.NATU', 'V V K8', nbexcp, inatu)
!
    if (ltable) then
!
        spvare = spectr//'.VARE'
        spvate = spectr//'.VATE'
        call jeveuo(spvare, 'L', ispre)
        call jeveuo(spvate, 'L', ispte)
!
        do iex = 1, nbexcp
            zi(inuno+iex-1) = char8_to_int(zk8(ispno+iex-1))
            theta = zr(ispre+iex-1)*dgrd
            zr(iteta+2*(iex-1)) = dble(cos(theta))
            zr(iteta+2*(iex-1)+1) = dble(sin(theta))
            zk8(inatu+iex-1) = zk16(ispte+4+iex-1) (1:8)
        end do
!
    else
!
        nuno = char8_to_int(zk8(ispno))
        do iex = 1, nbexcp
            zi(inuno+iex-1) = nuno
            zr(iteta+2*(iex-1)) = 1.d0
            zr(iteta+2*(iex-1)+1) = 1.d0
            inat = iex-int(iex/2)*2
            if (inat .eq. 0) then
                zk8(inatu+iex-1) = 'MOMENT'
            else
                zk8(inatu+iex-1) = 'FORCE'
            end if
        end do
!
    end if
!
! --- 3.RECUPERATION DES DEFORMEES MODALES OU/ET DES DERIVEES DES
! ---   DEFORMEES MODALES AUX NOEUDS D'APPLICATION
!
    call wkvect('&&SCALEP.TEMP.DEFM', 'V V R', 2*nbexcp*nbmr, idefm)
!
    do imod = imodi, imodf
!
        write (chvale, '(A8,A5,2I3.3,A5)') base(1:8), '.C01.', nuor(imod), &
            iv, '.VALE'
        call jeveuo(chvale, 'L', ivale)
        imr = imod-imodi+1
!
        do iex = 1, nbexcp
            idec = 2*nbexcp*(imr-1)+2*(iex-1)
            if (zk8(inatu+iex-1) (1:1) .eq. 'F') then
                zr(idefm+idec) = zr(ivale+6*(zi(inuno+iex-1)-1)+idir1-1)
                zr(idefm+idec+1) = zr(ivale+6*(zi(inuno+iex-1)-1)+idir2-1)
            else
                zr(idefm+idec) = zr(ivale+6*(zi(inuno+iex-1)-1)+irot1-1)
                zr(idefm+idec+1) = zr(ivale+6*(zi(inuno+iex-1)-1)+irot2-1)
            end if
        end do
!
        call jelibe(chvale)
!
    end do
!
! --- 4.CALCUL DES PRODUITS SCALAIRES
!
    do imr = 1, nbmr
        do iex = 1, nbexcp
            idec = 2*nbexcp*(imr-1)+2*(iex-1)
            scal(iex, imr) = zr(idefm+idec)*zr(iteta+2*(iex-1))+zr(idefm+idec+1)*zr(iteta+2*(iex&
                            &-1)+1)
        end do
    end do
!
    call jedetr('&&SCALEP.TEMP.NUNO')
    call jedetr('&&SCALEP.TEMP.TETA')
    call jedetr('&&SCALEP.TEMP.NATU')
    call jedetr('&&SCALEP.TEMP.DEFM')
    call jedema()
end subroutine
