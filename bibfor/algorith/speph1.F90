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
subroutine speph1(intphy, intmod, nomu, cham, specmr, &
                  specmi, nnoe, nomcmp, nbmode, nbn, &
                  nbpf)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: intphy, intmod
    integer(kind=8) :: nbmode, nbn, nbpf
    real(kind=8) :: cham(nbn, *), specmr(nbpf, *), specmi(nbpf, *)
    character(len=8) :: nomu, nnoe(*), nomcmp(*)
!  RESTITUTION SUR BASE PHYSIQUE D'UNE TABL_INTSP DE REPONSE MODALE
!  LA BASE MODALE EST DEFINIE PAR UN CONCEPT MODE_MECA
!-----------------------------------------------------------------------
! IN  : INTPHY : BOOLEEN
!                CARACTERISE LE CONTENU DE LA TABLE D'INTERSPECTRES DE
!                REPONSE PHYSIQUE A CALCULER
!       INTPHY = .TRUE.  TOUS LES INTERSPECTRES SERONT CALCULES
!       INTPHY = .FALSE. SEULS LES AUTOSPECTRES SERONT CALCULES
! IN  : INTMOD : BOOLEEN
!                CARACTERISE LE CONTENU DE LA TABLE D'INTERSPECTRES DE
!                REPONSE MODALE (DONNEE DU CALCUL)
!       INTMOD = .TRUE.  TOUS LES INTERSPECTRES ONT ETE CALCULES
!       INTMOD = .FALSE. SEULS LES AUTOSPECTRES ONT ETE CALCULES
! IN  : NOMU   : NOM UTILISATEUR DU CONCEPT TABL_INTSP DE REPONSE
!                PHYSIQUE : A PRODUIRE
! IN  : CHAM   : CHAMP DE GRANDEURS MODALES AUX NOEUDS DE REPONSE
! IN  : SPECMR : VECTEUR DE TRAVAIL
! IN  : SPECMI : VECTEUR DE TRAVAIL
! IN  : NNOE   : LISTE DES NOEUDS OU LA REPONSE EST CALCULEE
! IN  : NBMODE : NBR. DE MODES PRIS EN COMPTE
! IN  : NBN    : NBR. DE NOEUDS DE REPONSE
! IN  : NBPF   : NBR. DE POINTS DE LA DISCRETISATION FREQUENTIELLE
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nbpar, inj, idebn, ini, il, im2, idebm, im1, ism
    integer(kind=8) :: nbabs, ispec, mxval, lnoei, lnoej, lcmpi, lcmpj, ij
    parameter(nbpar=5)
    real(kind=8) :: specr, speci
    character(len=24) :: kval(nbpar)
    character(len=24) :: chnoei, chnoej, chcmpi, chcmpj, chvals
!
!-----------------------------------------------------------------------
    call jemarq()
!
!    --- CREATION ET REMPLISSAGE DES FONCTIONS - SPECTRES REPONSES
!
    chnoei = nomu//'.NOEI'
    chnoej = nomu//'.NOEJ'
    chcmpi = nomu//'.CMPI'
    chcmpj = nomu//'.CMPJ'
    chvals = nomu//'.VALE'
!
    if (intphy) then
        mxval = nbn*(nbn+1)/2
    else
        mxval = nbn
    end if
!
    call wkvect(chnoei, 'G V K8', mxval, lnoei)
    call wkvect(chnoej, 'G V K8', mxval, lnoej)
    call wkvect(chcmpi, 'G V K8', mxval, lcmpi)
    call wkvect(chcmpj, 'G V K8', mxval, lcmpj)
    call jecrec(chvals, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mxval)
!
    ij = 0
    do inj = 1, nbn
!
        kval(3) = nnoe(inj)
        kval(4) = nomcmp(inj)
!
        idebn = inj
        if (intphy) idebn = 1
!
        do ini = idebn, inj
!
            ij = ij+1
            kval(1) = nnoe(ini)
            kval(2) = nomcmp(ini)
!
            zk8(lnoei-1+ij) = kval(1) (1:8)
            zk8(lnoej-1+ij) = kval(3) (1:8)
            zk8(lcmpi-1+ij) = kval(2) (1:8)
            zk8(lcmpj-1+ij) = kval(4) (1:8)
!
            if ((kval(1) .eq. kval(3)) .and. (kval(2) .eq. kval(4))) then
                nbabs = nbpf
            else
                nbabs = 2*nbpf
            end if
!
            call jecroc(jexnum(chvals, ij))
            call jeecra(jexnum(chvals, ij), 'LONMAX', nbabs)
            call jeecra(jexnum(chvals, ij), 'LONUTI', nbabs)
            call jeveuo(jexnum(chvals, ij), 'E', ispec)
!
            do il = 1, nbpf
!
                specr = 0.d0
                speci = 0.d0
!
                do im2 = 1, nbmode
!
                    idebm = im2
                    if (intmod) idebm = 1
!
                    do im1 = idebm, im2
                        ism = (im2*(im2-1))/2+im1
!
                        if (im1 .eq. im2) then
!                 --------------------
!
                            specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)
!
!
                        else
!                 ----
!
                            specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)+cham(i&
                                    &ni, im2)*cham(inj, im1)*specmr(il, ism)
                            speci = speci+cham(ini, im1)*cham(inj, im2)*specmi(il, ism)-cham(i&
                                    &ni, im2)*cham(inj, im1)*specmi(il, ism)
!
                        end if
!                 -----
!
                    end do
                end do
!
                if ((kval(1) .eq. kval(3)) .and. (kval(2) .eq. kval(4))) then
                    zr(ispec-1+il) = specr
                else
                    zr(ispec+2*(il-1)) = specr
                    zr(ispec+2*(il-1)+1) = speci
                end if
            end do
!
        end do
!
    end do
!
    call jedema()
end subroutine
