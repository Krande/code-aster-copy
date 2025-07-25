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
subroutine spephy(ioptch, intphy, intmod, nomu, table, &
                  freq, cham, specmr, specmi, disc, &
                  nnoe, nomcmp, nuor, nbmr, nbn, &
                  imod1, nbpf, nbm, ivitef)
    implicit none
!  RESTITUTION SUR BASE PHYSIQUE D'UNE TABL_INTSP DE REPONSE MODALE
!  LA BASE MODALE EST DEFINIE PAR UN CONCEPT MELASFLU
!  LE CONCEPT PRODUIT EST UNE TABL_INTSP
!  LE CONCEPT TABL_INTSP SE COMPOSE :
!        D'UNE STRUCTURE TABLE QUI POINTE SUR UNE TABLE DE FONCTIONS
!        COMPLEXES
!-----------------------------------------------------------------------
! IN  : IOPTCH : INDICE DONNANT LA NATURE DES INTERSPECTRES A CALCULER
!       IOPTCH = 1 : INTERSPECTRES DE DEPLACEMENTS
!       IOPTCH = 2 : INTERSPECTRES DE VITESSES
!       IOPTCH = 3 : INTERSPECTRES D' ACCELERATIONS
!       IOPTCH = 4 : INTERSPECTRES DE CONTRAINTES
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
! IN  : TABLE  : NOM UTILISATEUR DU CONCEPT TABL_INTSP DE REPONSE
!                MODALE : DONNEE DU CALCUL
! IN  : FREQ   : CARACT. MODALES DE LA BASE DE CONCEPT MELASFLU
! IN  : CHAM   : CHAMP DE GRANDEURS MODALES AUX NOEUDS DE REPONSE
! IN  : SPECMR : VECTEUR DE TRAVAIL
! IN  : SPECMI : VECTEUR DE TRAVAIL
! I/O : DISC   : DISCRETISATION FREQUENTIELLE POUR CHAQUE VITESSE
! IN  : NNOE   : LISTE DES NOEUDS OU LA REPONSE EST CALCULEE
! IN  : NUOR   : LISTE DES NUMEROS D'ORDRE DES MODES PRIS EN COMPTE
! IN  : NBMR   : NBR. DE MODES PRIS EN COMPTE
! IN  : NBN    : NBR. DE NOEUDS DE REPONSE
! IN  : NBFO1  : NBR. DE SPECTRES D EXCITATION DEFINIS
! IN  : IMOD1  : INDICE DU PREMIER MODE PRIS EN COMPTE DANS LA BASE DE
!                CONCEPT MELASFLU
! IN  : NBPF   : NBR. DE POINTS DE LA DISCRETISATION FREQUENTIELLE
! IN  : NBM    : NBR. DE MODES DE LA BASE DE CONCEPT MELASFLU
! IN  : IVITEF : INDICE VITESSE DE FLUIDE
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: intphy, intmod, exiind
    integer(kind=8) :: ioptch, nbmr, nbn, imod1, nbpf, nbm, ivitef
    integer(kind=8) :: nuor(nbmr), lnumi, lnumj, ij
    real(kind=8) :: cham(nbn, nbmr), specmr(nbpf, *), specmi(nbpf, *)
    real(kind=8) :: disc(*), freq(2, nbm, *)
    character(len=8) :: nomu, table, nnoe(nbn), nomcmp
!
    integer(kind=8) :: ival(3)
    integer(kind=8) :: nbabs, ispec, mxval, mxvals, lnoei, lnoej, lcmpi, lcmpj
    real(kind=8) :: pi
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, i2, ideb, idebm, idebn, if1, ifon
    integer(kind=8) :: il, im1, im2, imi, imj, ini, inj, lfreq, lfreqs
    integer(kind=8) :: isj, ism, iv
!
    real(kind=8) :: speci, specr
!
    character(len=24) :: kval(5), valk(2), chvals, chfreq
    character(len=24) :: chnumi, chnumj, chvale, chnoei, chnoej, chcmpi, chcmpj
!
!-----------------------------------------------------------------------
    call jemarq()
    pi = r8pi()
!
    chnumi = table//'.NUMI'
    chnumj = table//'.NUMJ'
    chfreq = table//'.DISC'
    chvale = table//'.VALE'
    call jeveuo(chnumi, 'L', lnumi)
    call jeveuo(chnumj, 'L', lnumj)
    call jeveuo(chfreq, 'L', lfreq)
    call jelira(chnumi, 'LONMAX', mxval)
!
!
! --- POUR CHAQUE PAS EN VITESSE ON CALCULE L INTERSPECTRE DE REPONSE
!
    iv = ivitef
!
    ival(1) = iv
!
! ---   TEST POUR DETECTER UN EVENTUEL PROBLEME DE CONVERGENCE EN AMONT
! ---   DANS L'OPERATEUR CALC_FLUI_STRU POUR CALCULER LES PARAMETRES
! ---   MODAUX A LA VITESSE D'ECOULEMENT CONSIDEREE
! ---   DANS CE CAS LES INTERSPECTRES DE REPONSE MODALE N'ONT PAS ETE
! ---   CALCULES PAR L'OPERATEUR DYNA_SPEC_MODAL
! ---   => ON PASSE A LA VITESSE SUIVANTE
!
    ival(2) = nuor(1)
    ival(3) = nuor(1)
    exiind = .false.
    do i1 = 1, mxval
        if ((zi(lnumi-1+i1) .eq. ival(2)) .and. (zi(lnumj-1+i1) .eq. ival(3))) exiind = &
            .true.
    end do
!
    if (.not. exiind) goto 20
!
!     --- RECUPERATION DES FONCTIONS (SPECTRES) ET STOCKAGE DANS
!     ---          SPECMR,SPECMI
!
    do imj = 1, nbmr
!
        ival(2) = nuor(imj)
!
        ideb = imj
        if (intmod) ideb = 1
!
        do imi = ideb, imj
!
            ival(3) = nuor(imi)
!
            exiind = .false.
            do i1 = 1, mxval
                if ((zi(lnumi-1+i1) .eq. ival(2)) .and. (zi(lnumj-1+i1) .eq. ival(3))) then
                    exiind = .true.
                    call jeveuo(jexnum(chvale, i1), 'L', ifon)
                end if
            end do
!
            if (.not. exiind) then
                valk(1) (1:10) = 'INTE_SPEC'
                valk(2) (1:8) = table
                call utmess('F', 'MODELISA2_91', nk=2, valk=valk)
            end if
!
            isj = (imj*(imj-1))/2+imi
            if (isj .eq. 1) then
                do if1 = 1, nbpf
                    disc(if1) = zr(lfreq+(if1-1))
                end do
            end if
!
            do if1 = 1, nbpf
                if (ival(2) .eq. ival(3)) then
                    specmr(if1, isj) = zr(ifon+(if1-1))
                    specmi(if1, isj) = 0.d0
                else
                    specmr(if1, isj) = zr(ifon+(if1-1)*2)
                    specmi(if1, isj) = zr(ifon+(if1-1)*2+1)
                end if
            end do
        end do
!
    end do
!
!    --- CREATION ET REMPLISSAGE DES FONCTIONS - SPECTRES REPONSES
!
    chnoei = nomu//'.NOEI'
    chnoej = nomu//'.NOEJ'
    chcmpi = nomu//'.CMPI'
    chcmpj = nomu//'.CMPJ'
    chvals = nomu//'.VALE'
    call wkvect(nomu//'.DISC', 'G V R', nbpf, lfreqs)
    do il = 1, nbpf
        zr(lfreqs+il-1) = disc(il)
    end do
!
    if (intphy) then
        mxvals = nbn*(nbn+1)/2
    else
        mxvals = nbn
    end if
!
    call wkvect(chnoei, 'G V K8', mxvals, lnoei)
    call wkvect(chnoej, 'G V K8', mxvals, lnoej)
    call wkvect(chcmpi, 'G V K8', mxvals, lcmpi)
    call wkvect(chcmpj, 'G V K8', mxvals, lcmpj)
    call jecrec(chvals, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mxvals)
!
    ij = 0
    do inj = 1, nbn
!
        kval(3) = nnoe(inj)
        kval(4) = nomcmp
!
        idebn = inj
        if (intphy) idebn = 1
!
        do ini = idebn, inj
!
            ij = ij+1
            kval(1) = nnoe(ini)
            kval(2) = nomcmp
!
            zk8(lnoei-1+ij) = kval(1)
            zk8(lnoej-1+ij) = kval(3)
            zk8(lcmpi-1+ij) = kval(2)
            zk8(lcmpj-1+ij) = kval(4)
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
                do im2 = 1, nbmr
!
                    idebm = im2
                    if (intmod) idebm = 1
!
                    do im1 = idebm, im2
                        i1 = imod1+im1-1
                        i2 = imod1+im2-1
                        ism = (im2*(im2-1))/2+im1
!
                        if (im1 .eq. im2) then
!                 --------------------
!
                            if (ioptch .eq. 1 .or. ioptch .eq. 4) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)
!
                            else if (ioptch .eq. 2) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*freq(1, i1, iv)*freq(&
                                        &1, i2, iv)*specmr(il, ism)
!
                            else if (ioptch .eq. 3) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*freq(1, i1, iv)*freq(&
                                        &1, i2, iv)*freq(1, i1, iv)*freq(1, i2, iv)*specmr(il, ism)
!
                            end if
!
                        else
!                 ----
!
                            if (ioptch .eq. 1 .or. ioptch .eq. 4) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)+ch&
                                        &am(ini, im2)*cham(inj, im1)*specmr(il, ism)
                                speci = speci+cham(ini, im1)*cham(inj, im2)*specmi(il, ism)-ch&
                                        &am(ini, im2)*cham(inj, im1)*specmi(il, ism)
!
                            else if (ioptch .eq. 2) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)*freq&
                                        &(1, i1, iv)*freq(1, i2, iv)+cham(ini, im2)*cham(inj, im1)*&
                                        &freq(1, i1, iv)*freq(1, i2, iv)*specmr(il, ism)
                                speci = speci+cham(ini, im1)*cham(inj, im2)*specmi(il, ism)*freq&
                                        &(1, i1, iv)*freq(1, i2, iv)-cham(ini, im2)*cham(inj, im1)*&
                                        &freq(1, i1, iv)*freq(1, i2, iv)*specmi(il, ism)
!
                            else if (ioptch .eq. 3) then
                                specr = specr+cham(ini, im1)*cham(inj, im2)*specmr(il, ism)*freq&
                                      &(1, i1, iv)*freq(1, i2, iv)*freq(1, i1, iv)*freq(1, i2, iv)+&
                                      & cham(ini, im2)*cham(inj, im1)*freq(1, i1, iv)*freq(1, i2, i&
                                        &v)*specmr(il, ism)*freq(1, i1, iv)*freq(1, i2, iv)
                                speci = speci+cham(ini, im1)*cham(inj, im2)*specmi(il, ism)*freq&
                                      &(1, i1, iv)*freq(1, i2, iv)*freq(1, i1, iv)*freq(1, i2, iv)-&
                                      & cham(ini, im2)*cham(inj, im1)*freq(1, i1, iv)*freq(1, i2, i&
                                        &v)*specmi(il, ism)*freq(1, i1, iv)*freq(1, i2, iv)
!
                            end if
!
                        end if
!                 -----
!
                    end do
                end do
!
                if (ioptch .eq. 2) then
                    specr = specr*4.d0*pi*pi
                    speci = speci*4.d0*pi*pi
                else if (ioptch .eq. 3) then
                    specr = specr*16.d0*pi*pi*pi*pi
                    speci = speci*16.d0*pi*pi*pi*pi
                end if
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
20  continue
!
    call jedema()
end subroutine
