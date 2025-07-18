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
subroutine rc36si(noma, nbma, listma)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/codent.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ordis.h"
#include "asterfort/rc36cm.h"
#include "asterfort/rc36th.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbma, listma(*)
    character(len=8) :: noma
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!     RECUPERATION DES DONNEES DE "SITUATION"
!
! IN  : NOMA   : MAILLAGE
! IN  : NBMA   : NOMBRE DE MAILLES D'ANALYSE
! IN  : LISTMA : LISTE DES MAILLES D'ANALYSE
!     ------------------------------------------------------------------
!
    integer(kind=8) :: n1, nbsitu, iocc, jmomea, jmomeb, ii, nocc, jreth, jnbocc
    integer(kind=8) :: jnumgr, jpresa, jpresb, nbchar, jchar1, jnsitu, jcombi
    integer(kind=8) :: jpassa, ig, numpas(2), nscy, nbgr, numgr, nbsigr, jnsg
    integer(kind=8) :: nbth, jseigr, jchth, nume, nbm, nbp12, nbp23, nbp13, jsp12, jsp23
    integer(kind=8) :: jsp13, nbsg1, nbsg2, nbsg3, jsigr, vali(3), nbgrt, numg1, numg2
    integer(kind=8) :: jspas, ing, jnbvg, nbvg, ndim, numgs, nbseis
    aster_logical :: yapass
    character(len=8) :: k8b, ouinon
    character(len=16) :: motcl1, motcl2
    character(len=24) :: chmome
    integer(kind=8), pointer :: char_meca(:) => null()
    integer(kind=8), pointer :: nume_group(:) => null()
! DEB ------------------------------------------------------------------
!
    motcl1 = 'SITUATION'
    motcl2 = 'SEISME'
!
    call getfac(motcl1, nbsitu)
    call getfac(motcl2, nbseis)
!
    ndim = nbsitu+nbseis
    AS_ALLOCATE(vi=nume_group, size=ndim)
    call wkvect('&&RC32SI.SITU_GROUP', 'V V I', 2*ndim, jsigr)
!
    call wkvect('&&RC3600.SITU_NUMERO', 'V V I', ndim, jnsitu)
    call wkvect('&&RC3600.SITU_NB_OCCUR', 'V V I', 2*ndim, jnbocc)
    call wkvect('&&RC3600.SITU_PRES_A', 'V V R', nbsitu, jpresa)
    call wkvect('&&RC3600.SITU_PRES_B', 'V V R', nbsitu, jpresb)
    call wkvect('&&RC3600.SITU_COMBINABLE', 'V V L', ndim, jcombi)
    call wkvect('&&RC3600.SITU_PASSAGE', 'V V I', 2*nbsitu, jpassa)
    call wkvect('&&RC3600.SITU_MOMENT_A', 'V V K24', ndim, jmomea)
    call wkvect('&&RC3600.SITU_MOMENT_B', 'V V K24', nbsitu, jmomeb)
    call jecrec('&&RC3600.SITU_THERMIQUE', 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                ndim)
    call wkvect('&&RC3600.CHAM_THER', 'V V K24', ndim, jchth)
!
    call wkvect('&&RC32SI.PASSAGE_1_2', 'V V I', ndim, jsp12)
    call wkvect('&&RC32SI.PASSAGE_2_3', 'V V I', ndim, jsp23)
    call wkvect('&&RC32SI.PASSAGE_1_3', 'V V I', ndim, jsp13)
    call jeecra('&&RC32SI.PASSAGE_1_2', 'LONUTI', 0)
    call jeecra('&&RC32SI.PASSAGE_2_3', 'LONUTI', 0)
    call jeecra('&&RC32SI.PASSAGE_1_3', 'LONUTI', 0)
!
    nbgr = 0
    yapass = .false.
!
    do iocc = 1, nbsitu, 1
!
        call codent(iocc, 'D0', k8b)
!
! ------ LE NUMERO DE SITUATION:
!        -----------------------
        call getvis(motcl1, 'NUME_SITU', iocc=iocc, scal=zi(jnsitu+iocc-1), nbret=n1)
!
! ------ LE NOMBRE D'OCCURRENCE:
!        -----------------------
        call getvis(motcl1, 'NB_OCCUR', iocc=iocc, scal=nocc, nbret=n1)
        zi(jnbocc+2*iocc-2) = nocc
!
! ------ LES PRESSIONS:
!        --------------
        call getvr8(motcl1, 'PRES_A', iocc=iocc, scal=zr(jpresa+iocc-1), nbret=n1)
        call getvr8(motcl1, 'PRES_B', iocc=iocc, scal=zr(jpresb+iocc-1), nbret=n1)
!
! ------ LES NUMEROS DE GROUPE:
!        ----------------------
        call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            nbvg = -n1
            call wkvect('&&RC36SI.VALE_GR', 'V V I', nbvg, jnbvg)
            call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=nbvg, vect=zi(jnbvg), &
                        nbret=n1)
            do ing = 1, nbvg
                numgr = zi(jnbvg+ing-1)
                if (numgr .le. 0) then
                    call utmess('F', 'POSTRCCM_12')
                end if
                do ig = 1, nbgr
                    if (nume_group(ig) .eq. numgr) goto 21
                end do
                nbgr = nbgr+1
                nume_group(nbgr) = numgr
21              continue
            end do
            if (nbvg .eq. 1) then
                zi(jsigr+2*iocc-2) = zi(jnbvg)
                zi(jsigr+2*iocc-1) = zi(jnbvg)
            else
                zi(jsigr+2*iocc-2) = zi(jnbvg)
                zi(jsigr+2*iocc-1) = zi(jnbvg+1)
            end if
            call jedetr('&&RC36SI.VALE_GR')
        end if
!
! ------ LES NUMEROS DE PASSAGE:
!        -----------------------
        call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=2, vect=numpas, &
                        nbret=n1)
            if (numpas(1) .le. 0) then
                call utmess('F', 'POSTRCCM_12')
            end if
            if (numpas(2) .le. 0) then
                call utmess('F', 'POSTRCCM_12')
            end if
            if (numpas(1) .gt. 3) then
                call utmess('F', 'POSTRCCM_12')
            end if
            if (numpas(2) .gt. 3) then
                call utmess('F', 'POSTRCCM_12')
            end if
            yapass = .true.
            zi(jsigr+2*iocc-2) = min(numpas(1), numpas(2))
            zi(jsigr+2*iocc-1) = max(numpas(1), numpas(2))
            numgr = numpas(1)
            do ig = 1, nbgr
                if (nume_group(ig) .eq. numgr) goto 23
            end do
            nbgr = nbgr+1
            nume_group(nbgr) = numgr
23          continue
            numgr = numpas(2)
            do ig = 1, nbgr
                if (nume_group(ig) .eq. numgr) goto 25
            end do
            nbgr = nbgr+1
            nume_group(nbgr) = numgr
25          continue
        end if
!
! ------ COMBINABLE DANS SON GROUPE:
!        ---------------------------
        call getvtx(motcl1, 'COMBINABLE', iocc=iocc, scal=ouinon, nbret=n1)
        if (ouinon(1:3) .eq. 'OUI') then
            zl(jcombi+iocc-1) = .true.
        else
            zl(jcombi+iocc-1) = .false.
        end if
!
! ------ ETAT DE CHARGEMENT POUR "A":
!        ----------------------------
        call getvis(motcl1, 'CHAR_ETAT_A', iocc=iocc, nbval=0, nbret=n1)
        nbchar = -n1
        call wkvect('&&RC36SI.CHAR_ETAT', 'V V I', nbchar, jchar1)
        call getvis(motcl1, 'CHAR_ETAT_A', iocc=iocc, nbval=nbchar, vect=zi(jchar1), &
                    nbret=n1)
!
        chmome = '&&RC36SI_A'//k8b
        call rc36cm(iocc, 'A', nbma, listma, nbchar, &
                    zi(jchar1), chmome)
        zk24(jmomea+iocc-1) = chmome
        call jedetr('&&RC36SI.CHAR_ETAT')
!
! ------ ETAT DE CHARGEMENT POUR "B":
!        ----------------------------
        call getvis(motcl1, 'CHAR_ETAT_B', iocc=iocc, nbval=0, nbret=n1)
        nbchar = -n1
        AS_ALLOCATE(vi=char_meca, size=nbchar)
        call getvis(motcl1, 'CHAR_ETAT_B', iocc=iocc, nbval=nbchar, vect=char_meca, &
                    nbret=n1)
!
        chmome = '&&RC36SI_B'//k8b
        call rc36cm(iocc, 'B', nbma, listma, nbchar, &
                    char_meca, chmome)
        zk24(jmomeb+iocc-1) = chmome
        AS_DEALLOCATE(vi=char_meca)
!
! ------ TRANSITOIRE THERMIQUE ASSOCIE A LA SITUATION:
!        ---------------------------------------------
        call getvis(motcl1, 'NUME_RESU_THER', iocc=iocc, nbval=0, nbret=n1)
        nbth = -n1
        call jecroc(jexnum('&&RC3600.SITU_THERMIQUE', iocc))
        nbm = max(1, nbth)
        call jeecra(jexnum('&&RC3600.SITU_THERMIQUE', iocc), 'LONMAX', nbm)
!
        if (nbth .eq. 0) then
            call jeecra(jexnum('&&RC3600.SITU_THERMIQUE', iocc), 'LONUTI', 0)
        else
            call jeecra(jexnum('&&RC3600.SITU_THERMIQUE', iocc), 'LONUTI', nbth)
            call jeveuo(jexnum('&&RC3600.SITU_THERMIQUE', iocc), 'E', jreth)
            call getvis(motcl1, 'NUME_RESU_THER', iocc=iocc, nbval=nbth, vect=zi(jreth), &
                        nbret=n1)
!     ------------------------------------------------------------------
!                   RESULTATS DES CALCULS THERMIQUES
!     ------------------------------------------------------------------
            call rc36th(noma, nbma, listma, zk24(jchth), iocc, &
                        nbth, zi(jreth))
        end if
!
    end do
!
    do iocc = 1, nbseis, 1
!
        call codent(nbsitu+iocc, 'D0', k8b)
!
        call getvis(motcl2, 'NUME_GROUPE', iocc=iocc, scal=nume, nbret=n1)
        zi(jsigr+2*(nbsitu+iocc)-2) = nume
        zi(jsigr+2*(nbsitu+iocc)-1) = nume
!
        zl(jcombi+nbsitu+iocc-1) = .true.
!
! ------ LE NUMERO DE SITUATION:
!        -----------------------
        call getvis(motcl2, 'NUME_SITU', iocc=iocc, scal=zi(jnsitu+nbsitu+iocc-1), nbret=n1)
!
! ------ LE NOMBRE D'OCCURRENCE:
!        -----------------------
        call getvis(motcl2, 'NB_OCCUR', iocc=iocc, scal=nocc, nbret=n1)
        zi(jnbocc+2*(nbsitu+iocc)-2) = nocc
        call getvis(motcl2, 'NB_CYCL_SEISME', iocc=iocc, scal=nscy, nbret=n1)
        zi(jnbocc+2*(nbsitu+iocc)-1) = nscy
!
! ------ ETAT DE CHARGEMENT:
!        -------------------
        call getvis(motcl2, 'CHAR_ETAT', iocc=iocc, nbval=0, nbret=n1)
        nbchar = -n1
        call wkvect('&&RC36SI.CHAR_ETAT', 'V V I', nbchar, jchar1)
        call getvis(motcl2, 'CHAR_ETAT', iocc=iocc, nbval=nbchar, vect=zi(jchar1), &
                    nbret=n1)
!
        chmome = '&&RC36SI_A'//k8b
        call rc36cm(iocc, 'S', nbma, listma, nbchar, &
                    zi(jchar1), chmome)
        zk24(jmomea+nbsitu+iocc-1) = chmome
        call jedetr('&&RC36SI.CHAR_ETAT')
!
! ------ TRANSITOIRE THERMIQUE ASSOCIE A LA SITUATION:
!        ---------------------------------------------
        nbth = 0
        call jecroc(jexnum('&&RC3600.SITU_THERMIQUE', nbsitu+iocc))
        nbm = max(1, nbth)
        call jeecra(jexnum('&&RC3600.SITU_THERMIQUE', nbsitu+iocc), 'LONMAX', nbm)
!
        call jeecra(jexnum('&&RC3600.SITU_THERMIQUE', nbsitu+iocc), 'LONUTI', 0)
!
    end do
!
    call ordis(nume_group, nbgr)
!
    if (nbgr .gt. 3 .and. yapass) then
        call utmess('F', 'POSTRCCM_34')
    end if
!
!     ------------------------------------------------------------------
! --- ON AJOUTE 1 GROUPE POUR LES SITUATIONS DE PASSAGE
    if (yapass) nbgr = nbgr+1
!
!     ------------------------------------------------------------------
! --- DEFINITION DES GROUPES
    call wkvect('&&RC3600.SITU_NUME_GROUP', 'V V I', nbgr, jnumgr)
    call wkvect('&&RC3600.SITU_SEISME', 'V V I', nbgr, jseigr)
    call jecrec('&&RC3600.LES_GROUPES', 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbgr)
!
    if (yapass) then
        nbgrt = nbgr-1
    else
        nbgrt = nbgr
    end if
    do ig = 1, nbgrt, 1
!
        numgr = nume_group(ig)
!
        zi(jnumgr+ig-1) = numgr
!
! ------ ON COMPTE LES SITUATIONS DU GROUPE
        nbsigr = 0
        do iocc = 1, nbsitu, 1
            call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                nbvg = -n1
                call wkvect('&&RC36SI.VALE_GR', 'V V I', nbvg, jnbvg)
                call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=nbvg, vect=zi(jnbvg), &
                            nbret=n1)
                do ing = 1, nbvg
                    if (zi(jnbvg+ing-1) .eq. numgr) nbsigr = nbsigr+1
                end do
                call jedetr('&&RC36SI.VALE_GR')
            end if
            call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=2, vect=numpas, &
                            nbret=n1)
                if (numpas(1) .eq. numgr) nbsigr = nbsigr+1
                if (numpas(2) .eq. numgr) nbsigr = nbsigr+1
            end if
        end do
!
! ------ ON COMPTE LES SITUATIONS DE SEISME
        do iocc = 1, nbseis, 1
            call getvis(motcl2, 'NUME_GROUPE', iocc=iocc, scal=numgs, nbret=n1)
            if (numgs .eq. numgr) nbsigr = nbsigr+1
        end do
!
! ------ ON STOCKE LE NUMERO DE L'OCCURRENCE
        call jecroc(jexnum('&&RC3600.LES_GROUPES', numgr))
        call jeecra(jexnum('&&RC3600.LES_GROUPES', numgr), 'LONMAX', nbsigr)
        call jeveuo(jexnum('&&RC3600.LES_GROUPES', numgr), 'E', jnsg)
        ii = 0
        do iocc = 1, nbsitu, 1
            call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                nbvg = -n1
                call wkvect('&&RC36SI.VALE_GR', 'V V I', nbvg, jnbvg)
                call getvis(motcl1, 'NUME_GROUPE', iocc=iocc, nbval=nbvg, vect=zi(jnbvg), &
                            nbret=n1)
                do ing = 1, nbvg
                    if (zi(jnbvg+ing-1) .eq. numgr) then
                        ii = ii+1
                        zi(jnsg+ii-1) = iocc
                    end if
                end do
                call jedetr('&&RC36SI.VALE_GR')
            end if
            call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvis(motcl1, 'NUME_PASSAGE', iocc=iocc, nbval=2, vect=numpas, &
                            nbret=n1)
                if (numpas(1) .eq. numgr) then
                    ii = ii+1
                    zi(jnsg+ii-1) = iocc
                end if
                if (numpas(2) .eq. numgr) then
                    ii = ii+1
                    zi(jnsg+ii-1) = iocc
                end if
            end if
        end do
        do iocc = 1, nbseis, 1
            call getvis(motcl2, 'NUME_GROUPE', iocc=iocc, scal=numgs, nbret=n1)
            if (numgs .eq. numgr) then
                ii = ii+1
                zi(jnsg+ii-1) = nbsitu+iocc
                if (zi(jseigr+ig-1) .ne. 0) then
                    vali(1) = numgr
                    vali(2) = iocc
                    vali(3) = zi(jseigr+ig-1)
                    call utmess('F', 'POSTRCCM_26', ni=3, vali=vali)
                end if
                zi(jseigr+ig-1) = iocc
            end if
        end do
!
    end do
!     ------------------------------------------------------------------
! --- TRAITEMENT DES SITUATIONS DE PASSAGE
    if (yapass) then
!
        call wkvect('&&RC32SI.PASSAGE_SIT', 'V V I', 3, jspas)
!
        nbsg1 = 0
        nbsg2 = 0
        nbsg3 = 0
        nbp12 = 0
        nbp23 = 0
        nbp13 = 0
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 1 .and. numg2 .eq. 1) then
                nbsg1 = nbsg1+1
            else if (numg1 .eq. 1 .and. numg2 .eq. 2) then
                nbsg1 = nbsg1+1
                nbp12 = nbp12+1
                zi(jsp12+nbp12-1) = iocc
            else if (numg1 .eq. 2 .and. numg2 .eq. 2) then
                nbsg2 = nbsg2+1
            else if (numg1 .eq. 2 .and. numg2 .eq. 3) then
                nbsg2 = nbsg2+1
                nbp23 = nbp23+1
                zi(jsp23+nbp23-1) = iocc
            else if (numg1 .eq. 3 .and. numg2 .eq. 3) then
                nbsg3 = nbsg3+1
            else if (numg1 .eq. 1 .and. numg2 .eq. 3) then
                nbsg3 = nbsg3+1
                nbp13 = nbp13+1
                zi(jsp13+nbp13-1) = iocc
            end if
        end do
        call jeecra('&&RC32SI.PASSAGE_1_2', 'LONUTI', nbp12)
        call jeecra('&&RC32SI.PASSAGE_2_3', 'LONUTI', nbp23)
        call jeecra('&&RC32SI.PASSAGE_1_3', 'LONUTI', nbp13)
        zi(jspas) = nbsg1
        zi(jspas+1) = nbsg2
        zi(jspas+2) = nbsg3
!
        zi(jnumgr+nbgr-1) = -nbgr
        call jecroc(jexnum('&&RC3600.LES_GROUPES', nbgr))
        call jeecra(jexnum('&&RC3600.LES_GROUPES', nbgr), 'LONMAX', ndim)
        call jeveuo(jexnum('&&RC3600.LES_GROUPES', nbgr), 'E', jnsg)
!
        ii = 0
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 1 .and. numg2 .eq. 1) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 1 .and. numg2 .eq. 2) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 2 .and. numg2 .eq. 2) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 2 .and. numg2 .eq. 3) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 3 .and. numg2 .eq. 3) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        do iocc = 1, ndim, 1
            numg1 = zi(jsigr+2*iocc-2)
            numg2 = zi(jsigr+2*iocc-1)
            if (numg1 .eq. 1 .and. numg2 .eq. 3) then
                ii = ii+1
                zi(jnsg+ii-1) = iocc
            end if
        end do
        call jeecra(jexnum('&&RC3600.LES_GROUPES', nbgr), 'LONUTI', ii)
    end if
!
    AS_DEALLOCATE(vi=nume_group)
!
end subroutine
