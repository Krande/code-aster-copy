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

subroutine exphgl(nomres, typsd, modcyc, profno, indirf, &
                  mailsk, nbsec, numdia, nbmode)
    implicit none
!
!  BUT:
!
!  RESTITUER LES RESULTATS ISSUS D'UN CALCUL CYCLIQUE
!     => RESULTAT COMPOSE DEJA ALLOUE PAR LA
!        ROUTINE APPELLANTE
!
!  DONNEES DU NUMEEQUA DEJA CONSTITUE ET DE LA TABLE INDIRECTION
!  DES NUMEROS EQUATIONS CORRESPONDANTES (COLLECTION NUMEROTEE
!  POINTEE PAR LES NUMEROS DE SECTEUR)
!-----------------------------------------------------------------------
!
! NOMRES  /I/: NOM UT DU CONCEPT RESULTAT A REMPLIR
! MODCYC  /I/: NOM UT DU RESULTAT ISSU DU CALCUL CYCLIQUE
! PROFNO  /I/: NOM K19 DU PROFIL CHAMNO DEJA CONSTITUE
! INDIRF  /I/: NOM K24 DE LA FAMILLE DES INDIRECTIONS
! MAILSK  /I/: NOM K8 DU MAILLAGE SKELETTE
! TYPSD   /I/: NOM DU TYPE DE STRUCTURE DE DONNEES RESULTAT
! NBSEC   /I/: NBRE DE SECTEUR
! NUMDIA  /I/: NUMERO DU DIAMETRE
!
!
!
!
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rotchm.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslipa.h"
#include "asterfort/rsnoch.h"
#include "asterfort/vtcrea.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: nomres, modcyc, mailsk, k8b, modcys
    character(len=16) :: depl, typsd
    character(len=19) :: chamva, profno, chamno
    character(len=24) :: indirf, crefe(2), nomchc, pfchno, nomchs
    real(kind=8) :: depi, genek, beta
    integer(kind=8) :: nbmode, ibid, iret, neqsec, llfreq, ltveco, ldfreq, ldkge
    integer(kind=8) :: ldmge, ldom2, ldomo, nbnot, nbcmp, nbsec, neq, ires2
    integer(kind=8) :: numdia
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icomp, ieqf, ieqi, ier, j, k
    integer(kind=8) :: ldtyd, ltinds, n1, nddcou
    real(kind=8), pointer :: teta_secteur(:) => null()
    integer(kind=8), pointer :: skeleton(:) => null()
    real(kind=8), pointer :: nllcham(:) => null()
    real(kind=8), pointer :: nltvesi(:) => null()
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
    depi = r8depi()
!
!-----REMPLISSAGE DU CREFE POUR CREATION CHAMNO-------------------------
!
    crefe(1) = mailsk
    crefe(2) = profno
!
!-----RECUPERATION DU NOMBRE DE DDL PHYSIQUES DU SECTEUR----------------
!
    call rsexch('F', modcyc, 'DEPL', 1, chamno, &
                ier)
    call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=pfchno)
!
    call dismoi('NB_EQUA', pfchno, 'NUME_EQUA', repi=neqsec)
!     -- QUESTION "POURRIE" :
    call dismoi('NOM_GD ', pfchno, 'NUME_EQUA', repi=ibid, repk=k8b)
    call dismoi('NB_CMP_MAX', k8b, 'GRANDEUR', repi=nbcmp)
!
!-----RECUPERATION DU NOMBRE DE DDL PHYSIQUES GLOBAUX-------------------
!
    call jelira(profno//'.DEEQ', 'LONMAX', neq)
    neq = neq/2
!
!-----RECUPERATION DES FREQUENCES---------------------------------------
!
    if ((typsd(1:9) .eq. 'MODE_MECA') .or. (typsd(1:4) .eq. 'BASE')) then
        call rslipa(modcyc, 'FREQ', '&&EXPHGL.LIR8', llfreq, n1)
    else
        call rslipa(modcyc, 'INST', '&&EXPHGL.LIR8', llfreq, n1)
    end if
!
!-----ALLOCATION DES VECTEURS DE TRAVAIL--------------------------------
!
    call wkvect('&&EXPHGL.VEC.REEL', 'V V R', neqsec, ltveco)
!
!-----CALCUL DU TETA DE CHAQUE SECTEUR----------------------------------
!
    AS_ALLOCATE(vr=teta_secteur, size=nbsec)
    do i = 1, nbsec
        teta_secteur(i) = depi*(i-1)/nbsec
    end do
!
!-----RECUPERATION DE L'INDIRECTION SQUELETTE---------------------------
!
    call jeveuo(mailsk//'.INV.SKELETON', 'L', vi=skeleton)
    call dismoi('NB_NO_MAILLA', mailsk, 'MAILLAGE', repi=nbnot)
!
!***********************************************************************
!
    call getvid('CYCLIQUE', 'RESULTAT2', iocc=1, scal=modcys, nbret=ires2)
!
    icomp = 0
!
!  CALCUL DU DEPHASAGE INTER-SECTEUR
!
    beta = numdia*(depi/nbsec)
!
!  BOUCLE SUR LES MODES PROPRES DU DIAMETRE COURANT
!
    do i = 1, nbmode
        icomp = icomp+1
        call rsexch('F', modcyc, 'DEPL', i, nomchc, &
                    iret)
        call jeveuo(nomchc(1:19)//'.VALE', 'L', ltveco)
        if (ires2 .ne. 0) then
            call rsexch('F', modcys, 'DEPL', i, nomchs, &
                        iret)
            call jeveuo(nomchs(1:19)//'.VALE', 'L', vr=nltvesi)
        end if
!
!
!***********************************************************************
!
        call rsexch(' ', nomres, depl, i, chamva, &
                    iret)
        call vtcrea(chamva, crefe, 'G', 'R', neq)
        call rsnoch(nomres, depl, i)
        call jeveuo(chamva//'.VALE', 'E', vr=nllcham)
!
!  COMMUN POUR MODE_MECA ET BASE_MODALE
!
        if ((typsd(1:9) .eq. 'MODE_MECA')) then
            call rsadpa(nomres, 'E', 1, 'FREQ', i, &
                        0, sjv=ldfreq, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'RIGI_GENE', i, &
                        0, sjv=ldkge, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'MASS_GENE', i, &
                        0, sjv=ldmge, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'OMEGA2', i, &
                        0, sjv=ldom2, styp=k8b)
            call rsadpa(nomres, 'E', 1, 'NUME_MODE', i, &
                        0, sjv=ldomo, styp=k8b)
            genek = (zr(llfreq+icomp-1)*depi)**2
            zr(ldfreq) = zr(llfreq+icomp-1)
            zr(ldkge) = genek
            zr(ldmge) = 1.d0
            zr(ldom2) = genek
            zi(ldomo) = i
!
!  SPECIFIQUE A BASE_MODALE
!
            call rsadpa(nomres, 'E', 1, 'TYPE_DEFO', i, &
                        0, sjv=ldtyd, styp=k8b)
            zk16(ldtyd) = 'PROPRE          '
        else
            call rsadpa(nomres, 'E', 1, 'INST', i, &
                        0, sjv=ldfreq, styp=k8b)
            zr(ldfreq) = zr(llfreq+icomp-1)
        end if
!
!  BOUCLE SUR LES SECTEURS
!
        do k = 1, nbsec
            call jeveuo(jexnum(indirf, k), 'L', ltinds)
            call jelira(jexnum(indirf, k), 'LONMAX', nddcou)
            nddcou = nddcou/2
            do j = 1, nddcou
                ieqi = zi(ltinds+(j-1)*2)
                ieqf = zi(ltinds+(j-1)*2+1)
                if (ires2 .ne. 0) then
                    nllcham(ieqf) = sin( &
                                        (k-1)*beta)*zr(ltveco+ieqi-1)+cos((k-1)*beta)*nltvesi(1&
                                        &+ieqi-1 &
                                        )
                else
                    nllcham(ieqf) = zr(ltveco+ieqi-1)
                end if
            end do
        end do
!
!  PRISE EN COMPTE ROTATION SUR CHAQUE SECTEUR
!
        call rotchm(profno, nllcham, teta_secteur, nbsec, skeleton, &
                    nbnot, nbcmp, 3)
!
        call jelibe(nomchc(1:19)//'.VALE')
        if (ires2 .ne. 0) then
            call jelibe(nomchs(1:19)//'.VALE')
        end if
    end do
!
    call jedetr('&&EXPHGL.VEC.REEL')
    call jedetr('&&EXPHGL.ORDRE.FREQ')
    AS_DEALLOCATE(vr=teta_secteur)
    call jedetr('&&EXPHGL.TETGD')
    call jedetr('&&EXPHGL.LIR8')
!
    call jedema()
end subroutine
