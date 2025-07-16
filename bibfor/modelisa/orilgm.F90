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

subroutine orilgm(noma)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/orilma.h"
#include "asterfort/ornorm.h"
#include "asterfort/orvlse.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: noma
! ======================================================================
!
!     ORILGM  --  LE BUT EST DE REORIENTER, SI C'EST NECESSAIRE,
!                 LES MAILLES DE PEAU DE GROUPES DE MAILLES
!                 DONNES SOUS LES MOTS CLES :
!                 'ORIE_PEAU_2D' EN 2D (ou 3D)
!                 ET 'ORIE_NORM_COQUE' EN 3D
!                 DE TELLE FACON A CE QUE LA NORMALE A LA MAILLE DE
!                 PEAU SOIT EXTERIEURE AU VOLUME.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    MODELZ         IN    K*      NOM DU MODELE
!
! ========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: ifm, niv, nbf1, jjj, n1, n2
    integer(kind=8) :: n3, noeud, iocc, ier, ndim, igr, ng, nbmail, norit, norien
    integer(kind=8) :: ntrait, jjv, nbmato, ima, nbmavi, jmavi, k
    integer(kind=8) :: ncf3, ngs, jgs, nbmasu, jmafr, nconex
    real(kind=8) :: vect(3)
    aster_logical :: reorie, orivec
    character(len=8) :: k8b
    character(len=16) :: mofac, mofc3d
    character(len=24) :: grmama, nnoeud, gmat
    character(len=24) :: valk(2)
    integer(kind=8), pointer :: listCellNume(:) => null()
!
! ========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
    call infniv(ifm, niv)
!
! --- INITIALISATIONS :
!     ---------------
!
    norit = 0
    reorie = .true.
    vect(:) = 0.d0
!
    mofac = 'ORIE_PEAU'
    mofc3d = 'ORIE_LIGNE'
!
    call getfac(mofac, nbf1)
    call getfac(mofc3d, ncf3)
!
! --- RECUPERATION DU MAILLAGE ASSOCIE AU MODELE :
!     ------------------------------------------
    grmama = noma//'.GROUPEMA'
!
! --- RECUPERATION DE LA DIMENSION (2 OU 3) DU PROBLEME :
!     -------------------------------------------------
    call dismoi('Z_CST', noma, 'MAILLAGE', repk=k8b)
    if (k8b(1:3) .eq. 'OUI') then
        ndim = 2
    else
        ndim = 3
    end if

! --- TRAITEMENT DE 'ORIE_PEAU' :
!     ----------------------------
!
    do iocc = 1, nbf1
        call getvem(noma, 'GROUP_MA', mofac, 'GROUP_MA_PEAU', iocc, &
                    0, k8b, ng)
        ng = -ng
        call wkvect('&&ORILGM.WORK', 'V V K24', ng, jjj)
        call getvem(noma, 'GROUP_MA', mofac, 'GROUP_MA_PEAU', iocc, &
                    ng, zk24(jjj), ng)
!        PRESENCE DE GROUP_MA_INTERNE ?
!        ---------------------------
        call getvtx(mofac, 'GROUP_MA_INTERNE', iocc=iocc, nbval=0, nbret=ngs)
        if (ngs .ne. 0) then
            ngs = -ngs
            call wkvect('&&ORILGM.WORK2', 'V V K24', ngs, jgs)
            call getvem(noma, 'GROUP_MA', mofac, 'GROUP_MA_INTERNE', iocc, &
                        ngs, zk24(jgs), ngs)
            call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmato)
            call wkvect('&&ORILGM.WORK3', 'V V I', nbmato, jjv)
            do ima = 1, nbmato
                zi(jjv+ima-1) = 0
            end do
            do igr = 1, ngs
                gmat = zk24(jgs+igr-1)
                call jelira(jexnom(grmama, gmat), 'LONMAX', nbmavi)
                call jeveuo(jexnom(grmama, gmat), 'L', jmavi)
                do ima = 1, nbmavi
                    zi(jjv+zi(jmavi+ima-1)-1) = 1
                end do
            end do
!          NOMBRE DE MAILLES 'INTERNES' (SANS DOUBLON) : NBMASU
            nbmasu = 0
            do ima = 1, nbmato
                nbmasu = nbmasu+zi(jjv+ima-1)
            end do
!          LISTE DES MAILLES 'INTERNES' (SANS DOUBLON) : ZI(jmafr)
            if (nbmasu .ne. 0) then
                call wkvect('&&ORILGM.GROUP_MA_FRONT', 'V V I', nbmasu, jmafr)
                k = 0
                do ima = 1, nbmato
                    if (zi(jjv+ima-1) .eq. 1) then
                        k = k+1
                        zi(jmafr+k-1) = ima
                    end if
                end do
            end if
            call jedetr('&&ORILGM.WORK3')
            call jedetr('&&ORILGM.WORK2')
        else
            nbmasu = 0
            call wkvect('&&ORILGM.GROUP_MA_FRONT', 'V V I', 1, jmafr)
        end if
!
        do igr = 1, ng
            gmat = zk24(jjj+igr-1)
            call jelira(jexnom(grmama, gmat), 'LONUTI', nbmail)
            call jeveuo(jexnom(grmama, gmat), 'L', vi=listCellNume)
            write (ifm, 1000) gmat, nbmail
            norien = 0
            call orilma(noma, ndim, listCellNume, nbmail, norien, &
                        ntrait, reorie, nbmasu, zi(jmafr))
            norit = norit+norien
            write (ifm, 1100) norien
            if (ntrait .ne. 0) write (ifm, 1110) ntrait
        end do
        call jedetr('&&ORILGM.WORK')
        call jedetr('&&ORILGM.GROUP_MA_FRONT')
    end do

!
!
! --- TRAITEMENT DE 'ORIE_LIGNE':
!     ------------------------------
!
    do iocc = 1, ncf3
        orivec = .false.
        call getvr8(mofc3d, 'VECT_TANG', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            orivec = .true.
            call getvr8(mofc3d, 'VECT_TANG', iocc=iocc, nbval=-n1, vect=vect, &
                        nbret=n1)
            call getvtx(mofc3d, 'NOEUD', iocc=iocc, nbval=0, nbret=n2)
            if (n2 .ne. 0) then
                call getvtx(mofc3d, 'NOEUD', iocc=iocc, scal=nnoeud, nbret=n2)
                noeud = char8_to_int(nnoeud)

                if (noeud .eq. 0) then
                    call utmess('F', 'MODELISA5_97', sk=nnoeud)
                end if
            else
                call getvtx(mofc3d, 'GROUP_NO', iocc=iocc, scal=nnoeud, nbret=n3)
                call utnono(' ', noma, 'NOEUD', nnoeud, k8b, &
                            ier)
                if (ier .eq. 10) then
                    call utmess('F', 'MODELISA8_75', sk=nnoeud)
                else if (ier .eq. 1) then
                    valk(1) = nnoeud
                    valk(2) = k8b
                    call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
                end if
                noeud = char8_to_int(k8b)
            end if
        end if
        call getvem(noma, 'GROUP_MA', mofc3d, 'GROUP_MA', iocc, &
                    0, k8b, ng)
        ng = -ng
        call wkvect('&&ORILGM.WORK', 'V V K24', ng, jjj)
        call getvem(noma, 'GROUP_MA', mofc3d, 'GROUP_MA', iocc, &
                    ng, zk24(jjj), ng)
        if (orivec) then
            do igr = 1, ng
                gmat = zk24(jjj+igr-1)
                call jelira(jexnom(grmama, gmat), 'LONUTI', nbmail)
                call jeveuo(jexnom(grmama, gmat), 'L', vi=listCellNume)
                write (ifm, 1000) gmat, nbmail
                norien = 0
                call orvlse(noma, listCellNume, nbmail, norien, vect, &
                            noeud)
                norit = norit+norien
                write (ifm, 1100) norien
            end do
        else
            do igr = 1, ng
                gmat = zk24(jjj+igr-1)
                call jelira(jexnom(grmama, gmat), 'LONUTI', nbmail)
                call jeveuo(jexnom(grmama, gmat), 'L', vi=listCellNume)
                write (ifm, 1000) gmat, nbmail
                norien = 0
                call ornorm(noma, listCellNume, nbmail, reorie, norien, nconex, &
                            onlySkin1D_=ASTER_TRUE)
                if (nconex .gt. 1) then
                    call utmess('F', 'MESH3_6')
                end if
                norit = norit+norien
                write (ifm, 1100) norien
            end do
        end if
        call jedetr('&&ORILGM.WORK')
    end do
!
    if (norit .ne. 0) write (ifm, 1010) norit
!
1000 format('TRAITEMENT DU GROUP_MA: ', a24, ' DE ', i7, ' MAILLES')
1100 format(24x, i7, ' MAILLE(S) ONT ETE ORIENTEE(S)')
1110 format(24x, i7, ' MAILLE(S) N''ONT PAS ETE TRAITEE(S) ')
1010 format('AU TOTAL ', i7, ' MAILLE(S) ORIENTEE(S) ')
!
    call jedema()
end subroutine
