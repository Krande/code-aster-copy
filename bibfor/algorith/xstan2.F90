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

subroutine xstan2(noma, modele, crit2, lfiss)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/cncinv.h"
#include "asterfort/cnscno.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref2.h"
#include "asterfort/infdbg.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "asterfort/xcrvol.h"
#include "asterfort/isnomi.h"
!
    character(len=8) :: modele, noma
    real(kind=8) :: crit2(2)
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (PREPARATION)
!
! REATION LISTE DE NOEUDS OU IL FAUDRA ANNULER LES DDLS HEAVISIDE
! POUR DES RAISONS DE CONDITIONNEMENT DE MATRICE ET POUR EVITER
! DES PIVOTS NULS DANS LA MATRICE DE RAIDEUR
! (VOIR BOOK V 15/04/05)
!
! ----------------------------------------------------------------------
!
!
!  IN  CRIT2 : CRITERE (RAPPORT MAXIMUM ENTRE LES VOLUMES)
!  IN  NOMA   : NOM DE L'OBJET MAILLAGE
!  I/O MODELE   : NOM DE LA SD MODELE_XFEM
!
!
!
!
    character(len=24) :: geom
    character(len=19) :: ces(8), cel(8), cnxinv, noxfem, cns2, ligrel
    character(len=16) :: notype
    character(len=8) :: typma, lirefe(10), elrefp
    character(len=2) :: ch2
    real(kind=8) :: crit, vhea, vtot, crimaxi
    integer(kind=8) :: jcesd(8), jcesl(8), jcesv(8), iad
    integer(kind=8) :: jnoxfl, itypma, nncp, ibid, ier
    integer(kind=8) :: ifm, niv, jpint, jcnse
    integer(kind=8) :: jconx2, adrma, ndime, ndim, nbno
    integer(kind=8) :: nbmano, nbnoma, nuno, ino, ino2, numa, numa2, ima
    integer(kind=8) :: itypel, nbelr, igeom, nuno2, inoloc, cpt
    integer(kind=8) :: i, j, nheav, iheav, nfiss, nfiss2, ifiss, nse, nnose
    aster_logical :: lelim, lfiss
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: nbsp(:) => null()
    integer(kind=8), pointer :: maille(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: nbsp2(:) => null()
    integer(kind=8), pointer :: tmdim(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
    write (ifm, *) 'RECHERCHE DES DDLS HEAVISIDE A ANNULER'
!
!     RECUPERATION DES DONNEES SUR LE MAILLAGE
!
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=tmdim)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
    call jeveuo(modele//'.MAILLE', 'L', vi=maille)
!     CONNECTIVITE INVERSEE
    cnxinv = '&&XSTAN2.CNCINV'
    call cncinv(noma, [ibid], 0, 'V', cnxinv)
    noxfem = modele//'.NOXFEM'
    ligrel = modele//'.MODELE'
    cns2 = '&&XCONNO.CNS2'
    geom = '&&XSTAN2.GEOM'
!      CALL CNOCNS(MODELE//'.NOXFEM','V',NOXFEM)
    call jeveuo(cns2//'.CNSL', 'E', jnoxfl)
    call jeveuo(cns2//'.CNSV', 'E', vi=cnsv)
!
!     RECUPERATION DES DONNEES ELEMENTAIRES XFEM
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
    call jeveuo('&&XTYELE.NBSP2', 'L', vi=nbsp2)
    cel(1) = modele//'.TOPOSE.LON'
    cel(2) = modele//'.STNO'
    cel(3) = modele//'.TOPOSE.HEA'
    cel(4) = modele//'.TOPOSE.CNS'
    cel(5) = modele//'.TOPOSE.PIN'
    cel(6) = modele//'.TOPOSE.PMI'
    cel(7) = modele//'.FISSNO'
    cel(8) = modele//'.TOPONO.HNO'
    do i = 1, 8
        call codent(i, 'G', ch2)
        ces(i) = '&&XSTAN2.CES'//ch2
        call jeexin(cel(i)//'.CELD', ier)
        if (ier .eq. 0) goto 10
        call celces(cel(i), 'V', ces(i))
        call jeveuo(ces(i)//'.CESD', 'L', jcesd(i))
        call jeveuo(ces(i)//'.CESL', 'L', jcesl(i))
        call jeveuo(ces(i)//'.CESV', 'E', jcesv(i))
10      continue
    end do
!
    cpt = 0
!     BOUCLE SUR LES NOEUDS DU MAILLAGE
    do nuno = 1, nbno
        if (.not. zl(jnoxfl-1+2*nuno)) goto 20
! --- RECUP DU NUMERO LOCAL INO DU NOEUD NUNO DANS LA MAILLE X-FEM NUMA
        numa = cnsv(2*(nuno-1)+1)
        ino = cnsv(2*(nuno-1)+2)
        nfiss = nbsp(numa)
        ASSERT(nfiss .ge. 1)
        nheav = max(1, nbsp2(numa))
!
!       RECUPERATION DES MAILLES CONTENANT LE NOEUD
        call jelira(jexnum(cnxinv, nuno), 'LONMAX', nbmano)
        call jeveuo(jexnum(cnxinv, nuno), 'L', adrma)
!
!       BOUCLE SUR LES DDL HEAVISIDE
        do iheav = 1, nheav
            if (nfiss .eq. 1) then
                ifiss = 1
            else
                call cesexi('S', jcesd(7), jcesl(7), numa, ino, &
                            iheav, 1, iad)
                ifiss = zi(jcesv(7)-1+iad)
            end if
            call cesexi('S', jcesd(2), jcesl(2), numa, ino, &
                        ifiss, 1, iad)
            if (zi(jcesv(2)-1+iad) .ne. 1 .and. zi(jcesv(2)-1+iad) .ne. 3) goto 40
!
!         INITIALISATION DES VOLUMES
            vhea = 0
            vtot = 0
!         BOUCLE SUR LES MAILLES SUPPORT DU NOEUD
            do ima = 1, nbmano
                numa2 = zi(adrma-1+ima)
                ndime = tmdim(typmail(numa2))
                if (ndime .lt. ndim) goto 50
                itypma = typmail(numa2)
                call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
                if (.not. ismali(typma)) then
                    if (ndim .eq. 2) then
                        nnose = 6
                    else
                        nnose = 10
                    end if
                else
                    nnose = ndim+1
                end if
!
!         1ER ELEMENT DE REFERENCE ASSOCIE A LA MAILLE
                itypel = maille(numa2)
                call jenuno(jexnum('&CATA.TE.NOMTE', itypel), notype)
                call elref2(notype, 10, lirefe, nbelr)
                elrefp = lirefe(1)
!           NOMBRE DE NOEUDS
                nbnoma = zi(jconx2+numa2)-zi(jconx2+numa2-1)
!           RECUP DU NUM DE FISS CORRESPONDANT À IHEAV DE NUNO DS NUMA2
!           ET DU NUMERO LOCALE INOLOC DANS LA MAILLE
                do ino2 = 1, nbnoma
                    if (connex(zi(jconx2+numa2-1)+ino2-1) .eq. nuno) then
                        inoloc = ino2
                        goto 200
                    end if
                end do
                ASSERT(.false.)
200             continue
                nfiss2 = nbsp(numa2)
                if (nfiss2 .eq. 1) then
                    ifiss = 1
                else
                    call cesexi('S', jcesd(7), jcesl(7), numa2, inoloc, &
                                iheav, 1, iad)
                    ifiss = zi(jcesv(7)-1+iad)
                end if
!           CREATION DE VECTEUR DES COORDONNÉES DE LA MAILLE IMA
!           AVEC DES VALEURS CONTIGUES
                call wkvect(geom, 'V V R', ndim*nbnoma, igeom)
                do ino2 = 1, nbnoma
                    nuno2 = connex(zi(jconx2+numa2-1)+ino2-1)
                    do j = 1, ndim
                        zr(igeom-1+ndim*(ino2-1)+j) = vale(3*( &
                                                           nuno2-1)+j)
                    end do
                end do
!
!           RECUPERATION DU NOMBRE TOTAL DE SOUS ELEMENTS
                call cesexi('S', jcesd(1), jcesl(1), numa2, 1, &
                            1, 1, iad)
                nse = zi(jcesv(1)-1+iad)
!           POINTEUR DE CONNECTIVITÉ DU SOUS ELEMENT
                call cesexi('S', jcesd(4), jcesl(4), numa2, 1, &
                            1, 1, iad)
                jcnse = jcesv(4)-1+iad
!           POINTEUR DE COORDONNÉES DES POINTS D'INTERSECTIONS
                call cesexi('S', jcesd(5), jcesl(5), numa2, 1, &
                            1, 1, iad)
                jpint = jcesv(5)-1+iad
                call xcrvol(nse, ndim, jcnse, nnose, jpint, &
                            igeom, elrefp, inoloc, nbnoma, jcesd(3), &
                            jcesl(3), jcesv(3), numa2, iheav, nfiss2, vhea, &
                            jcesd(8), jcesl(8), jcesv(8), lfiss, vtot)
                call jedetr(geom)
50              continue
            end do
!         CALCUL DU CRITERE
            if (.not. isnomi(elrefp, inoloc)) then
                crimaxi = crit2(1)**ndim
            else
                crimaxi = crit2(2)**ndim
            end if
            crit = vhea/vtot
            if (crit .lt. crimaxi) then
                cpt = cpt+1
!           BOUCLE SUR LES MAILLES SUPPORT DU NOEUD
                do ima = 1, nbmano
                    numa2 = zi(adrma-1+ima)
!             MISE À ZÉRO DU STATUT DANS TOUS LES ÉLÉMENTS DU SUPPORT
                    nbnoma = zi(jconx2+numa2)-zi(jconx2+numa2-1)
                    do ino2 = 1, nbnoma
                        if (connex(zi(jconx2+numa2-1)+ino2-1) .eq. nuno) then
                            if (nbsp(numa2) .eq. 1) then
                                ifiss = 1
                            else if (nbsp(numa2) .eq. 0) then
                                goto 150
                            else
                                call cesexi('S', jcesd(7), jcesl(7), numa2, ino2, &
                                            iheav, 1, iad)
                                ifiss = zi(jcesv(7)-1+iad)
                            end if
                            call cesexi('S', jcesd(2), jcesl(2), numa2, ino2, &
                                        ifiss, 1, iad)
                            if (zi(jcesv(2)-1+iad) .eq. 3) then
                                zi(jcesv(2)-1+iad) = 2
                            elseif (zi(jcesv(2)-1+iad) .eq. 1) then
                                zi(jcesv(2)-1+iad) = 0
                            end if
                            goto 150
                        end if
                    end do
                    ASSERT(.false.)
150                 continue
                end do
            end if
!
40          continue
        end do
!       ELIMINATION DE LA LISTE DES NOEUDS XFEM SI NECESSAIRE
        lelim = .true.
        do iheav = 1, nheav
            nfiss = nbsp(numa)
            if (nfiss .eq. 1) then
                ifiss = 1
            else
                call cesexi('S', jcesd(7), jcesl(7), numa, ino, &
                            iheav, 1, iad)
                ifiss = zi(jcesv(7)-1+iad)
            end if
            call cesexi('S', jcesd(2), jcesl(2), numa, ino, &
                        ifiss, 1, iad)
            if (abs(zi(jcesv(2)-1+iad)) .ge. 1) lelim = .false.
        end do
        if (lelim) then
            zl(jnoxfl-1+2*(nuno-1)+1) = .false.
            zl(jnoxfl-1+2*(nuno-1)+2) = .false.
        end if
20      continue
    end do
!
    write (ifm, *) 'NOMBRE DE NOEUDS OU LES DDLS H SONT MIS A ZERO :', cpt
!
! --- CONVERSION CHAM_NO_S -> CHAM_NO POUR MODELE.NOXFEM
!
    call cnscno(cns2, noxfem(1:13)//".NUMEQ", 'NON', 'G', noxfem, &
                'F', ibid)
!
! --- CONVERSION CHAM_ELEM_S -> CHAM_ELEM POUR MODELE.STNO
!
    call cescel(ces(2), ligrel, 'INI_XFEM_ELNO', 'PSTANO', 'OUI', &
                nncp, 'G', cel(2), 'F', ibid)
!
! --- MENAGE
!
    call jedetr(cnxinv)
    call detrsd('CHAM_NO_S', cns2)
    do i = 1, 8
        call jeexin(ces(i)//'.CESD', ier)
        if (ier .eq. 0) goto 130
        call detrsd('CHAM_ELEM_S', ces(i))
130     continue
    end do
    call jedema()
end subroutine
