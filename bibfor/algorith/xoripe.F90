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
subroutine xoripe(modele)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/panbno.h"
#include "asterfort/provec.h"
#include "asterfort/utmasu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xelfis_lists.h"
#include "blas/ddot.h"
!
    character(len=8) :: modele
!
!
! person_in_charge: samuel.geniaut at edf.fr
!
!        ORIENTER LES SOUS-ELEMENTS DE PEAU DES ELEMENTS X-FEM
!      (ET CALCUL DE HEAV SUR LES BORDS COINCIDANT AVEC INTERACE)
!
!  IN/OUT         MODELE    : NOM DE L'OBJET MODELE
!
!     ------------------------------------------------------------------
!
    real(kind=8) :: gbo(3), gpr(3), next(3), norme, lsn
    real(kind=8) :: co(3, 3), ab(3), ac(3), n2d(3), a(3), b(3), c(3)
    integer(kind=8) :: ima, nbma, j, kk, i, ifis, nfis, iad2
    integer(kind=8) :: nbmail, iret, jm3d, ibid, jvecno
    integer(kind=8) :: numapr, numab, nbnopr, nbnobo, nbnose, nbnott(3)
    integer(kind=8) :: nbnos, id4, id6
    integer(kind=8) :: jconx2, ino, nuno
    integer(kind=8) :: ich, jcesd(5), jcesv(5), jcesl(5), iad, nse, ise, in
    integer(kind=8) :: ndime, icmp, ndim, id(3), intemp, nseori, ifm, niv, nncp
    integer(kind=8) :: s1, s2, jgrp, nmaenr, jlsnd, jlsnl
    integer(kind=8) :: nsignp, nsignm, nsignz, ihe, he, itypma
    integer(kind=8) :: ifiss, nfiss, mailvo(1)
    character(len=8) :: noma, typbo, fiss
    character(len=6) :: nompro
    parameter(nompro='XORIPE')
    character(len=2) :: kdim
    character(len=19) :: ligrel, chs(5), chlsn
    character(len=24) :: grmape, nomob, vecnor, grp(3)
    character(len=24) :: elfis_heav, elfis_ctip, elfis_hect
    character(len=19) :: pintto, cnseto, loncha, heav
    aster_logical :: quadratique
    integer(kind=8) :: itypbo
    integer(kind=8), pointer :: vnfis(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: tmdim(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    character(len=8), pointer :: vfiss(:) => null()
    integer(kind=8), pointer :: listCellNume(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
!     INITIALISATION DU NOMBRE DE SOUS-ELEMENTS RE-ORIENTES
    nseori = 0
!
    ligrel = modele//'.MODELE'
!
!     1.RECUPERATION D'INFORMATIONS DANS MODELE
!
    call jeveuo(modele//'.NFIS', 'L', vi=vnfis)
    nfis = vnfis(1)
    call jeveuo(modele//'.FISS', 'L', vk8=vfiss)
!
!     RECUPERATION DU MAILLAGE ASSOCIE AU MODELE :
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=tmdim)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
    chlsn = '&&XORIPE.CHLSN'
    call celces(modele//'.LNNO', 'V', chlsn)
    call jeveuo(chlsn//'.CESL', 'L', jlsnl)
    call jeveuo(chlsn//'.CESD', 'L', jlsnd)
    call jeveuo(chlsn//'.CESV', 'L', vr=cesv)
!
!     ------------------------------------------------------------------
!     I°) CREATION DE LA LISTE DES NUMEROS DES MAILLES DE PEAU ENRICHIES
!     ------------------------------------------------------------------
!
    grmape = '&&XORIPE.GRMAPE'
    call wkvect(grmape, 'V V I', nbma, vi=listCellNume)
!
!     INITIALISATION DU NOMBRE DE MAILLES DE LA LISTE
    nbmail = 0
!
    do ifis = 1, nfis
!
        fiss = vfiss(ifis)
!
!       REMPLISSAGE DE LA LISTE
        elfis_heav = '&&'//nompro//'.ELEMFISS.HEAV'
        elfis_ctip = '&&'//nompro//'.ELEMFISS.CTIP'
        elfis_hect = '&&'//nompro//'.ELEMFISS.HECT'
        call xelfis_lists(fiss, modele, elfis_heav, elfis_ctip, elfis_hect)
        grp(1) = elfis_heav
        grp(2) = elfis_ctip
        grp(3) = elfis_hect
!
!       BOUCLE SUR LES 3 GROUPES : HEAV, CTIP ET HECT
        do kk = 1, 3
!
            call jeexin(grp(kk), iret)
            if (iret .ne. 0) then
                call jeveuo(grp(kk), 'L', jgrp)
                call jelira(grp(kk), 'LONMAX', nmaenr)
!
!           BOUCLE SUR LES MAILLES DE CHAQUE GROUPE
                do i = 1, nmaenr
                    ima = zi(jgrp-1+i)
!             NDIME : DIMENSION TOPOLOGIQUE DE LA MAILLE
                    ndime = tmdim(typmail(ima))
                    if (ndim .eq. ndime+1) then
                        nbmail = nbmail+1
                        listCellNume(nbmail) = ima
                    end if
                end do
!               menage
                call jedetr(grp(kk))
!
            end if
        end do
!
    end do
    if (nbmail .eq. 0) goto 999
!
!     ------------------------------------------------------------------
!     II°) RECHERCHE DES MAILLES SUPPORT
!     ------------------------------------------------------------------
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    if (ndim .eq. 2) kdim = '2D'
    if (ndim .eq. 3) kdim = '3D'
!
    nomob = '&&XORIPE.NU_MAILLE_3D'
!
    call utmasu(noma, kdim, nbmail, listCellNume, nomob, &
                vale, 0, mailvo, .false._1)
    call jeveuo(nomob, 'L', jm3d)
!
    do ima = 1, nbmail
        if (zi(jm3d-1+ima) .eq. 0) then
            call utmess('F', 'XFEM2_59')
        end if
    end do
!
!     ------------------------------------------------------------------
!     III°) CREATION DU VECTEUR DES NORMALES SORTANTES
!     ------------------------------------------------------------------
!
    vecnor = '&&XORIPE.VECNOR'
    call wkvect(vecnor, 'V V R', nbmail*ndim, jvecno)
!
    do ima = 1, nbmail
!       NUMEROS DES MAILLES PRINCIPALE ET DE BORD
        numab = listCellNume(ima)
        numapr = zi(jm3d-1+ima)
!
!       NOMBRES DE NOEUDS DES MAILLES PRINCIPALE ET DE BORD
        nbnobo = zi(jconx2+numab)-zi(jconx2+numab-1)
        nbnopr = zi(jconx2+numapr)-zi(jconx2+numapr-1)
!
!       GBO : CENTRE DE GRAVITÉ DE LA MAILLE DE BORD
        gbo = 0.d0
        do ino = 1, nbnobo
            nuno = connex(zi(jconx2+numab-1)+ino-1)
            do j = 1, ndim
                gbo(j) = gbo(j)+vale(3*(nuno-1)+j)/nbnobo
            end do
        end do
!
!     GPR : CENTRE DE GRAVITÉ DE LA MAILLE PRICIPALE
        gpr = 0.d0
        do ino = 1, nbnopr
            nuno = connex(zi(jconx2+numapr-1)+ino-1)
            do j = 1, ndim
                gpr(j) = gpr(j)+vale(3*(nuno-1)+j)/nbnopr
            end do
        end do
!
!       NORMALE EXTERIEURE : NEXT = GBO - GPR
        next = gbo-gpr
        call normev(next, norme)
!
        do j = 1, ndim
            zr(jvecno-1+ndim*(ima-1)+j) = next(j)
        end do
!
    end do
!
!
!
!
!     ------------------------------------------------------------------
!     IV°) ORIENTATION DES SOUS-ELEMENTS
!     ------------------------------------------------------------------
!
    chs(1) = '&&XORIPE.PINTTO'
    chs(2) = '&&XORIPE.CNSETO'
    chs(3) = '&&XORIPE.LONCHA'
    chs(4) = '&&XORIPE.HEAV'
!
    pintto = modele//'.TOPOSE.PIN'
    cnseto = modele//'.TOPOSE.CNS'
    loncha = modele//'.TOPOSE.LON'
    heav = modele//'.TOPOSE.HEA'
!
    call celces(pintto, 'V', chs(1))
    call celces(cnseto, 'V', chs(2))
    call celces(loncha, 'V', chs(3))
    call celces(heav, 'V', chs(4))
!
    do ich = 1, 4
        call jeveuo(chs(ich)//'.CESD', 'L', jcesd(ich))
        call jeveuo(chs(ich)//'.CESV', 'E', jcesv(ich))
        call jeveuo(chs(ich)//'.CESL', 'L', jcesl(ich))
    end do
!
    do ima = 1, nbmail
        do j = 1, ndim
            next(j) = zr(jvecno-1+ndim*(ima-1)+j)
        end do
!
        numab = listCellNume(ima)
! --- CA NE SERT A RIEN DE RECUPERER NDIME CAR ON A SELECTIONNÉ NUMAB
! --- TEL QUE NDIME = NDIM-1 (BOUCLE 120)
        ndime = tmdim(typmail(numab))
        nfiss = zi(jcesd(4)-1+5+4*(numab-1)+2)
        numapr = zi(jm3d-1+ima)
        nbnopr = zi(jconx2+numapr)-zi(jconx2+numapr-1)
!
        itypma = typmail(numapr)
        call panbno(itypma, nbnott)
        quadratique = .false.
        if (ndim .eq. 2) then
            itypbo = typmail(numab)
            call jenuno(jexnum('&CATA.TM.NOMTM', itypbo), typbo)
            nbnobo = zi(jconx2+numab)-zi(jconx2+numab-1)
            nbnose = nbnobo
            nbnos = nbnobo
        else
            itypbo = typmail(numab)
            call jenuno(jexnum('&CATA.TM.NOMTM', itypbo), typbo)
            nbnobo = zi(jconx2+numab)-zi(jconx2+numab-1)
            if (nbnobo .gt. 4) then
                nbnose = 6
                quadratique = .true.
            else
                nbnose = 3
                quadratique = .false.
            end if
            nbnos = 3
        end if
!
!       RECUPERATION DE LA SUBDIVISION LA MAILLE DE PEAU EN NIT
!       SOUS-ELEMENTS
        call cesexi('S', jcesd(3), jcesl(3), numab, 1, &
                    1, 1, iad)
        nse = zi(jcesv(3)-1+iad)
!
!         BOUCLE SUR LES NSE SOUS-ELEMENTS
        do ise = 1, nse
!
!         CO(J,IN) : JEME COORDONNEE DU INEME SOMMET DU SOUS-ELEMENT
            do in = 1, nbnos
                icmp = nbnose*(ise-1)+in
                call cesexi('S', jcesd(2), jcesl(2), numab, 1, &
                            1, icmp, id(in))
                ino = zi(jcesv(2)-1+id(in))
                if (ino .lt. 1000) then
                    nuno = connex(zi(jconx2+numab-1)+ino-1)
!
                    do j = 1, ndim
                        co(j, in) = vale(3*(nuno-1)+j)
                    end do
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    do j = 1, ndim
                        icmp = ndim*(ino-1000-1)+j
                        call cesexi('S', jcesd(1), jcesl(1), numab, 1, &
                                    1, icmp, iad)
                        co(j, in) = zr(jcesv(1)-1+iad)
                    end do
                end if
            end do
!
            do j = 1, ndim
                a(j) = co(j, 1)
                b(j) = co(j, 2)
                if (ndim .eq. 3) c(j) = co(j, 3)
            end do
            if (ndim .eq. 2) then
                a(3) = 0.d0
                b(3) = 0.d0
            end if
!
!         NORMALE AU SOUS-ELEMENT 2D
            ab = 0.d0
            ab = b-a
            n2d = 0.d0
            if (ndim .eq. 3) then
                ac = 0.d0
                ac = c-a
                call provec(ab, ac, n2d)
            else if (ndim .eq. 2) then
                n2d(1) = ab(2)
                n2d(2) = -ab(1)
                n2d(3) = 0
            end if
            call normev(n2d, norme)
!
!
!         PRODUIT SCALAIRE DES NORMALES : N2D.NEXT
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            if (ddot(b_n, n2d, b_incx, next, b_incy) .lt. 0.d0) then
!           ON INVERSE LES SOMMETS S1 ET S2
!           (ON INVERSE 1 ET 2 EN 2D
                s1 = ndim-1
                s2 = ndim
!            ON INVERSE 2 ET 3 EN 3D)
                nseori = nseori+1
                intemp = zi(jcesv(2)-1+id(s1))
                zi(jcesv(2)-1+id(s1)) = zi(jcesv(2)-1+id(s2))
                zi(jcesv(2)-1+id(s2)) = intemp
                if (quadratique) then
!            ON INVERSE LES NOEUDS MILIEUX 4 ET 6 EN 3D)
!               RECUPERATION DE LA CONNECTIVITE :
                    icmp = nbnose*(ise-1)+4
                    call cesexi('S', jcesd(2), jcesl(2), numab, 1, &
                                1, icmp, id4)
                    icmp = nbnose*(ise-1)+6
                    call cesexi('S', jcesd(2), jcesl(2), numab, 1, &
                                1, icmp, id6)
!               INVERSION  DE LA CONNECTIVITE :
                    intemp = zi(jcesv(2)-1+id4)
                    zi(jcesv(2)-1+id4) = zi(jcesv(2)-1+id6)
                    zi(jcesv(2)-1+id6) = intemp
                end if
            end if
!
!         ON MODIFIE HEAVISIDE SI BORD COINCIDANT AVEC INTERFACE
!         RECUPERATION DE LA VALEUR DE LA FONCTION HEAVISIDE
            do ifiss = 1, nfiss
                ihe = ise
                call cesexi('S', jcesd(4), jcesl(4), numab, 1, &
                            ifiss, ihe, iad)
                he = zi(jcesv(4)-1+iad)
                ASSERT(he .eq. -1 .or. he .eq. 1 .or. he .eq. 0 .or. he .eq. 99)
                if (he .eq. 99) then
!             VERIF QUE C'EST NORMAL
                    if ((typbo(1:4) .eq. 'TRIA') .or. (typbo(1:3) .eq. 'SEG')) then
                        ASSERT(nse .eq. 1)
                    else if (typbo(1:4) .eq. 'QUAD') then
                        ASSERT(nse .eq. 2)
                    else
                        ASSERT(.false.)
                    end if
!             SIGNE LEVEL SET SUR LA MAILLE PRINCIPALE
                    nsignp = 0
                    nsignm = 0
                    nsignz = 0
!             LSN SUR LES NOEUDS SOMMETS DE LA MAILLE PRINCIPALE
                    do ino = 1, nbnott(1)
                        call cesexi('S', jlsnd, jlsnl, numapr, ino, &
                                    ifiss, 1, iad2)
!                NUNO=ZI(JCONX1-1+ZI(JCONX2+NUMAPR-1)+INO-1)
                        lsn = cesv(iad2)
                        if (lsn .gt. 0.d0) nsignp = nsignp+1
                        if (lsn .eq. 0.d0) nsignz = nsignz+1
                        if (lsn .lt. 0.d0) nsignm = nsignm+1
                    end do
                    ASSERT(nsignz .ne. 0)
!             REMARQUE : LES DEUX TESTS SUIVANTS NE SONT PAS CORRECTS
!             VOIR FICHE 13265
!              ASSERT(NSIGNP+NSIGNM.NE.0)
!              ASSERT(NSIGNP*NSIGNM.EQ.0)
!             ON ECRIT HE
                    if (nsignp .gt. 0) zi(jcesv(4)-1+iad) = 1
                    if (nsignm .gt. 0) zi(jcesv(4)-1+iad) = -1
                end if
            end do
!
        end do
!
    end do
!
!     ON SAUVE LES NOUVEAUX CHAM_ELEM MODIFIES A LA PLACE DES ANCIENS
    call cescel(chs(2), ligrel, 'TOPOSE', 'PCNSETO', 'OUI', &
                nncp, 'G', cnseto, 'F', ibid)
    call cescel(chs(4), ligrel, 'TOPOSE', 'PHEAVTO', 'OUI', &
                nncp, 'G', heav, 'F', ibid)
!     ------------------------------------------------------------------
!     FIN
!     ------------------------------------------------------------------
!
    call jedetr('&&XORIPE.NU_MAILLE_3D')
    call jedetr('&&XORIPE.VECNOR')
!
999 continue
!
    write (ifm, *) 'NOMBRE DE SOUS-ELEMENTS DE PEAU RE-ORIENTES :', nseori
!
    call jedetr('&&XORIPE.GRMAPE')
!
    call jedema()
end subroutine
