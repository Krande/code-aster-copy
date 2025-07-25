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

subroutine pj5dco(mo1, mo2, corres)
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/dismoi.h"
#include "asterfort/exmano.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/pacoa2.h"
#include "asterfort/pj3da4.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=16) :: corres
    character(len=8) :: mo1, mo2
!     BUT :
!       CREER UNE SD CORRESP_2_MAILLA
!       DONNANT LA CORRESPONDANCE ENTRE LES NOEUDS DU MAILLAGE M1
!       ET CEUX DE M2 (DANS LE CAS DE MAILLAGE EN SEG2)

!  IN/JXIN   MO1      K8  : NOM DU MODELE INITIAL
!  IN/JXIN   MO2      K8  : NOM DU MODELE SUR LEQUEL ON VEUT PROJETER
!                           DES CHAMPS
!  IN/JXOUT  CORRES  K16 : NOM DE LA SD CORRESP_2_MAILLA

! ----------------------------------------------------------------------

    integer(kind=8) :: nbmail, nbdim, nno1, nno2, nbno
    integer(kind=8) :: llin1, llin2, inode, nbtr, lno1, lno2, lco1, lco2, idecal
    integer(kind=8) :: out1, out2, jcoo1, jcoo2, ilcnx1, numnoe
    integer(kind=8) :: i, nbmano, ima, imail, ino, j2xxk1, i2conb, i2conu, i2cocf
    parameter(nbmail=10)
    parameter(nbdim=3)
    integer(kind=8) :: numano(nbmail), nunoe(nbmail)
    real(kind=8) :: a(nbdim), b(nbdim), m(nbdim), un, deux
    real(kind=8) :: dpmin, l1, l2, xabs, dp, am(nbdim), bm(nbdim), a1, a2, dist
    character(len=8) :: m1, m2
    character(len=16) :: lisin1, lisin2, lisou1, lisou2
    character(len=16) :: noeud1, noeud2, cobar1, cobar2
    character(len=24) :: coormo, coorme
    integer(kind=8), pointer :: connex(:) => null()

! DEB ------------------------------------------------------------------
    call jemarq()

    un = 1.d0
    deux = 2.d0

    call dismoi('NOM_MAILLA', mo1, 'MODELE', repk=m1)
    call dismoi('NOM_MAILLA', mo2, 'MODELE', repk=m2)

    call dismoi('NB_NO_MAILLA', m1, 'MAILLAGE', repi=nno1)
    call dismoi('NB_NO_MAILLA', m2, 'MAILLAGE', repi=nno2)

    if (nno2 .eq. 0) then
        call utmess('F', 'PROJECTION4_54')
    end if

!     DETERMINATION DE LA DIMENSION DE L'ESPACE :
!     --------------------------------------------------------

!     Initialisation des A, B et M
    do i = 1, nbdim
        a(i) = 0.d0
        b(i) = 0.d0
        m(i) = 0.d0
        am(i) = 0.d0
        bm(i) = 0.d0
    end do

    coormo = m1//'.COORDO    .VALE'
    call jeveuo(coormo, 'L', jcoo1)

    coorme = m2//'.COORDO    .VALE'
    call jeveuo(coorme, 'L', jcoo2)

!     1. RECHECHE DES MAILLES LIEES AU NOEUD LE PLUS PROCHE
!     ------------------------------------------------
    lisin1 = 'NOEUD_MODELE'
    call wkvect(lisin1, 'V V K8', nno1, llin1)

    lisin2 = 'NOEUD_MESURE'
    call wkvect(lisin2, 'V V K8', nno2, llin2)

    lisou1 = 'NOEUD_MODELE_VIS'
    lisou2 = 'NOEUD_MESURE_VIS'

    do inode = 1, nno1
        zk8(llin1-1+inode) = int_to_char8(inode)
    end do
    do inode = 1, nno2
        zk8(llin2-1+inode) = int_to_char8(inode)
    end do

    call pacoa2(lisin1, lisin2, nno1, nno2, m1, &
                m2, lisou1, lisou2, nbtr)
    if (nbtr .ne. nno2) then
        call utmess('F', 'ALGORITH9_91')
    end if

    noeud1 = 'NOEUD_DEBUT'
    call wkvect(noeud1, 'V V I', nbtr, lno1)

    noeud2 = 'NOEUD_FIN'
    call wkvect(noeud2, 'V V I', nbtr, lno2)

    cobar1 = 'COEFF_DEBUT'
    call wkvect(cobar1, 'V V R', nbtr, lco1)

    cobar2 = 'COEFF_FIN'
    call wkvect(cobar2, 'V V R', nbtr, lco2)

    call jeveuo(lisou1, 'L', out1)
    call jeveuo(lisou2, 'L', out2)

    call jeveuo(m1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(m1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)

!     2. RECHECHE DE LA MAILLE LE PLUS PROCHE DU NOEUD MESURE
!     ------------------------------------------------
    do inode = 1, nbtr
        numnoe = char8_to_int(zk8(out2-1+inode))
        do i = 1, nbdim
            m(i) = zr(jcoo2-1+(numnoe-1)*nbdim+i)
        end do

        numnoe = char8_to_int(zk8(out1-1+inode))
        call exmano(m1, numnoe, numano, nbmano)
        if (nbmano .eq. 0) then
            call utmess('F', 'ALGORITH9_92')
        end if
        dpmin = r8maem()
        do ima = 1, nbmano
            imail = numano(ima)
            nbno = zi(ilcnx1-1+imail+1)-zi(ilcnx1-1+imail)
!     THEORIQUEMENT NBNO = 2 (POUR SEG2)
            do ino = 1, nbno
                nunoe(ino) = connex(zi(ilcnx1-1+imail)-1+ino)
            end do
            do i = 1, nbdim
                a(i) = zr(jcoo1-1+(nunoe(1)-1)*nbdim+i)
                b(i) = zr(jcoo1-1+(nunoe(2)-1)*nbdim+i)
            end do

!     3. CALCUL DE LA DISTANCE NOEUD-MAILLE (AM + BM)
!     ------------------------------------------------
            do i = 1, nbdim
                am(i) = m(i)-a(i)
                bm(i) = m(i)-b(i)
            end do

            a1 = 0.d0
            a2 = 0.d0
            do i = 1, nbdim
                a1 = a1+am(i)*am(i)
                a2 = a2+bm(i)*bm(i)
            end do
            dist = sqrt(a1)+sqrt(a2)

            if (dist .lt. dpmin) then
                dpmin = dist

!     4. CALCUL DES COORDONNEES BARYCENTRIQUES
!     ------------------------------------------------
                call pj3da4(m, a, b, l1, l2, &
                            dp)
                zi(lno1-1+inode) = nunoe(1)
                zi(lno2-1+inode) = nunoe(2)
                zr(lco1-1+inode) = l1
                zr(lco2-1+inode) = l2
            end if
        end do

!     5. APPLICATION FONCTION DE FORME (ELEMENT ISOPARAMETRIQUE)
!     ---------------------------------------------------
!     XABS : ABSCISSE DU POINT DANS L'ELEMENT DE REFERENCE (SEG2)
        xabs = -un*zr(lco1-1+inode)+un*zr(lco2-1+inode)
        zr(lco1-1+inode) = (un-xabs)/deux
        zr(lco2-1+inode) = (un+xabs)/deux
    end do

!     6. CREATION DE LA SD CORRESP_2_MAILLA : CORRES
!     ---------------------------------------------------
    call wkvect(corres//'.PJXX_K1', 'V V K24', 5, j2xxk1)
    call wkvect(corres//'.PJEF_NB', 'V V I', nno2, i2conb)
    call wkvect(corres//'.PJEF_NU', 'V V I', nno2*2, i2conu)
    call wkvect(corres//'.PJEF_CF', 'V V R', nno2*2, i2cocf)

    zk24(j2xxk1-1+1) = m1
    zk24(j2xxk1-1+2) = m2
    zk24(j2xxk1-1+3) = 'COLLOCATION'

    do ino = 1, nno2
        zi(i2conb-1+ino) = 2
    end do

    idecal = 0
    do ino = 1, nno2
        zi(i2conu-1+idecal+1) = zi(lno1-1+ino)
        zi(i2conu-1+idecal+2) = zi(lno2-1+ino)
        zr(i2cocf-1+idecal+1) = zr(lco1-1+ino)
        zr(i2cocf-1+idecal+2) = zr(lco2-1+ino)
        idecal = idecal+zi(i2conb-1+ino)
    end do

    call jedetr(lisin1)
    call jedetr(lisin2)
    call jedetr(noeud1)
    call jedetr(noeud2)
    call jedetr(cobar1)
    call jedetr(cobar2)

    call jedema()
end subroutine
