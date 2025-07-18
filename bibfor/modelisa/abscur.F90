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
subroutine abscur(ma)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/arcseg34.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/imprsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
#include "asterfort/sdmail.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: ma
!-----------------------------------------------------------------------
!     calcul d'une abscisse curviligne pour un groupe de mailles :
!     toutes les mailles doivent etre du type 'poi1' ou 'seg2,3,4'
!
!     arguments en entree
!     ------------------
!     ma   : nom du maillage
!
!     en sortie
!     ---------
!     creation d'une carte contenant l'abscisse curviligne
!-----------------------------------------------------------------------
!
    character(len=8) :: typm, noma, nono
    character(len=24) :: mesmai, mesnoe
    character(len=24) :: cooval, coodsc, grpnoe
    character(len=24) :: gpptnn, grpmai, gpptnm, connex, titre, typmai, adapma
    character(len=16) :: motcle(3), typmcl(3)
    integer(kind=8) :: adrm, iseg1, iseg2, isegprev, jtmp, kseg, nbextr, nbnot
    integer(kind=8) :: iab1, iab2, iadr2, numa2, nuno1, nuno2
    integer(kind=8) :: icoo1, icoo2, icoo3, icoo4, icor2, kma, nunosuiv, vuorig
    integer(kind=8) :: jpoi, jseg, ino
    integer(kind=8) :: ipoi1, iseg, itypm, jcoor
    integer(kind=8) :: mi, n, n1, n2, nunorig
    integer(kind=8) :: nbpoi1, nbma, nbseg, nbno
    integer(kind=8) :: jmesma, jmesno, numa, iexi, nbnoseg
    real(kind=8) :: s, stot
    real(kind=8) :: s13, s32, s34, s42, abscurv(4), coor(3, 4)
    aster_logical, parameter :: dbg = .false.
    integer(kind=8), pointer :: icoseg(:) => null()
    integer(kind=8), pointer :: nu2seg(:) => null()
    integer(kind=8), pointer :: segordo(:) => null()
    real(kind=8), pointer :: valv(:) => null()
    character(len=8), pointer :: ncmp(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    call sdmail(ma, cooval, coodsc, &
                grpnoe, gpptnn, grpmai, gpptnm, &
                connex, titre, typmai, adapma)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbnot)
    call jeveuo(cooval, 'L', jcoor)
!
!
!   -- 1. numero du noeud origine :
!   -------------------------------
    mesnoe = '&&ABSCUR.MESNOE'
    motcle(1) = 'GROUP_NO_ORIG'
    motcle(2) = 'NOEUD_ORIG'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
    call reliem(' ', ma, 'NU_NOEUD', 'ABSC_CURV', 1, &
                2, motcle, typmcl, mesnoe, nbno)
    if (nbno .ne. 1) call utmess('F', 'MODELISA_5', sk=motcle(1))
    call jeveuo(mesnoe, 'L', jmesno)
    nunorig = zi(jmesno)
!
!
!   -- 2. liste des mailles a traiter :
!   -----------------------------------
    mesmai = '&&ABSCUR.MESMAI'
    motcle(1) = 'TOUT'
    motcle(2) = 'GROUP_MA'
    motcle(3) = 'MAILLE'
    typmcl(1) = 'TOUT'
    typmcl(2) = 'GROUP_MA'
    typmcl(3) = 'MAILLE'
    call reliem(' ', ma, 'NU_MAILLE', 'ABSC_CURV', 1, &
                3, motcle, typmcl, mesmai, nbma)
    call jeveuo(mesmai, 'L', jmesma)
!
!
!   --  3. creation d'objets temporaires :
!   ---------------------------------------
    call wkvect('&&ABSCUR.IPOI1', 'V V I', nbma, jpoi)
    call wkvect('&&ABSCUR.SEG', 'V V I', nbma, jseg)
    call wkvect('&&ABSCUR.AB1  ', 'V V R', nbma, iab1)
    call wkvect('&&ABSCUR.AB2  ', 'V V R', nbma, iab2)
    call wkvect('&&ABSCUR.COR2 ', 'V V I', nbma, icor2)
!
!
!   --  4. tri des mailles poi1 et seg
!   -----------------------------------
    call jeveuo(typmai, 'L', itypm)
    nbseg = 0
    nbpoi1 = 0
    do kma = 1, nbma
        numa = zi(jmesma-1+kma)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(itypm+numa-1)), typm)
        if (typm .eq. 'SEG2' .or. typm .eq. 'SEG3' .or. typm .eq. 'SEG4') then
            nbseg = nbseg+1
            zi(jseg+nbseg-1) = numa
        else if (typm .eq. 'POI1') then
            nbpoi1 = nbpoi1+1
            zi(jpoi+nbpoi1-1) = numa
        else
            call utmess('F', 'MODELISA_2')
        end if
    end do
!
!
!   --  5. Les segments doivent former une ligne ouverte avec
!          deux extremites.
!   --------------------------------------------------------
!   -- 5.1 : on note les noeuds extremites des segments :
    AS_ALLOCATE(vi=icoseg, size=nbnot)
    AS_ALLOCATE(vi=nu2seg, size=2*nbnot)
    do iseg = 1, nbseg
        numa = zi(jseg-1+iseg)
        call jeveuo(jexnum(connex, numa), 'L', jtmp)
        nuno1 = zi(jtmp-1+1)
        nuno2 = zi(jtmp-1+2)
        icoseg(nuno1) = icoseg(nuno1)+1
        icoseg(nuno2) = icoseg(nuno2)+1
        if (icoseg(nuno1) .le. 2) nu2seg(2*(nuno1-1)+icoseg(nuno1)) = iseg
        if (icoseg(nuno2) .le. 2) nu2seg(2*(nuno2-1)+icoseg(nuno2)) = iseg
    end do
!
!   -- 5.2 : on verifie que les segments forment une ligne ouverte
!            et que l'une des 2 extremites est nunorig :
    vuorig = 0
    nbextr = 0
    do ino = 1, nbnot
        if (icoseg(ino) .eq. 2) then
        else if (icoseg(ino) .eq. 1) then
            nbextr = nbextr+1
            if (ino .eq. nunorig) then
                vuorig = vuorig+1
            end if
        else if (icoseg(ino) .eq. 0) then
        else
            nono = int_to_char8(ino)
            call utmess('F', 'INTEMAIL_36', sk=nono)
        end if
    end do
    if (nbextr .ne. 2) call utmess('F', 'INTEMAIL_37', si=nbextr)
    if (vuorig .ne. 1) call utmess('F', 'INTEMAIL_38')
!
!   -- 5.3 : on etablit la liste ordonnee des segments
!   ---------------------------------------------------
    AS_ALLOCATE(vi=segordo, size=nbseg)
!
!   -- on initialise la recherche (1er segment):
    kseg = 1
    nunosuiv = nunorig
    iseg = nu2seg(2*nunosuiv-1)
    ASSERT(iseg .gt. 0)
    numa = zi(jseg-1+iseg)
    call jeveuo(jexnum(connex, numa), 'L', jtmp)
    ASSERT(nunosuiv .eq. zi(jtmp) .or. nunosuiv .eq. zi(jtmp+1))
    if (nunosuiv .eq. zi(jtmp)) then
        segordo(kseg) = +iseg
        nunosuiv = zi(jtmp+1)
    else
        segordo(kseg) = -iseg
        nunosuiv = zi(jtmp)
    end if
    isegprev = iseg
!
!   -- on met les segments bout a bout :
    do while (kseg .lt. nbseg)
        kseg = kseg+1
        iseg1 = nu2seg(2*nunosuiv-1)
        iseg2 = nu2seg(2*nunosuiv)
        if (iseg1 .eq. isegprev) then
            iseg = iseg2
        else if (iseg2 .eq. isegprev) then
            iseg = iseg1
        else
            ASSERT(.false.)
        end if
        numa = zi(jseg-1+iseg)
        call jeveuo(jexnum(connex, numa), 'L', jtmp)
        ASSERT(nunosuiv .eq. zi(jtmp) .or. nunosuiv .eq. zi(jtmp+1))
        if (nunosuiv .eq. zi(jtmp)) then
            segordo(kseg) = +iseg
            nunosuiv = zi(jtmp+1)
        else
            segordo(kseg) = -iseg
            nunosuiv = zi(jtmp)
        end if
        isegprev = iseg
    end do
!
!
!   -- 6. On verifie que les POI1 sont sur des segments :
!   --------------------------------------------------------
    do ipoi1 = 1, nbpoi1
        numa = zi(jpoi+ipoi1-1)
        call jeveuo(jexnum(connex, numa), 'L', adrm)
        n = zi(adrm)
        do iseg = 1, nbseg
            numa2 = zi(jseg-1+iseg)
            call jeveuo(jexnum(connex, numa2), 'L', iadr2)
            n1 = zi(iadr2)
            n2 = zi(iadr2+1)
            if (n1 .eq. n) then
                zi(icor2+ipoi1-1) = iseg
                goto 15
            else if (n2 .eq. n) then
                zi(icor2+ipoi1-1) = -iseg
                goto 15
            end if
        end do
        noma = int_to_char8(numa)
        call utmess('F', 'MODELISA_3', sk=noma)
15      continue
    end do
!
!
!
!   --  7. allocation de la carte :
!   ------------------------------
    call exisd('CHAMP', ma//'.ABSC_CURV', iexi)
    if (iexi .eq. 1) then
        call utmess('A', 'INTEMAIL_35')
        call detrsd('CHAMP', ma//'.ABSC_CURV')
    end if
    call alcart('G', ma//'.ABSC_CURV', ma, 'ABSC_R')
    call jeveuo(ma//'.ABSC_CURV .NCMP', 'E', vk8=ncmp)
    call jeveuo(ma//'.ABSC_CURV .VALV', 'E', vr=valv)
    ncmp(1) = 'ABSC1'
    ncmp(2) = 'ABSC2'
    ncmp(3) = 'ABSC3'
    ncmp(4) = 'ABSC4'
    stot = 0.d0
!
!
!   --  8. calcul de l'abscisse curviligne
!   ---------------------------------------
    do kseg = 1, nbseg
        mi = segordo(kseg)
        numa = zi(jseg-1+abs(mi))
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(itypm+numa-1)), typm)
        call jeveuo(jexnum(connex, numa), 'L', jtmp)
!
!       noeuds 1 et 2 (extremites) :
        icoo1 = 3*(zi(jtmp-1+1)-1)
        icoo2 = 3*(zi(jtmp-1+2)-1)
        coor(:, 1) = zr(jcoor+icoo1-1+1:jcoor+icoo1-1+3)
        coor(:, 2) = zr(jcoor+icoo2-1+1:jcoor+icoo2-1+3)
!
!       noeuds 3 et 4 (si necessaire) :
        if (typm .eq. 'SEG3') then
            nbnoseg = 3
            icoo3 = 3*(zi(jtmp-1+3)-1)
            coor(:, 3) = zr(jcoor+icoo3-1+1:jcoor+icoo3-1+3)
        else if (typm .eq. 'SEG4') then
            nbnoseg = 4
            icoo3 = 3*(zi(jtmp-1+3)-1)
            icoo4 = 3*(zi(jtmp-1+4)-1)
            coor(:, 3) = zr(jcoor+icoo3-1+1:jcoor+icoo3-1+3)
            coor(:, 4) = zr(jcoor+icoo4-1+1:jcoor+icoo4-1+3)
        else
            ASSERT(typm .eq. 'SEG2')
            nbnoseg = 2
        end if
!
!       -- calcul des abscisses curvilignes :
        call arcseg34(nbnoseg, coor, abscurv)
!
        if (nbnoseg .eq. 2) then
            s = abscurv(2)-abscurv(1)
        else if (nbnoseg .eq. 3) then
            s13 = abscurv(3)-abscurv(1)
            s32 = abscurv(2)-abscurv(3)
            s = s13+s32
        else if (nbnoseg .eq. 4) then
            s13 = abscurv(3)-abscurv(1)
            s34 = abscurv(4)-abscurv(3)
            s42 = abscurv(2)-abscurv(4)
            s = s13+s34+s42
        end if
!
        if (mi .gt. 0) then
            valv(1) = stot
            zr(iab1+abs(mi)-1) = stot
!
            if (nbnoseg .eq. 3) then
                valv(3) = stot+s13
            else if (nbnoseg .eq. 4) then
                valv(3) = stot+s13
                valv(4) = stot+s13+s34
            end if
!
            stot = stot+s
            valv(2) = stot
            zr(iab2+abs(mi)-1) = stot
        else
            valv(2) = stot
            zr(iab2+abs(mi)-1) = stot
!
            if (nbnoseg .eq. 3) then
                valv(3) = stot+s32
            else if (nbnoseg .eq. 4) then
!
                valv(4) = stot+s42
                valv(3) = stot+s42+s34
            end if
!
            stot = stot+s
            valv(1) = stot
            zr(iab1+abs(mi)-1) = stot
        end if
        call nocart(ma//'.ABSC_CURV', 3, nbnoseg, mode='NUM', nma=1, &
                    limanu=[numa])
!
    end do
!
!
!   --  cas des poi1 :
!   --------------------
    do ipoi1 = 1, nbpoi1
        numa = zi(jpoi+ipoi1-1)
        mi = zi(icor2+ipoi1-1)
        if (mi .gt. 0) then
            s = zr(iab1+mi-1)
        else
            s = zr(iab2-mi-1)
        end if
        valv(1) = s
        call nocart(ma//'.ABSC_CURV', 3, 1, mode='NUM', nma=1, &
                    limanu=[numa])
    end do
!
!
    if (dbg) call imprsd('CHAMP', ma//'.ABSC_CURV', 6, 'ABSC_CURV')
!
!
!   -- menage
!   ---------
    AS_DEALLOCATE(vi=icoseg)
    AS_DEALLOCATE(vi=nu2seg)
    AS_DEALLOCATE(vi=segordo)
    call jedetr('&&ABSCUR.AB1')
    call jedetr('&&ABSCUR.AB2')
    call jedetr('&&ABSCUR.COR2')
    call jedetr('&&ABSCUR.IPOI1')
    call jedetr('&&ABSCUR.SEG')
    call jedetr(mesnoe)
    call jedetr(mesmai)
    call jedetr(ma//'.ABSC_CURV .NCMP')
    call jedetr(ma//'.ABSC_CURV .VALV')
!
    call jedema()
end subroutine
