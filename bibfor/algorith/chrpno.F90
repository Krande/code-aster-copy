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
subroutine chrpno(champ1, repere, nom_cham, type)
! aslint: disable=W1501
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/getfac.h"
#include "asterfort/angvxy.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/reliem.h"
#include "asterfort/selectCompN.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvgl.h"
#include "blas/ddot.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: champ1, repere, type, nom_cham
! ----------------------------------------------------------------------
!
!     BUT : CHANGEMENT DE REPERE DANS LE CAS D'UN CHAM_NO
! ----------------------------------------------------------------------
!     ARGUMENTS :
!     CHAMP1   IN  K16  : NOM DU CHAMP A TRAITER
!     REPERE   IN  K16  : TYPE DE REPERE (UTILISATEUR OU CYLINDRIQUE)
! ----------------------------------------------------------------------
! ---------------------------------------------------------------------
!
    integer(kind=8) :: i, nbno, ino, ibid, nbcmp, ndim_type
    integer(kind=8) :: ii, nbma, ipt2, inel
    integer(kind=8) :: jcnsv, jconx1, jconx2, nbpt
    integer(kind=8) :: ipt, inot, ndim, licmpu(6), jcnsl
    integer(kind=8) :: nbn, idnoeu, nbnoeu, inoe
    real(kind=8) :: angnot(3), pgl(3, 3), valer(6), valed(6)
    real(kind=8) :: valr, valei(6)
    real(kind=8) :: orig(3), axez(3), axer(3), axet(3)
    real(kind=8) :: epsi, prosca, xnormr, valet(6)
    real(kind=8) :: vectx(3), vecty(3)
    complex(kind=8) :: valetc(6)
    character(len=3) :: tsca
    character(len=8) :: ma, k8b, typmcl(4), nomgd
    character(len=16) :: motcle(4)
    character(len=19) :: chams1, chams0
    character(len=24) :: mesnoe
    character(len=24) :: valk
    integer(kind=8), parameter :: nbCmpMax = 8
    character(len=8) :: nom_cmp(nbCmpMax)
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8) :: iocc, nocc
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
    epsi = 1.0d-6
    motcle(1) = 'GROUP_NO'
    typmcl(1) = 'GROUP_NO'
    motcle(2) = 'NOEUD'
    typmcl(2) = 'NOEUD'
    motcle(3) = 'GROUP_MA'
    typmcl(3) = 'GROUP_MA'
    motcle(4) = 'MAILLE'
    typmcl(4) = 'MAILLE'
    mesnoe = '&&CHRPEL.MES_NOEUDS'
!
! ----- DEFINITION ET CREATION DU CHAM_NO SIMPLE CHAMS1
! ----- A PARTIR DU CHAM_NO CHAMP1
!
    chams0 = '&&CHRPNO.CHAMS0'
    chams1 = '&&CHRPNO.CHAMS1'
    call cnocns(champ1, 'V', chams0)
!   sélection des composantes :
    call selectCompN(chams0, nom_cham, type, nbcmp, nom_cmp, &
                     ndim_type)
    call cnsred(chams0, 0, [0], nbcmp, nom_cmp, &
                'V', chams1)
    call detrsd('CHAM_NO_S', chams0)
    call jeveuo(chams1//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(chams1//'.CNSD', 'L', vi=cnsd)
    call jeveuo(chams1//'.CNSC', 'L', vk8=cnsc)
    ma = cnsk(1)
    nomgd = cnsk(2)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
    nbno = cnsd(1)
    call jeveuo(chams1//'.CNSV', 'E', jcnsv)
    call jeveuo(chams1//'.CNSL', 'E', jcnsl)
    call dismoi('Z_CST', ma, 'MAILLAGE', repk=k8b)
    ndim = 3
    if (k8b .eq. 'OUI') ndim = 2
!
    if (ndim .gt. ndim_type) then
        call utmess('F', 'ALGORITH12_45', sk=type)
    else if (ndim .lt. ndim_type) then
        call utmess('A', 'ALGORITH12_44', sk=type)
        ndim = 3
    end if
!
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
!
!
!
!  -- Le mot-clé AFFE définit les caractéristiques du nouveau repère
!     On peut définir un repère variable en définissant ces paramètres
!     par mailles/groupes de mailles
    call getfac('AFFE', nocc)
    do iocc = 1, nocc
! Construction de la liste des numéros de noeuds
! sélectionnées par les mots-clés GROUP_NO et NOEUD
        call reliem(' ', ma, 'NU_NOEUD', 'AFFE', iocc, &
                    4, motcle, typmcl, mesnoe, nbn)
!
        if (nbn .gt. 0) then
            nbnoeu = nbn
            call jeveuo(mesnoe, 'L', idnoeu)
        else
            nbnoeu = nbno
        end if
!
!    SI UNE DES COMPOSANTES EST ABSENTE SUR UN NOEUD
!    IL EST IGNORE
!
        do ino = 1, nbnoeu
            if (nbn .ne. 0) then
                inoe = zi(idnoeu+ino-1)
            else
                inoe = ino
            end if
            do ii = 1, nbcmp
                if (.not. zl(jcnsl-1+(inoe-1)*nbcmp+ii)) goto 112
            end do
            goto 110
112         continue
            do ii = 1, nbcmp
                zl(jcnsl-1+(inoe-1)*nbcmp+ii) = .false.
            end do
110         continue
        end do
!
        do i = 1, 6
            valed(i) = 0.0d0
            valer(i) = 0.0d0
            valet(i) = 0.0d0
        end do
        do i = 1, 3
            axer(i) = 0.0d0
            axet(i) = 0.0d0
            axez(i) = 0.0d0
            orig(i) = 0.0d0
            angnot(i) = 0.0d0
        end do
        licmpu(1) = 1
        licmpu(2) = 2
        licmpu(3) = 3
        licmpu(4) = 4
        licmpu(5) = 5
        licmpu(6) = 6
!
! ----- CHANGEMENT DE REPERE SUIVANT LE CHOIX UTILISATEUR
!
        if (repere(1:11) .eq. 'UTILISATEUR') then
!        SI LE NOUVEAU REPERE EST DONNE VIA DES VECTEURS
            call getvr8('AFFE', 'VECT_X', iocc=iocc, nbval=3, vect=vectx, &
                        nbret=ibid)
            if (ibid .ne. 0) then
                call getvr8('AFFE', 'VECT_Y', iocc=iocc, nbval=3, vect=vecty, &
                            nbret=ibid)
                if (ndim .ne. 3) then
                    call utmess('F', 'ALGORITH2_4')
                end if
                call angvxy(vectx, vecty, angnot)
            else
                if (ndim .eq. 3) then
                    call getvr8('AFFE', 'ANGL_NAUT', iocc=iocc, nbval=3, vect=angnot, &
                                nbret=ibid)
                    if (ibid .ne. 3) then
                        call utmess('F', 'ALGORITH2_7')
                    end if
                else
                    call getvr8('AFFE', 'ANGL_NAUT', iocc=iocc, scal=angnot(1), nbret=ibid)
                    if (ibid .ne. 1) then
                        valr = angnot(1)
                        call utmess('A', 'ALGORITH12_43', sr=valr)
                    end if
                end if
                angnot(1) = angnot(1)*r8dgrd()
                angnot(2) = angnot(2)*r8dgrd()
                angnot(3) = angnot(3)*r8dgrd()
            end if
            call matrot(angnot, pgl)
            if (type(1:4) .eq. 'TENS') then
                do ino = 1, nbnoeu
                    if (nbn .ne. 0) then
                        inoe = zi(idnoeu+ino-1)
                    else
                        inoe = ino
                    end if
                    if (.not. zl(jcnsl-1+(inoe-1)*nbcmp+1)) goto 10
                    if (tsca .eq. 'R') then
! CHAMP REEL
                        do ii = 1, nbcmp
                            valet(ii) = zr(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        valed(1) = valet(1)
                        valed(2) = valet(4)
                        valed(3) = valet(2)
                        valed(4) = valet(5)
                        valed(5) = valet(6)
                        valed(6) = valet(3)
                        call utpsgl(1, 3, pgl, valed, valet)
                        valer(1) = valet(1)
                        valer(2) = valet(3)
                        valer(3) = valet(6)
                        valer(4) = valet(2)
                        valer(5) = valet(4)
                        valer(6) = valet(5)
                        do ii = 1, nbcmp
                            zr(jcnsv-1+(inoe-1)*nbcmp+ii) = valer(ii)
                        end do
                    else
! CHAMP COMPLEXE
                        do ii = 1, nbcmp
                            valetc(ii) = zc(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        valed(1) = dble(valetc(1))
                        valed(2) = dble(valetc(4))
                        valed(3) = dble(valetc(2))
                        valed(4) = dble(valetc(5))
                        valed(5) = dble(valetc(6))
                        valed(6) = dble(valetc(3))
                        call utpsgl(1, 3, pgl, valed, valet)
                        valer(1) = valet(1)
                        valer(2) = valet(3)
                        valer(3) = valet(6)
                        valer(4) = valet(2)
                        valer(5) = valet(4)
                        valer(6) = valet(5)
!
                        valed(1) = dimag(valetc(1))
                        valed(2) = dimag(valetc(4))
                        valed(3) = dimag(valetc(2))
                        valed(4) = dimag(valetc(5))
                        valed(5) = dimag(valetc(6))
                        valed(6) = dimag(valetc(3))
                        call utpsgl(1, 3, pgl, valed, valet)
                        valei(1) = valet(1)
                        valei(2) = valet(3)
                        valei(3) = valet(6)
                        valei(4) = valet(2)
                        valei(5) = valet(4)
                        valei(6) = valet(5)
                        do ii = 1, nbcmp
                            zc(jcnsv-1+(inoe-1)*nbcmp+ii) = dcmplx(valer(ii), valei(ii))
                        end do
!
                    end if
!
10                  continue
                end do
            else
! VECTEUR
                do ino = 1, nbnoeu
                    if (nbn .ne. 0) then
                        inoe = zi(idnoeu+ino-1)
                    else
                        inoe = ino
                    end if
                    if (.not. zl(jcnsl-1+(inoe-1)*nbcmp+1)) goto 13
                    if (tsca .eq. 'R') then
! VECTEUR REEL
                        do ii = 1, nbcmp
                            valed(ii) = zr(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        if (ndim .eq. 3) then
                            call utpvgl(1, 3, pgl, valed, valer)
                            if (nbcmp .gt. 3) call utpvgl(1, 3, pgl, valed(4), valer(4))
                        else
                            call ut2vgl(1, 2, pgl, valed, valer)
                        end if
                        do ii = 1, nbcmp
                            zr(jcnsv-1+(inoe-1)*nbcmp+ii) = valer(ii)
                        end do
                    else
! VECTEUR COMPLEXE
                        do ii = 1, nbcmp
                            valetc(ii) = zc(jcnsv-1+(inoe-1)*nbcmp+ii)
                            valed(ii) = dble(valetc(ii))
                            valet(ii) = dimag(valetc(ii))
                        end do
                        if (ndim .eq. 3) then
                            call utpvgl(1, 3, pgl, valed, valer)
                            call utpvgl(1, 3, pgl, valet, valei)
                            if (nbcmp .gt. 3) then
                                call utpvgl(1, 3, pgl, valed(4), valer(4))
                                call utpvgl(1, 3, pgl, valet(4), valei(4))
                            end if
                        else
                            call ut2vgl(1, 2, pgl, valed, valer)
                            call ut2vgl(1, 2, pgl, valet, valei)
                        end if
                        do ii = 1, nbcmp
                            zc(jcnsv-1+(inoe-1)*nbcmp+ii) = dcmplx(valer(ii), valei(ii))
                        end do
!
                    end if
13                  continue
                end do
            end if
        else if (repere(1:5) .eq. 'COQUE') then
            call utmess('F', 'ALGORITH2_3')
        else
! REPERE CYLINDRIQUE
            if (ndim .eq. 3) then
                call getvr8('AFFE', 'ORIGINE', iocc=iocc, nbval=3, vect=orig, &
                            nbret=ibid)
                if (ibid .ne. 3) then
                    call utmess('F', 'ALGORITH2_8')
                end if
                call getvr8('AFFE', 'AXE_Z', iocc=iocc, nbval=3, vect=axez, &
                            nbret=ibid)
                if (ibid .eq. 0) then
                    call utmess('F', 'ALGORITH2_9')
                end if
            else
                call getvr8('AFFE', 'ORIGINE', iocc=iocc, nbval=2, vect=orig, &
                            nbret=ibid)
                if (ibid .ne. 2) then
                    call utmess('A', 'ALGORITH2_10')
                end if
                call getvr8('AFFE', 'AXE_Z', iocc=iocc, nbval=0, nbret=ibid)
                if (ibid .ne. 0) then
                    call utmess('A', 'ALGORITH2_11')
                end if
                axez(1) = 0.0d0
                axez(2) = 0.0d0
                axez(3) = 1.0d0
            end if
            xnormr = 0.0d0
            call normev(axez, xnormr)
            call jeveuo(ma//'.COORDO    .VALE', 'L', vr=vale)
            if (type(1:4) .eq. 'TENS') then
                if (ndim .eq. 2) then
                    licmpu(1) = 1
                    licmpu(2) = 2
                    licmpu(3) = 3
                    licmpu(4) = 5
                end if
!
                do ino = 1, nbnoeu
                    if (nbn .ne. 0) then
                        inoe = zi(idnoeu+ino-1)
                    else
                        inoe = ino
                    end if
                    axer(1) = vale(1+3*(inoe-1))-orig(1)
                    axer(2) = vale(1+3*(inoe-1)+1)-orig(2)
                    if (ndim .eq. 3) then
                        axer(3) = vale(1+3*(inoe-1)+2)-orig(3)
                    else
                        axer(3) = 0.0d0
                    end if
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    prosca = ddot(b_n, axer, b_incx, axez, b_incy)
                    axer(1) = axer(1)-prosca*axez(1)
                    axer(2) = axer(2)-prosca*axez(2)
                    if (ndim .eq. 3) then
                        axer(3) = axer(3)-prosca*axez(3)
                    else
                        axer(3) = 0.0d0
                    end if
                    xnormr = 0.0d0
                    call normev(axer, xnormr)
                    if (xnormr .lt. epsi) then
                        k8b = int_to_char8(inoe)
                        call utmess('A', 'ALGORITH2_13')
                        call jeveuo(ma//'.CONNEX', 'L', jconx1)
                        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', jconx2)
                        ipt2 = 0
                        do inel = 1, nbma
                            nbpt = zi(jconx2-1+inel+1)-zi(jconx2-1+inel)
                            do ipt = 1, nbpt
                                inot = zi(jconx1-1+zi(jconx2-1+inel)+ipt-1)
                                if (inot .eq. inoe) then
                                    axer(1) = 0.0d0
                                    axer(2) = 0.0d0
                                    axer(3) = 0.0d0
                                    do ipt2 = 1, nbpt
                                        inot = zi(jconx1-1+zi(jconx2-1+inel)+ipt2-1)
                                        axer(1) = axer(1)+vale(1+3*(inot-1))
                                        axer(2) = axer(2)+vale(1+3*(inot-1)+1)
                                        if (ndim .eq. 3) then
                                            axer(3) = axer(3)+vale(1+3*(inot-1)+2)
                                        end if
                                    end do
                                    axer(1) = axer(1)/nbpt
                                    axer(2) = axer(2)/nbpt
                                    axer(3) = axer(3)/nbpt
                                    goto 17
                                end if
                            end do
17                          continue
                        end do
!
!                LE NOEUD SUR L'AXE N'APPARTIENT A AUCUNE MAILLE
                        if (ipt2 .eq. 0) then
                            do ii = 1, nbcmp
                                zl(jcnsl-1+(inoe-1)*nbcmp+ii) = .false.
                            end do
                            goto 16
                        end if
!
                        axer(1) = axer(1)-orig(1)
                        axer(2) = axer(2)-orig(2)
                        axer(3) = axer(3)-orig(3)
                        b_n = to_blas_int(3)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        prosca = ddot(b_n, axer, b_incx, axez, b_incy)
                        axer(1) = axer(1)-prosca*axez(1)
                        axer(2) = axer(2)-prosca*axez(2)
                        if (ndim .eq. 3) then
                            axer(3) = axer(3)-prosca*axez(3)
                        else
                            axer(3) = 0.0d0
                        end if
                        xnormr = 0.0d0
                        call normev(axer, xnormr)
                        if (xnormr .lt. epsi) then
                            k8b = int_to_char8(inoe)
                            valk = k8b
                            call utmess('F', 'ALGORITH14_91', sk=valk)
                        end if
                    end if
                    call provec(axez, axer, axet)
                    xnormr = 0.0d0
                    call normev(axet, xnormr)
                    do i = 1, 3
                        pgl(1, i) = axer(i)
                        pgl(2, i) = axez(i)
                        pgl(3, i) = axet(i)
                    end do
                    if (.not. zl(jcnsl-1+(inoe-1)*nbcmp+1)) goto 16
                    if (tsca .eq. 'R') then
! CHAMP REEL
                        do ii = 1, nbcmp
                            valet(ii) = zr(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        valed(1) = valet(1)
                        valed(2) = valet(4)
                        valed(3) = valet(2)
                        valed(4) = valet(5)
                        valed(5) = valet(6)
                        valed(6) = valet(3)
                        call utpsgl(1, 3, pgl, valed, valet)
                        valer(1) = valet(1)
                        valer(2) = valet(3)
                        valer(3) = valet(6)
                        valer(4) = valet(2)
                        valer(5) = valet(4)
                        valer(6) = valet(5)
                        do ii = 1, nbcmp
                            zr(jcnsv-1+(inoe-1)*nbcmp+ii) = valer(licmpu(ii))
                        end do
                    else
! CHAMP COMPLEXE
                        do ii = 1, nbcmp
                            valetc(ii) = zc(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        valed(1) = dble(valetc(1))
                        valed(2) = dble(valetc(4))
                        valed(3) = dble(valetc(2))
                        valed(4) = dble(valetc(5))
                        valed(5) = dble(valetc(6))
                        valed(6) = dble(valetc(3))
                        call utpsgl(1, 3, pgl, valed, valet)
                        valer(1) = valet(1)
                        valer(2) = valet(3)
                        valer(3) = valet(6)
                        valer(4) = valet(2)
                        valer(5) = valet(4)
                        valer(6) = valet(5)
!
                        valed(1) = dimag(valetc(1))
                        valed(2) = dimag(valetc(4))
                        valed(3) = dimag(valetc(2))
                        valed(4) = dimag(valetc(5))
                        valed(5) = dimag(valetc(6))
                        valed(6) = dimag(valetc(3))
                        call utpsgl(1, 3, pgl, valed, valet)
                        valei(1) = valet(1)
                        valei(2) = valet(3)
                        valei(3) = valet(6)
                        valei(4) = valet(2)
                        valei(5) = valet(4)
                        valei(6) = valet(5)
                        do ii = 1, nbcmp
                            zc(jcnsv-1+(inoe-1)*nbcmp+ii) = dcmplx(valer(ii), valei(ii))
                        end do
!
                    end if
16                  continue
                end do
            else
! VECTEUR
!
                do ino = 1, nbnoeu
                    if (nbn .ne. 0) then
                        inoe = zi(idnoeu+ino-1)
                    else
                        inoe = ino
                    end if
                    axer(1) = vale(1+3*(inoe-1))-orig(1)
                    axer(2) = vale(1+3*(inoe-1)+1)-orig(2)
                    if (ndim .eq. 3) then
                        axer(3) = vale(1+3*(inoe-1)+2)-orig(3)
                    else
                        axer(3) = 0.0d0
                    end if
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    prosca = ddot(b_n, axer, b_incx, axez, b_incy)
                    axer(1) = axer(1)-prosca*axez(1)
                    axer(2) = axer(2)-prosca*axez(2)
                    if (ndim .eq. 3) then
                        axer(3) = axer(3)-prosca*axez(3)
                    else
                        axer(3) = 0.0d0
                    end if
                    xnormr = 0.0d0
                    call normev(axer, xnormr)
                    if (xnormr .lt. epsi) then
                        k8b = int_to_char8(inoe)
                        call utmess('A', 'ALGORITH2_13')
                        call jeveuo(ma//'.CONNEX', 'L', jconx1)
                        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', jconx2)
                        ipt2 = 0
                        do inel = 1, nbma
                            nbpt = zi(jconx2-1+inel+1)-zi(jconx2-1+inel)
                            do ipt = 1, nbpt
                                inot = zi(jconx1-1+zi(jconx2-1+inel)+ipt-1)
                                if (inot .eq. inoe) then
                                    axer(1) = 0.0d0
                                    axer(2) = 0.0d0
                                    axer(3) = 0.0d0
                                    do ipt2 = 1, nbpt
                                        inot = zi(jconx1-1+zi(jconx2-1+inel)+ipt2-1)
                                        axer(1) = axer(1)+vale(1+3*(inot-1))
                                        axer(2) = axer(2)+vale(1+3*(inot-1)+1)
                                        if (ndim .eq. 3) then
                                            axer(3) = axer(3)+vale(1+3*(inot-1)+2)
                                        end if
                                    end do
                                    axer(1) = axer(1)/nbpt
                                    axer(2) = axer(2)/nbpt
                                    axer(3) = axer(3)/nbpt
                                    goto 24
                                end if
                            end do
24                          continue
                        end do
!                LE NOEUD SUR L'AXE N'APPARTIENT A AUCUNE MAILLE
                        if (ipt2 .eq. 0) then
                            do ii = 1, nbcmp
                                zl(jcnsl-1+(inoe-1)*nbcmp+ii) = .false.
                            end do
                            goto 23
                        end if
                        axer(1) = axer(1)-orig(1)
                        axer(2) = axer(2)-orig(2)
                        axer(3) = axer(3)-orig(3)
                        b_n = to_blas_int(3)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        prosca = ddot(b_n, axer, b_incx, axez, b_incy)
                        axer(1) = axer(1)-prosca*axez(1)
                        axer(2) = axer(2)-prosca*axez(2)
                        if (ndim .eq. 3) then
                            axer(3) = axer(3)-prosca*axez(3)
                        else
                            axer(3) = 0.0d0
                        end if
                        xnormr = 0.0d0
                        call normev(axer, xnormr)
                        if (xnormr .lt. epsi) then
                            k8b = int_to_char8(inoe)
                            valk = k8b
                            call utmess('F', 'ALGORITH14_91', sk=valk)
                        end if
                    end if
                    call provec(axez, axer, axet)
                    do i = 1, 3
                        pgl(1, i) = axer(i)
                        pgl(2, i) = axez(i)
                        pgl(3, i) = axet(i)
                    end do
                    if (.not. zl(jcnsl-1+(inoe-1)*nbcmp+1)) goto 23
                    if (tsca .eq. 'R') then
! VECTEUR REEL
                        do ii = 1, nbcmp
                            valed(ii) = zr(jcnsv-1+(inoe-1)*nbcmp+ii)
                        end do
                        if (ndim .eq. 3) then
                            call utpvgl(1, 3, pgl, valed, valer)
                            if (nbcmp .gt. 3) call utpvgl(1, 3, pgl, valed(4), valer(4))
                        else
                            call ut2vgl(1, 2, pgl, valed, valer)
                        end if
                        do ii = 1, nbcmp
                            zr(jcnsv-1+(inoe-1)*nbcmp+ii) = valer(ii)
                        end do
                    else
! VECTEUR COMPLEXE
                        do ii = 1, nbcmp
                            valetc(ii) = zc(jcnsv-1+(inoe-1)*nbcmp+ii)
                            valed(ii) = dble(valetc(ii))
                            valet(ii) = dimag(valetc(ii))
                        end do
                        if (ndim .eq. 3) then
                            call utpvgl(1, 3, pgl, valed, valer)
                            call utpvgl(1, 3, pgl, valet, valei)
                            if (nbcmp .gt. 3) then
                                call utpvgl(1, 3, pgl, valed(4), valer(4))
                                call utpvgl(1, 3, pgl, valet(4), valei(4))
                            end if
                        else
                            call ut2vgl(1, 2, pgl, valed, valer)
                            call ut2vgl(1, 2, pgl, valet, valei)
                        end if
                        do ii = 1, nbcmp
                            zc(jcnsv-1+(inoe-1)*nbcmp+ii) = dcmplx(valer(ii), valei(ii))
                        end do
!
                    end if
!
23                  continue
                end do
            end if
        end if
!
! Fin de la boucle sur les occcurrences du mot-clé AFFE
        call jedetr(mesnoe)
    end do
    call cnscno(chams1, ' ', 'NON', 'G', champ1, &
                'F', ibid)
    call detrsd('CHAM_NO_S', chams1)
!
    call jedema()
!
end subroutine
