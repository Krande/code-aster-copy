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
subroutine crprol()
    implicit none
!
!     COMMANDE:  PROL_RTZ
!
!
!
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/normev.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbextb.h"
#include "asterfort/tbexve.h"
#include "asterfort/tbtri.h"
#include "asterfort/tbutnu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: ibid, ndimf, nbnoi, nbnof, nbinst, iad
    integer(kind=8) :: vali, iret
    integer(kind=8) :: jinst, iord, jcnsvl, jcnsle, nbval
    integer(kind=8) :: axyzmf, jtbcor, jtbres
    integer(kind=8) :: imin, imax, inoi, inof, indice, jtbpdg
    integer(kind=8) :: jtbnoe, ordef, ino, inomin, inomax
    real(kind=8) :: xnormr, prec, rmin, rmax
    real(kind=8) :: rval, max, min, prosca, rref, rpro
    real(kind=8) :: lambda, orig(3), axer(3), axet(3), axez(3)
    real(kind=8) :: dinst
    complex(kind=8) :: cbid
    character(len=8) :: pdroit, pgauch, nommaf, k8b, nomgd
    character(len=8) :: nom1, crit, licmpr, table, resu
    character(len=16) :: motfac, typres, nomcmd
    character(len=19) :: tabcor, tabval, tabpdg, tabres, tabnoe
    character(len=19) :: cnoinr, cnsinr
    character(len=24) :: knum, tabl2
    character(len=24) :: valk(2)
    real(kind=8), pointer :: cnsv(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
    cbid = (0.d0, 0.d0)
    motfac = 'PROL_RTZ'
    call getres(resu, typres, nomcmd)
    if (typres .ne. 'EVOL_THER') then
        call utmess('F', 'ALGORITH2_44')
    end if
!
! --- RECUPERATION DES DONNEES UTILISATEUR :
!     ------------------------------------
!
    call getvid(motfac, 'MAILLAGE_FINAL', iocc=1, scal=nommaf, nbret=ibid)
    call getvid(motfac, 'TABLE', iocc=1, scal=table, nbret=ibid)
    call getvtx(motfac, 'PROL_DROITE', iocc=1, scal=pdroit, nbret=ibid)
    call getvtx(motfac, 'PROL_GAUCHE', iocc=1, scal=pgauch, nbret=ibid)
    call getvr8(motfac, 'ORIGINE', iocc=1, nbval=3, vect=orig, &
                nbret=ibid)
    call getvr8(motfac, 'AXE_Z', iocc=1, nbval=3, vect=axez, &
                nbret=ibid)
!
    call dismoi('NB_NO_MAILLA', nommaf, 'MAILLAGE', repi=nbnof)
    call dismoi('Z_CST', nommaf, 'MAILLAGE', repk=k8b)
    ndimf = 3
    if (k8b .eq. 'OUI') ndimf = 2
    if (ndimf .ne. 3) then
        valk(1) = nommaf
        call utmess('F', 'ALGORITH12_68', sk=valk(1))
    end if
    call jeveuo(nommaf//'.COORDO    .VALE', 'L', axyzmf)
!
    call tbexp2(table, 'INST')
    call tbexp2(table, 'COOR_X')
    call tbexp2(table, 'TEMP')
    call normev(axez, xnormr)
!
    knum = '&&RS1D3D.INSTANT'
    call tbutnu(motfac, 1, knum, nbinst, table, &
                prec, crit)
    call jeveuo(knum, 'L', jinst)
!
    call rscrsd('G', resu, typres, nbinst)
!
    tabcor = '&&RS1D3D.COORMA'
    tabval = '&&RS1D3D.VALTEM'
    tabres = '&&RS1D3D.COORES'
    tabpdg = '&&RS1D3D.TABPDG'
    tabnoe = '&&RS1D3D.TABNOE'
!
    call wkvect(tabpdg, 'V V R', nbnof, jtbpdg)
    call wkvect(tabnoe, 'V V I', nbnof, jtbnoe)
!
    tabl2 = '&&RS1D3D.TABL2'
    cnsinr = '&&RS1D3D.CNSINR'
!
    do iord = 1, nbinst
        call jemarq()
        call jerecu('V')
!
! ------ ON EXTRAIT LA SOUS-TABLE POUR L'INSTANT COURANT
!
        dinst = zr(jinst+iord-1)
        call tbextb(table, 'V', tabl2, 1, 'INST', &
                    'EQ', [ibid], [dinst], [cbid], k8b, &
                    [prec], crit, iret)
        if (iret .eq. 10) then
            valk(1) = 'INST'
            valk(2) = table
            call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
        else if (iret .eq. 20) then
            valk(1) = table
            valk(2) = 'INST'
            call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
        end if
!
! ------ ON RECUPERE LES COORCONNEES DES NOEUDS POUR L'INSTANT COURANT
!
        call tbexve(tabl2, 'COOR_X', tabcor, 'V', nbnoi, &
                    k8b)
        call jeveuo(tabcor, 'L', jtbcor)
        call wkvect(tabres, 'V V I', nbnoi, jtbres)
!
! ------ ON RECUPERE LES VALEURS DE TEMPERATURE AUX NOEUDS
!
        call tbexve(tabl2, 'TEMP', tabval, 'V', nbval, &
                    k8b)
        call jeveuo(tabval, 'L', jcnsvl)
!
! ------ ON VERIFIE QUE LE MAILLAGE 1D COMMENCE A 0.
!
        rmax = 0.0d0
        rmin = 1.d0/r8prem()
        do inoi = 1, nbnoi
            rval = zr(jtbcor-1+inoi)
            rmax = max(rval, rmax)
            rmin = min(rval, rmin)
            if (rmin .ne. 0.0d0) then
                call utmess('F', 'ALGORITH12_69')
            end if
        end do
!
! ------ ON TRIE PAR ORDRE CROISSANT
!
        call tbtri(nbnoi, zi(jtbres), tabchr=zr(jtbcor))
!
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        prosca = ddot(b_n, orig, b_incx, axez, b_incy)
        axer(1) = orig(1)-prosca*axez(1)
        axer(2) = orig(2)-prosca*axez(2)
        axer(3) = orig(3)-prosca*axez(3)
        call normev(axer, xnormr)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        rref = ddot(b_n, orig, b_incx, axer, b_incy)
        rref = abs(rref)
!
        nomgd = 'TEMP_R'
        licmpr = 'TEMP'
        call cnscre(nommaf, nomgd, 1, licmpr, 'V', &
                    cnsinr)
        call jeveuo(cnsinr//'.CNSV', 'E', vr=cnsv)
        call jeveuo(cnsinr//'.CNSL', 'E', jcnsle)
!
        indice = 0
        do inof = 1, nbnof
            zl(jcnsle-1+(inof-1)+1) = .false.
            axet(1) = zr(axyzmf+3*(inof-1)-1+1)
            axet(2) = zr(axyzmf+3*(inof-1)-1+2)
            axet(3) = zr(axyzmf+3*(inof-1)-1+3)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            prosca = ddot(b_n, axet, b_incx, axez, b_incy)
            axer(1) = axet(1)-prosca*axez(1)
            axer(2) = axet(2)-prosca*axez(2)
            axer(3) = axet(3)-prosca*axez(3)
            call normev(axer, xnormr)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            rpro = ddot(b_n, axet, b_incx, axer, b_incy)
            rval = rpro-rref
            if (rval .lt. 0.0d0) then
                indice = indice+1
                zr(jtbpdg-1+indice) = rval
                zi(jtbnoe-1+indice) = inof
                goto 3
            end if
            do inoi = 1, nbnoi
                if (rval .le. zr(jtbcor-1+zi(jtbres-1+inoi))) then
                    imin = zi(jtbres-1+inoi-1)
                    rmin = zr(jtbcor-1+imin)
                    imax = zi(jtbres-1+inoi)
                    rmax = zr(jtbcor-1+imax)
                    goto 5
                end if
            end do
            indice = indice+1
            zr(jtbpdg-1+indice) = rval
            zi(jtbnoe-1+indice) = inof
            goto 3
5           continue
            if ((rmax-rmin) .eq. 0.0d0) then
                call utmess('F', 'ALGORITH12_70')
            end if
            lambda = (rval-rmin)/(rmax-rmin)
            cnsv((inof-1)+1) = (1-lambda)*zr(jcnsvl-1+(imin-1)+1)+lambda*zr(jcnsvl-1+(imax-1)+1)
            zl(jcnsle-1+(inof-1)+1) = .true.
3           continue
        end do
        do ordef = 1, indice
            ino = zi(jtbnoe-1+ordef)
            rval = zr(jtbpdg-1+ordef)
            if (rval .lt. 0.0d0) then
                if (pgauch .eq. 'EXCLU') then
                    nom1 = int_to_char8(ino)
                    valk(1) = nom1
                    call utmess('F', 'ALGORITH12_71', sk=valk(1))
                else if (pgauch .eq. 'CONSTANT') then
                    inomin = zi(jtbres)
                    cnsv((ino-1)+1) = zr(jcnsvl-1+(inomin-1)+1)
                    zl(jcnsle-1+(ino-1)+1) = .true.
                else
                    inomin = zi(jtbres)
                    inomax = zi(jtbres+1)
                    rmin = zr(jtbcor-1+inomin)
                    rmax = zr(jtbcor-1+inomax)
                    lambda = (rmin-rval)/(rmax-rval)
                    cnsv((ino-1)+1) = ( &
                                      zr(jcnsvl-1+(inomin-1)+1)-zr(jcnsvl-1+(inomax-1)+1)*lambda &
                                      )/(1-lambda &
                                         )
                    zl(jcnsle-1+(ino-1)+1) = .true.
                end if
            else
                if (pdroit .eq. 'EXCLU') then
                    nom1 = int_to_char8(ino)
                    valk(1) = nom1
                    call utmess('F', 'ALGORITH12_72', sk=valk(1))
                else if (pdroit .eq. 'CONSTANT') then
                    inomax = zi(jtbres-1+nbnoi)
                    cnsv((ino-1)+1) = zr(jcnsvl-1+(inomax-1)+1)
                    zl(jcnsle-1+(ino-1)+1) = .true.
                else
                    inomin = zi(jtbres-1+nbnoi-1)
                    inomax = zi(jtbres-1+nbnoi)
                    rmin = zr(jtbcor-1+inomin)
                    rmax = zr(jtbcor-1+inomax)
                    lambda = (rmax-rmin)/(rval-rmin)
                    cnsv((ino-1)+1) = ( &
                                      zr( &
                                      jcnsvl-1+(inomax-1)+1)-zr(jcnsvl-1+(inomin-1)+1)*(1-lambda) &
                                      )/lambda
                    zl(jcnsle-1+(ino-1)+1) = .true.
                end if
            end if
        end do
!
        call rsexch(' ', resu, 'TEMP', iord, cnoinr, &
                    ibid)
        if (ibid .ne. 100) then
            valk(1) = resu
            vali = iord
            call utmess('F', 'ALGORITH12_73', sk=valk(1), si=vali)
        end if
        call cnscno(cnsinr, ' ', 'NON', 'G', cnoinr, &
                    'F', ibid)
        call rsnoch(resu, 'TEMP', iord)
        call detrsd('CHAM_NO_S', cnsinr)
!
        call rsadpa(resu, 'E', 1, 'INST', iord, &
                    0, sjv=iad, styp=k8b)
        zr(iad) = dinst
!
        call detrsd('TABLE', tabl2)
        call jedetr(tabcor)
        call jedetr(tabres)
        call jedetr(tabval)
!
        call jedema()
    end do
    call jedetr(knum)
    call jedetr(tabpdg)
    call jedetr(tabnoe)
! FIN ------------------------------------------------------------------
    call jedema()
end subroutine
