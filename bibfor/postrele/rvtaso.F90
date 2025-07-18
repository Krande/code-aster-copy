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
subroutine rvtaso(releve, nomcmp, nbcmp, nbco, nbsp, &
                  nomtab, iocc, ncheff, i1)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
!
    integer(kind=8) :: iocc, i1, nbcmp, nbco, nbsp
    real(kind=8) :: releve(*)
    character(len=8) :: nomcmp(*)
    character(len=16) :: ncheff
    character(len=19) :: nomtab
!     MISE EN TABLEAU POUR UNE SOMME
!     ------------------------------------------------------------------
! IN  : RELEVE : TABLE DES RELEVE DE VALEURS
! IN  : NOMCMP : NOM DES COMPOSANTES
! IN  : NBCMP  : NOMBRE DE NOMCMP
! IN  : NBCO   : NOMBRE DE COUCHES PAR POINT
! IN  : NBSP   : NOMBRE DE SOUS-PT PAR POINT
! IN  : NOMTAB : INTITULE DE LA TABLE
!     ------------------------------------------------------------------
    integer(kind=8) :: nbpar, ilign, ls, lc, isp, icp, ico, n1, adrval, adracc, jacc, ik, ir, ii
    integer(kind=8) :: valei(12), nbacc, nbpr, jaces, iac, iadr, iord(1)
    real(kind=8) :: prec, valer(10)
    complex(kind=8) :: c16b
    aster_logical :: exist
    character(len=3) :: typpar
    character(len=8) :: acces, nomres, ctype, crit, k8b
    character(len=16) :: intitu
    character(len=24) :: nomval, nomacc, nnores, nopara(18), nomjv
    character(len=80) :: valek(11)
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call getvtx('ACTION', 'INTITULE', iocc=iocc, scal=intitu, nbret=n1)
!
    call getvr8('ACTION', 'PRECISION', iocc=iocc, scal=prec, nbret=n1)
    call getvtx('ACTION', 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
!
    nomval = ncheff//'.VALACCE'
    nomacc = ncheff//'.TYPACCE'
    nnores = ncheff//'.NOMRESU'
    call jeveuo(nomacc, 'L', jacc)
!
    ik = 1
    ii = 0
    ir = 0
    nbpar = 1
    nopara(nbpar) = 'INTITULE'
    valek(ik) = intitu
!
    if (zk8(jacc) .eq. 'DIRECT  ') then
        call jeveuo(jexnum(ncheff//'.LSCHEFF', 1), 'L', jacc)
        nbpar = nbpar+1
        nopara(nbpar) = 'CHAM_GD'
        ik = ik+1
        valek(ik) = zk24(jacc) (1:8)
    else
        call jeveuo(nnores, 'L', jacc)
        nomres = zk16(jacc) (1:8)
        nbpar = nbpar+1
        nopara(nbpar) = 'RESU'
        ik = ik+1
        valek(ik) = nomres
        nbpar = nbpar+1
        nopara(nbpar) = 'NOM_CHAM'
        ik = ik+1
        valek(ik) = zk16(jacc+1)
        call jeveuo(nomacc, 'L', adracc)
        call jeveuo(nomval, 'L', adrval)
        acces = zk8(adracc)
        if (acces(1:1) .eq. 'O') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_ORDRE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
            nomjv = '&&RVRCCM.NOMS_ACCES'
            call rsnopa(nomres, 0, nomjv, nbacc, nbpr)
            if (nbacc .ne. 0) then
                call jeveuo(nomjv, 'L', jaces)
                do iac = 1, nbacc
                    call rsadpa(nomres, 'L', 1, zk16(jaces-1+iac), zi(adrval+i1-1), &
                                1, sjv=iadr, styp=ctype, istop=0)
                    call tbexip(nomtab, zk16(jaces-1+iac), exist, typpar)
                    if (.not. exist) then
                        call tbajpa(nomtab, 1, zk16(jaces-1+iac), ctype)
                    end if
                    nbpar = nbpar+1
                    nopara(nbpar) = zk16(jaces-1+iac)
                    if (ctype(1:1) .eq. 'I') then
                        ii = ii+1
                        valei(ii) = zi(iadr)
                    else if (ctype(1:1) .eq. 'R') then
                        ir = ir+1
                        valer(ir) = zr(iadr)
                    else if (ctype(1:3) .eq. 'K80') then
                        ik = ik+1
                        valek(ik) = zk80(iadr)
                    else if (ctype(1:3) .eq. 'K32') then
                        ik = ik+1
                        valek(ik) = zk32(iadr)
                    else if (ctype(1:3) .eq. 'K24') then
                        ik = ik+1
                        valek(ik) = zk24(iadr)
                    else if (ctype(1:3) .eq. 'K16') then
                        ik = ik+1
                        valek(ik) = zk16(iadr)
                    else if (ctype(1:2) .eq. 'K8') then
                        ik = ik+1
                        valek(ik) = zk8(iadr)
                    end if
                end do
                call jedetr(nomjv)
            end if
        else if (acces(1:1) .eq. 'M') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'NUME_MODE', zi(adrval+i1-1), 0.d0, k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_MODE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
        else if (acces(1:1) .eq. 'I') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'INST', 0, zr(adrval+i1-1), k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            nopara(nbpar) = 'INST'
            ir = ir+1
            valer(ir) = zr(adrval+i1-1)
        else if (acces(1:1) .eq. 'F') then
            nbpar = nbpar+1
            nopara(nbpar) = 'FREQ'
            ir = ir+1
            valer(ir) = zr(adrval+i1-1)
        end if
    end if
!
    if (nbco .gt. 1) then
        call tbexip(nomtab, 'NUME_COUCHE', exist, typpar)
        if (.not. exist) then
            call tbajpa(nomtab, 1, 'NUME_COUCHE', 'I')
        end if
        nbpar = nbpar+1
        nopara(nbpar) = 'NUME_COUCHE'
    end if
    if (nbsp .gt. 1) then
        call tbexip(nomtab, 'NUME_GAUSS', exist, typpar)
        if (.not. exist) then
            call tbajpa(nomtab, 1, 'NUME_GAUSS', 'I')
        end if
        nbpar = nbpar+1
        nopara(nbpar) = 'NUME_GAUSS'
    end if
!
    do icp = 1, nbcmp, 1
        nbpar = nbpar+1
        nopara(nbpar) = nomcmp(icp)
    end do
!
    ASSERT(nbpar .le. 15)
    ASSERT(ii+2 .le. 10)
    ASSERT(ir+nbcmp .le. 10)
    ASSERT(ik .le. 10)
!
    ls = nbcmp
    lc = nbsp*ls
    ilign = 0
!
    do ico = 1, nbco, 1
        valei(ii+1) = ico
!
        do isp = 1, nbsp, 1
            valei(ii+2) = isp
!
            do icp = 1, nbcmp, 1
                valer(ir+icp) = releve(icp+lc*(ico-1)+ls*(isp-1))
            end do
!
            call tbajli(nomtab, nbpar, nopara, valei, valer, &
                        [c16b], valek, ilign)
!
        end do
!
    end do
!
    call jedema()
end subroutine
