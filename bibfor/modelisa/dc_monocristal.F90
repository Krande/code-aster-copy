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
subroutine dc_monocristal(nboccm, sdcomp)
!
! person_in_charge: jean-luc.flejou at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbexlr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
!
    integer(kind=8), intent(out) :: nboccm
    character(len=8), intent(in) :: sdcomp
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_COMPOR
!
! MONOCRISTAL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_para
    parameter(nb_para=5)
!
    character(len=8) :: materi, typpar(nb_para), chaine
    character(len=16) :: nompar(nb_para), ecoule, ecrois, ecroci, elasti
    character(len=24) :: nomvar(200)
    character(len=16) :: fasygl, noms(6), comdes, rota, tbinte, systgl
    character(len=19) :: listr
    real(kind=8) :: ms(6), ng(3), q(3, 3), lg(3), pgl(3, 3)
    complex(kind=8) :: cbid
    integer(kind=8) :: ifm, niv, nbtbsg, nums(2), indvar
    integer(kind=8) :: iocc, nbmat, nbecou, nbecro, nbcine, nbelas, nbfasy
    integer(kind=8) :: i, j, nbela1, nbsys, nvi, imk, imi, ipr, itab, itsg, irra, irr2
    integer(kind=8) :: ncprr, ir, irota, iadlr, decal, nbrota, nbsyst, tabdes(13)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    cbid = (0.d0, 0.d0)
    comdes = '&&OP0059.TABLETX'
    call tbcrsd(comdes, 'V')
    nompar(1) = 'FAMI_SYST_GLIS'
    nompar(2) = 'MAT_SYST'
    nompar(3) = 'ECOULEMENT'
    nompar(4) = 'ECRO_ISOT'
    nompar(5) = 'ECRO_CINE'
    typpar(1) = 'K16'
    typpar(2) = 'K16'
    typpar(3) = 'K16'
    typpar(4) = 'K16'
    typpar(5) = 'K16'
    nbsyst = 0
    nbelas = 0
    nvi = 6
!     DEFORMATION PLASTIQUE CUMULEE MACROSCOPIQUE EQUIVALENTE
    nvi = nvi+1
    call tbajpa(comdes, nb_para, nompar, typpar)
    call getfac('MONOCRISTAL', nboccm)
    call wkvect(sdcomp//'.CPRK', 'G V K24', 5*nboccm+1, imk)
!     DIMENSION MAX DE CPRR : 1812 = 12+5*6*30+30*30
!     ORGANISATION DE CPRR :
!     1 : NB TABLES SYST GLIS
!     2 : POSITION DE LA TABLE D'INTERACTION. 0 SINON
!     3 : NOMBRE DE SYSTEMES TABLE  1
!     4 : POSITION DE LA TABLE  1,
!     5 : NOMBRE DE SYSTEMES TABLE 2
!     6 : POSITION DE LA TABLE  2
!     ...
!     CPRR(4)   : TABLE DE SYS GLIS1 (6*NSYS1 VALEURS)
!     CPRR(6)    : TABLE DE SYS GLIS2 (6*NSYS2 VALEURS)
!     ...
!     CPRR(2)    : TABLE D'INTERACTION (NBSYST*NBSYST VALEURS)
    ncprr = 1812
    call wkvect(sdcomp//'.CPRR', 'G V R', ncprr, ipr)
    nbtbsg = 0
    decal = 12
    do i = 1, 13
        tabdes(i) = 0
    end do
    irra = 0
    irr2 = 0
!
    do iocc = 1, nboccm
        call getvid('MONOCRISTAL', 'MATER', iocc=iocc, scal=materi, nbret=nbmat)
        call getvtx('MONOCRISTAL', 'ECOULEMENT', iocc=iocc, scal=ecoule, nbret=nbecou)
        call getvtx('MONOCRISTAL', 'ECRO_ISOT', iocc=iocc, scal=ecrois, nbret=nbecro)
        call getvtx('MONOCRISTAL', 'ECRO_CINE', iocc=iocc, scal=ecroci, nbret=nbcine)
        call getvtx('MONOCRISTAL', 'ELAS', iocc=iocc, scal=elasti, nbret=nbela1)
        if (nbela1 .gt. 0) then
            if (nbelas .eq. 0) then
                nbelas = 1
            else
                call utmess('F', 'MODELISA5_64')
            end if
        end if
!        CAS DES LOIS DD
        if (ecoule(1:7) .eq. 'MONO_DD') then
            ecrois = ecoule
            ecroci = ' '
        end if
!        CAS DES LOIS DD_CC_IRRA
        if (ecoule .eq. 'MONO_DD_CC_IRRA') then
            irra = irra+1
        end if
!        CAS DES LOIS DD_CFC_IRRA
        if (ecoule .eq. 'MONO_DD_CFC_IRRA') then
            irr2 = irr2+1
        end if
!
        call getvtx('MONOCRISTAL', 'FAMI_SYST_GLIS', iocc=iocc, scal=fasygl, nbret=nbfasy)
        noms(1) = fasygl
        noms(2) = materi
        noms(3) = ecoule
        noms(4) = ecrois
        noms(5) = ecroci
        if (fasygl .eq. 'UTILISATEUR') then
            call getvid('MONOCRISTAL', 'TABL_SYST_GLIS', iocc=iocc, scal=systgl, nbret=itsg)
            noms(1) = 'UTIL'
            nbtbsg = nbtbsg+1
            call codent(nbtbsg, 'G', noms(1) (5:5))
            noms(1) (6:8) = '___'
            noms(1) (9:16) = systgl(1:8)
            fasygl = noms(1)
            listr = '&&LCMMAT.TABL_SYSGL'
            call tbexlr(systgl, listr, 'V')
            call jeveuo(listr//'.VALE', 'L', iadlr)
            nbsys = nint(zr(iadlr+2))
!           VERIF QUE LA MATRICE EST CARREE
            if (6 .ne. zr(iadlr+1)) then
                call utmess('F', 'COMPOR2_19', sr=zr(iadlr+1))
            end if
            zr(ipr+2+2*(nbtbsg-1)) = nbsys
            zr(ipr+2+2*(nbtbsg-1)+1) = decal+1
            b_n = to_blas_int(6*nbsys)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(iadlr+3), b_incx, zr(ipr+decal), b_incy)
            tabdes(8+iocc) = nbsys
            call detrsd('LISTR8', listr)
!
            if (niv .eq. 2) then
                write (ifm, *) ' TABLE SYSTEMES DE GLISSEMENT FAMILLE', &
                    iocc
                write (ifm, *) ' NX     NY     NZ     MX     MY     MZ '
                do i = 1, nbsys
                    write (ifm, '(I2,6(1X,E11.4))') i, (zr(ipr-1+decal+6* &
                                                           (i-1)+j), j=1, 6)
                end do
            end if
!
            decal = decal+6*nbsys
!
        else
            ir = 0
            call lcmmsg(fasygl, nbsys, 0, pgl, ms, &
                        ng, lg, ir, q)
        end if
!
        call tbajli(comdes, 5, nompar, [0], [0.d0], &
                    [cbid], noms, 0)
        do j = 1, 5
            zk24(imk-1+(iocc-1)*5+j) = noms(j)
        end do
        ir = 0
        nvi = nvi+4*nbsys
        nbsyst = nbsyst+nbsys
!
    end do
!
    zr(ipr) = nbtbsg
    zr(ipr+1) = decal+1
!     INDICATEUR PLASTIQUE
    nvi = nvi+1
!     CONTRAINTE DE CLIVAGE MAX
    nvi = nvi+1
!     ROTATION DE RESEAU
    call getvtx(' ', 'ROTA_RESEAU', scal=rota, nbret=nbrota)
    irota = 0
    if (nbrota .ne. 0) then
        if (rota .ne. 'NON') then
            if (rota .eq. 'POST') irota = 1
            if (rota .eq. 'CALC') irota = 2
        end if
        if (irota .gt. 0) nvi = nvi+16
    end if
!     RHO_IRR
    if (irra .gt. 0) then
        nvi = nvi+12*nboccm
    end if
!     RHO_IRR
    if (irr2 .gt. 0) then
        nvi = nvi+24*nboccm
    end if
!
    noms(2) = ecoule
    nums(1) = nboccm
    nums(2) = nvi
    call utmess('I', 'COMPOR4_1')
    call utmess('I', 'COMPOR4_2', si=nboccm)
    call utmess('I', 'COMPOR4_5', sk='MONOCRISTAL')
    call utmess('I', 'COMPOR4_6', sk=ecoule)
    call utmess('I', 'COMPOR4_9', si=nvi)
!
    nomvar(1) = 'EPSPXX'
    nomvar(2) = 'EPSPYY'
    nomvar(3) = 'EPSPZZ'
    nomvar(4) = 'EPSPXY'
    nomvar(5) = 'EPSPXZ'
    nomvar(6) = 'EPSPYZ'
    do i = 1, nbsyst
        call codent(i, 'G', chaine)
        nomvar(6+3*i-2) = 'ALPHA'//chaine
        nomvar(6+3*i-1) = 'GAMMA'//chaine
        nomvar(6+3*i) = 'P'//chaine
    end do
    if (irra .gt. 0) then
        do i = 1, 12*nboccm
            call codent(i, 'G', chaine)
            nomvar(6+3*nbsyst+i) = 'RHO_IRRA_'//chaine
        end do
    end if
!
    if (irr2 .gt. 0) then
        do i = 1, 12*nboccm
            call codent(i, 'G', chaine)
            nomvar(6+3*nbsyst+i) = 'RHO_LOOPS_'//chaine
            call codent(i, 'G', chaine)
            nomvar(6+4*nbsyst+i) = 'PHI_VOIDS_'//chaine
        end do
    end if
!
    indvar = 6+3*nbsyst
    if (irra .gt. 0) indvar = indvar+12*nboccm
!
    if (irr2 .gt. 0) indvar = indvar+24*nboccm
!
    if (irota .gt. 0) then
        do i = 1, 16
            call codent(i, 'G', chaine)
            nomvar(indvar+i) = 'ROTA_'//chaine
        end do
        indvar = indvar+16
    end if
!
    do i = 1, nbsyst
        call codent(i, 'G', chaine)
        nomvar(indvar+i) = 'TAU_'//chaine
    end do
!
    nomvar(nvi-2) = 'SIGM_CLIV'
    nomvar(nvi-1) = 'EPSPEQ'
    nomvar(nvi) = 'NBITER'
!
    do i = 1, nvi
        call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=i)
    end do
!
!
    zk24(imk+5*nboccm) = elasti
    call getvid(' ', 'MATR_INTER', scal=tbinte, nbret=itab)
    if (itab .ne. 0) then
        listr = '&&LCMMAT.TABL_INTER'
        call tbexlr(tbinte, listr, 'V')
        call jeveuo(listr//'.VALE', 'L', iadlr)
!        VERIF QUE LA MATRICE EST CARREE
        if (zr(iadlr+1) .ne. zr(iadlr+2)) then
            call utmess('F', 'COMPOR2_15', nr=2, valr=zr(iadlr+1))
        end if
!        VERIF QUE LE NB DE SYST EST OK
        if (zr(iadlr+1) .ne. nbsyst) then
            call utmess('F', 'COMPOR2_17', si=nbsyst)
        end if
        b_n = to_blas_int(nbsyst*nbsyst)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(iadlr+3), b_incx, zr(ipr+decal), b_incy)
!        VERIF QUE LA MATRICE EST SYMETRIQUE
        do i = 1, nbsyst
            do j = 1, nbsyst
                if (zr(ipr-1+decal+nbsyst*(i-1)+j) .ne. zr(ipr-1+decal+nbsyst*(j-1)+i)) then
                    call utmess('F', 'COMPOR2_18')
                end if
            end do
        end do
        call jedetc('V', listr, 1)
!
        if (niv .eq. 2) then
            write (ifm, *) ' MATRICE INTERACTION UTILISATEUR'
            do i = 1, nbsyst
                write (ifm, '(I2,12(1X,E11.4))') i, (zr(ipr-1+decal+ &
                                                        nbsyst*(i-1)+j), j=1, nbsyst)
            end do
        end if
    else
        if (nboccm .gt. 1) then
            call utmess('F', 'COMPOR2_20', si=nbsyst)
        end if
    end if
    tabdes(1) = 1
    tabdes(2) = 1
    tabdes(3) = nvi
    tabdes(4) = itab
    tabdes(5) = nboccm
    tabdes(6) = irota
    tabdes(7) = nvi
    tabdes(8) = nbsyst
!
!     organisation de CPRI :
!     1 : TYPE =1 pour MONOCRISTAL
!     2 : NBPHAS=1 pour MONOCRISTAL
!     3 : NVI
!     4 : NOMBRE DE MONOCRISTAUX différents  =1
!     5 : NBFAMILLES DE SYS GLIS
!     6 : 1 si ROTA=POST, 2 si CALC, 0 sinon
!     7 : NVI
!     8 : NOMBRE DE SYSTEMES DE GLISSEMENT TOTAL
!
    call wkvect(sdcomp//'.CPRI', 'G V I', 13, imi)
    do i = 1, 13
        zi(imi+i-1) = tabdes(i)
    end do
    call detrsd('TABLE', comdes)
!
! FIN ------------------------------------------------------------------
!
    call jedema()
end subroutine
