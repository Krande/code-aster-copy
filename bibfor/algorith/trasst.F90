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
subroutine trasst(modgen, numsst, isst1, lisint, nbeq1, &
                  nbmod, nbint, solveu)
!
!
!-------------------------------------------------------------C
!--       ROUTINE XXXXX2           M. CORUS - AOUT 2011     --C
!--       CALCUL DES TRAVAUX DANS LES SOUS STRUCTURES       --C
!--                                                         --C
!-------------------------------------------------------------C
!--   VARIABLES E/S  :
!--   MODGEN   /IN/  : NOM DU MODELE GENERALISE
!--   NUMSST   /IN/  : NUMERO DE LA SOUS STRUCTURE TRAITEE
!--   ISST1    /IN/  : NUMERO DE LA SOUS STRUCTURE
!--   LISINT   /IN/  : LISTE DES NOMS D'INTERFACES
!--   NBEQ1    /IN/  : NB DE DDL DE LA SST
!--   NBMOD    /IN/  : NOMBRE DE MODE DU MODELE REDUIT
!--   NBINT    /IN/  : NB D'INTERFACES ASSOCIEES A LA SST
!
    implicit none
!
!
!
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lceqvn.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
#include "asterfort/zerlag.h"
#include "blas/daxpy.h"
#include "blas/ddot.h"
!
    character(len=1) :: listyp(2)
    character(len=4) :: k4bid
    character(len=8) :: modgen, rest1, mraid, mmass
    character(len=19) :: imped, lismat(2), nume91, solveu
    character(len=24) :: indin1
    integer(kind=8) :: i1, ibid, iret, j1, k1, l1, nbeq1, nbmod, isst1, llint1, nbddl1
    integer(kind=8) :: tach1, lmod1, lbid, leff1
    integer(kind=8) :: lintf, nbint, lcopy1, lsecme, limped, unit, numsst
    real(kind=8) :: travm, travk, traint, comlin(2), shift, pi
    character(len=24) :: lisint
    real(kind=8), pointer :: mode_sst1_eff2(:) => null()
    real(kind=8), pointer :: pulsa_propres(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: matrice_mass(:) => null()
    real(kind=8), pointer :: trav_sst(:) => null()
    integer(kind=8), pointer :: matrice_raid(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call getvis(' ', 'UNITE', scal=unit, nbret=ibid)
    i1 = numsst
    pi = r8pi()
!
!-- RECHERCHE DU MACRO ELEMENT ASSOCIE A LA SST
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
    call jeveuo(zk8(ibid)//'      .NUME.DEEQ', 'L', vi=deeq)
!
!------------------------------------------------------------C
!--                                                        --C
!-- CONSTRUCTION DES MATRICES D'IMPEDANCE DYNAMIQUE K+MU*M --C
!--                   POUR L'ENRICHISSEMENT                --C
!--                                                        --C
!------------------------------------------------------------C
!
    call codent(numsst, 'D0', k4bid)
    imped = '&&OP0091.IMPED'//k4bid
!
    call jeveuo(jexnum(modgen//'      .MODG.SSME', isst1), 'L', ibid)
!
    call jeveuo(zk8(ibid)//'.MAEL_MASS_REFE', 'L', lbid)
    mmass = zk24(lbid+1) (1:8)
    call jeveuo(zk8(ibid)//'.MAEL_RAID_REFE', 'L', lbid)
    mraid = zk24(lbid+1) (1:8)
    call mtdefs(imped, mmass, 'V', ' ')
    lismat(1) = mraid
    lismat(2) = mmass
!
    call dismoi('NOM_NUME_DDL', mraid, 'MATR_ASSE', repk=nume91)
!
    call getvr8(' ', 'SHIFT', scal=shift, nbret=ibid)
    comlin(1) = 1.d0
    comlin(2) = -((shift*2.d0*pi)**2)
    listyp(1) = 'R'
    listyp(2) = 'R'
    call mtcmbl(2, listyp, comlin, lismat, imped, &
                ' ', nume91, 'ELIM1')
    call mtdscr(imped)
    call jeveuo(imped(1:19)//'.&INT', 'E', limped)
!
!
    call preres(solveu, 'V', iret, '&&OP0091.MATPRE', imped, &
                ibid, -9999)
    if (iret .eq. 2) then
        call utmess('F', 'ALGELINE4_37', sk=imped)
    end if
!
    rest1 = '&&91'//k4bid
    call jeveuo(jexnum(rest1//'           .TACH', 1), 'L', tach1)
    call jeveuo('&&OP0091.MODE_SST1', 'E', lmod1)
    call jeveuo('&&OP0091.MODE_SST1_EFF1', 'E', leff1)
    call jeveuo('&&OP0091.MODE_SST1_EFF2', 'E', vr=mode_sst1_eff2)
    call jeveuo('&&OP0091.MODE_SST1_COPY', 'E', lcopy1)
    call jeveuo(lisint, 'L', lintf)
!
    call jeveuo('&&OP0091.MATRICE_MASS', 'L', vi=matrice_mass)
    call jeveuo('&&OP0091.MATRICE_RAID', 'L', vi=matrice_raid)
    call jeveuo('&&OP0091.TRAV_SST', 'E', vr=trav_sst)
    call jeveuo('&&OP0091.PULSA_PROPRES', 'L', vr=pulsa_propres)
    call jeveuo('&&OP0091.MODE_INTF_DEPL', 'E', lsecme)
!
!-- BOUCLE SUR LES MODES
    do j1 = 1, nbmod
        call jeveuo(zk24(tach1+j1-1) (1:19)//'.VALE', 'L', ibid)
!
!-- RECOPIE DANS UN VECTEUR DE TRAVAIL
        call lceqvn(nbeq1, zr(ibid), zr(lcopy1))
!
!-- ANNULATION DES DDL DE LAGRANGE
        call zerlag(nbeq1, deeq, vectr=zr(lcopy1))
!
!-- NOUVELLE COPIE
        call lceqvn(nbeq1, zr(lcopy1), zr(lmod1))
!
!-- ANNULATION DES COMPOSANTES ASSOCIEES AUX INTERFACES
        do k1 = 1, nbint
            indin1 = '&&VEC_DDL_INTF_'//zk8(lintf+k1-1)
            call jeveuo(indin1, 'L', llint1)
            call jelira(indin1, 'LONMAX', nbddl1)
            do l1 = 1, nbddl1
                if (zi(llint1+l1-1) .gt. 0) then
                    zr(lmod1+zi(llint1+l1-1)-1) = 0
                end if
            end do
        end do
!
!-- CALCUL DES TRAVAUX
        call mrmult('ZERO', matrice_mass(isst1), zr(lcopy1), zr(leff1), 1, &
                    .true._1)
!
        b_n = to_blas_int(nbeq1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travm = ddot(b_n, zr(lmod1), b_incx, zr(leff1), b_incy)
        call mrmult('ZERO', matrice_raid(isst1), zr(lcopy1), mode_sst1_eff2, 1, &
                    .true._1)
        b_n = to_blas_int(nbeq1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        travk = ddot(b_n, zr(lmod1), b_incx, mode_sst1_eff2, b_incy)
        traint = travk-(pulsa_propres(j1)**2)*travm
        if (pulsa_propres(j1) .gt. 1) traint = traint/pulsa_propres(j1)
        write (unit, *) 'MODE ', j1, ' -  TRAVAIL SST =', traint
        trav_sst(1+nbmod*(i1-1)+j1-1) = traint
!
!--
!-- CALCUL DU SECOND MEMBRE ET DES ENRICHISSEMENTS
!--
        b_n = to_blas_int(nbeq1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -(pulsa_propres(j1)**2), zr(leff1), b_incx, mode_sst1_eff2, &
                   b_incy)
        call zerlag(nbeq1, deeq, vectr=zr(leff1))
        lbid = lsecme
        call lceqvn(nbeq1, zr(leff1), zr(lsecme+nbeq1*(j1-1)))
!
!-- DIFFERENTIATION DES SECONDS MEMBRES : INTERFACE / INTERIEUR
!
        do k1 = 1, nbint
            indin1 = '&&VEC_DDL_INTF_'//zk8(lintf+k1-1)
            call jeveuo(indin1, 'L', llint1)
            call jelira(indin1, 'LONMAX', nbddl1)
            do l1 = 1, nbddl1
                ibid = zi(llint1+l1-1)
                if (ibid .gt. 0) then
                    zr(lsecme+nbeq1*(j1-1)+ibid-1) = 0
                    zr(lsecme+nbeq1*(nbmod+j1-1)+ibid-1) = zr(leff1+ibid-1)
                end if
            end do
        end do
!
!
    end do
!
end subroutine
