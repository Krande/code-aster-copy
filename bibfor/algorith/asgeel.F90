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
subroutine asgeel(nomres, option, nugene)
!
    implicit none
!    M. CORUS     DATE 26/02/10
!-----------------------------------------------------------------------
!  BUT:      < ASSEMBLAGE GENERALISEE DE MATRICE PLEINE ET ELIMINATION >
!
!  ASSEMBLER UNE MATRICE A PARTIR D'UNE NUMEROTATION GENERALISEE
!  ET D'UNE OPTION (RIGI_GENE,MASS_GENE,AMOR_GENE)
!
! REMARQUE : L'ASSEMBLAGE DONNE UNE MATRICE ASSEMBLEE PLEINE
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM K8 DE LA MATRICE GENERALISEE RESULTAT
! OPTION   /I/: OPTION DE CALCUL (RIGI_GENE,MASS_GENE)
! NUGENE   /I/: NOM K14 DE LA NUMEROTATION GENERALISEE
! STOPLE   /I/: NOM K19 DU STOCKAGE DE LA MATRICE (PLEIN)
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
!
!
!
    character(len=8) :: nomres, modgen
    character(len=9) :: rigopt, masopt, amoopt
    character(len=14) :: nugene
    character(len=19) :: prgene
    character(len=10) :: adnom
    character(len=11) :: option
    character(len=24) :: k24bid, seliai, sizlia, sst, nomsst
    real(kind=8) :: temp
!
    integer(kind=8) :: ibid, i1, j1, k1, l1, n1, lres, neq, lneq, lproj, nbddl, indsst
    integer(kind=8) :: lsst, llref, nbmacr, exist, lmacr, indmcr, lresi
    integer(kind=8) :: decal, ddlcp, lselia, nlt, ind, ldconl, iret, nmaddl, ltemp
    integer(kind=8) :: jrefa, nbsst, lsilia, iadesc, ntria, lmacri
    character(len=8) :: kbid, k8bid
    integer(kind=8), pointer :: indices_macro(:) => null()
    character(len=8), pointer :: noms_macro(:) => null()
    integer(kind=8), pointer :: pointeurs_macro(:) => null()
    integer(kind=8), pointer :: pointeurs_macroi(:) => null()
    aster_logical :: lsym
    aster_logical, pointer :: syme_macro(:) => null()
    blas_int :: b_incx, b_incy, b_n
    data rigopt, masopt, amoopt/'RIGI_GENE',&
     &                        'MASS_GENE', 'AMOR_GENE'/
!-----------C
!--       --C
!-- DEBUT --C
!--       --C
!-----------C
!
    call jemarq()
!
!----------------------------------C
!--                              --C
!-- INITIALISATION DE LA MATRICE --C
!--                              --C
!----------------------------------C
!
    prgene = nugene//'.NUME'
!
!
    call wkvect(nomres//'           .REFA', 'G V K24', 20, jrefa)
    zk24(jrefa-1+11) = 'MPI_COMPLET'
    zk24(jrefa-1+1) = ' '
    zk24(jrefa-1+2) = nugene
    zk24(jrefa-1+8) = 'ASSE'
    zk24(jrefa-1+9) = 'MS'
    zk24(jrefa-1+10) = 'GENE'
    lsym = .true.
!
!--------------------RECUPERATION DU MODE_GENE AMONT--------------------
!
    k24bid = prgene//'.REFN'
    call jeveuo(k24bid, 'L', llref)
    modgen = zk24(llref) (1:8)
    nomsst = modgen//'      .MODG.SSNO'
!
!------------------RECUPERATION DU NOMBRE DE SOUS-STRUCTURE-------------
!
    call jelira(nomsst, 'NOMMAX', nbsst)
!
!--
!
!      SELIAI= '&&'//NUGENE(1:8)//'PROJ_EQ_LIAI'
!      SIZLIA='&&'//NUGENE(1:8)//'VECT_SIZE_SS'
!      SST=    '&&'//NUGENE(1:8)//'VECT_NOM_SS'
    seliai = nugene(1:14)//'.ELIM.BASE'
    sizlia = nugene(1:14)//'.ELIM.TAIL'
    sst = nugene(1:14)//'.ELIM.NOMS'
!
    k24bid = prgene//'.NEQU'
    call jeveuo(k24bid, 'L', lneq)
    neq = zi(lneq)
!
    call wkvect(nomres//'           .CONL', 'G V R', neq, ldconl)
    do i1 = 1, neq
        zr(ldconl+i1-1) = 1.d0
    end do
!
!-- RECUPERATION DES DIFFERENTS MACRO ELEMENTS
!
    if (option .eq. rigopt) then
        adnom = '.MAEL_RAID'
    else if (option .eq. masopt) then
        adnom = '.MAEL_MASS'
    else if (option .eq. amoopt) then
        adnom = '.MAEL_AMOR'
    end if
!
    AS_ALLOCATE(vk8=noms_macro, size=nbsst)
    AS_ALLOCATE(vi=indices_macro, size=nbsst*3)
    AS_ALLOCATE(vl=syme_macro, size=nbsst)
!
    call jeveuo(sizlia, 'L', lsilia)
    call jeveuo(sst, 'L', lsst)
!
    ddlcp = 0
    nlt = 0
    nmaddl = 0
    do i1 = 1, nbsst
        nbddl = zi(lsilia+i1-1)
        if (nbddl .gt. nmaddl) then
            nmaddl = nbddl
        end if
        nlt = nlt+nbddl
        ddlcp = ddlcp+((nbddl*(nbddl+1))/2)
        kbid = '        '
        indsst = i1
        call jenonu(jexnom(nomsst, zk8(lsst+i1-1)), indsst)
        call mgutdm(modgen, kbid, indsst, 'NOM_MACR_ELEM', ibid, &
                    k8bid)
        call jelira(k8bid//adnom//'_VALE', 'NMAXOC', ntria)
        syme_macro(i1) = .true.
        if (ntria .gt. 1) then
            syme_macro(i1) = .false.
            lsym = .false.
            zk24(jrefa-1+9) = 'MR'
        end if
!
        if (i1 .eq. 1) then
            nbmacr = 1
            indices_macro(1) = 1
            indices_macro(nbsst+1) = nbddl
            indices_macro(1+2*nbsst) = nlt-nbddl
            noms_macro(1) = k8bid
        else
            exist = 0
            do j1 = 1, nbmacr
                if (noms_macro(j1) .eq. k8bid) then
                    exist = 1
                    indices_macro(i1) = j1
                    indices_macro(1+nbsst+i1-1) = nbddl
                    indices_macro(1+2*nbsst+i1-1) = nlt-nbddl
                end if
            end do
            if (exist .eq. 0) then
                nbmacr = nbmacr+1
                indices_macro(i1) = nbmacr
                indices_macro(1+nbsst+i1-1) = nbddl
                indices_macro(1+2*nbsst+i1-1) = nlt-nbddl
                noms_macro(nbmacr) = k8bid
            end if
        end if
    end do
!
!
!-- RECUPERATION DES MATRICES DES DIFERENTS MACROS ELEMENTS
!
    AS_ALLOCATE(vi=pointeurs_macro, size=nbmacr)
    AS_ALLOCATE(vi=pointeurs_macroi, size=nbmacr)
    do i1 = 1, nbmacr
        call jeveuo(jexnum(noms_macro(i1)//adnom//'_VALE', 1), 'L', lmacr)
        pointeurs_macro(i1) = lmacr
        if (syme_macro(i1)) then
            pointeurs_macroi(i1) = lmacr
        else
            call jeveuo(jexnum(noms_macro(i1)//adnom//'_VALE', 2), 'L', lmacri)
            pointeurs_macroi(i1) = lmacri
        end if
    end do
!
! ----------- CREATION ET REMPLISSAGE DU .DESC ---------------
    call wkvect(nomres//'           .DESC', 'G V I', 3, iadesc)
    zi(iadesc) = 2
    zi(iadesc+1) = neq
    zi(iadesc+2) = 2
!
!-- .VALM NE DOIT PAS EXISTER :
    call jeexin(nomres//'           .VALM', iret)
    ASSERT(iret .eq. 0)
!
!-- ALLOCATION DE LA MATRICE PROJETEE
    call jecrec(nomres//'           .VALM', 'G V R', 'NU', 'DISPERSE', 'CONSTANT', &
                ntria)
    call jeecra(nomres//'           .VALM', 'LONMAX', int((neq*(neq+1))/2))
!
    call jecroc(jexnum(nomres//'           .VALM', 1))
    call jeveuo(jexnum(nomres//'           .VALM', 1), 'E', lres)
    if (.not. lsym) then
        call jecroc(jexnum(nomres//'           .VALM', 2))
        call jeveuo(jexnum(nomres//'           .VALM', 2), 'E', lresi)
    end if
!
!----------------------------------------C
!--                                    --C
!-- REMPLISSAGE DE LA MATRICE PROJETEE --C
!--                                    --C
!----------------------------------------C
!
    call jeveuo(seliai, 'L', lselia)
!
    call wkvect('&&ASGEEL.MATR_TEMP', 'V V R', nmaddl**2, ltemp)
    call wkvect('&&ASGEEL.PROJ_TEMP', 'V V R', nmaddl*neq, lproj)
!
    do n1 = 1, nbsst
        indmcr = indices_macro(n1)
        nbddl = indices_macro(1+nbsst+n1-1)
        decal = indices_macro(1+2*nbsst+n1-1)
        lmacr = pointeurs_macro(indmcr)
        lmacri = pointeurs_macroi(indmcr)
!
!-- ON RECOPIE LA MATRICE DU MACRO ELEMENT (STOCKAGE PLEIN)
        do j1 = 1, int(nbddl*(nbddl+1)/2)
            temp = (sqrt(1.d0+8.d0*j1)-1.d0)/2.d0
            l1 = int(temp)
            if (temp .gt. l1) then
                l1 = l1+1
            end if
            k1 = j1-int(l1*(l1-1)/2)
            zr(ltemp+nbddl*(k1-1)+l1-1) = zr(lmacr+j1-1)
            zr(ltemp+nbddl*(l1-1)+k1-1) = zr(lmacri+j1-1)
        end do
!
!-- ON FAIT K*T
        do i1 = 1, nbddl
            do j1 = 1, neq
                b_n = to_blas_int(nbddl)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(lproj+(j1-1)*nbddl+i1-1) = ddot( &
                                              b_n, zr(ltemp+(i1-1)*nbddl), b_incx, &
                                              zr(lselia+(j1-1)*nlt+decal), b_incy &
                                              )
            end do
        end do
!
!-- ON FAIT T^T*(K*T)
        do j1 = 1, neq
            do i1 = 1, j1
                ind = int(((j1-1)*j1)/2)+i1-1
                b_n = to_blas_int(nbddl)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(lres+ind) = zr(lres+ind)+ddot(b_n, zr(lproj+(j1-1)*nbddl), b_incx, zr(lselia+(&
                               &i1-1)*nlt+decal), b_incy)
                if (.not. lsym) then
                    b_n = to_blas_int(nbddl)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    zr(lresi+ind) = zr(lresi+ind)+ddot(b_n, zr(lproj+(i1-1)*nbddl), b_incx, zr(ls&
                                    &elia+(j1-1)*nlt+decal), b_incy)
                end if
            end do
        end do
!
    end do
!
    AS_DEALLOCATE(vi=pointeurs_macro)
    AS_DEALLOCATE(vi=pointeurs_macroi)
    AS_DEALLOCATE(vk8=noms_macro)
    AS_DEALLOCATE(vi=indices_macro)
    AS_DEALLOCATE(vl=syme_macro)
    call jedetr('&&ASGEEL.MATR_TEMP')
    call jedetr('&&ASGEEL.PROJ_TEMP')
!
!---------C
!--     --C
!-- FIN --C
!--     --C
!---------C
!      CALL MATIMP(NOMRES//'           ',8,'MATLAB')
    call jedema()
!
end subroutine
