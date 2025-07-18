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
subroutine cremnl(reprise, baseno, numrep, nbordr0, nbordr, &
                  nbpt, neq, nbhar, imat, numedd, &
                  parcho, nbchoc, vk8, modrep)
!
!
! IN REPRISE: REPRISE DE CALCUL (LOGICAL)
! IN BASENO : BASE NOM DES OBJETS DE TRAVAIL
! IN NBORDR : NOMBRE DE NUMEROS D ORDRE
! IN NBPT   : NOMBRE DE POINT DE DISCRETISATION
! IN NEQ    : NB EQUATIONS
! IN NBHAR  : NB HARMONIQUES
! IN IMAT   : DESCRIPTEUR DES MATRICES
!             IMAT(1) : MATRICE DE RIGIDITE
!             IAMT(2) : MATRICE DE MASSE
! IN NUMEDD : NUME_DDL DES MATRICES DU SYSTEME
! IN PARCHO : SD VOLATILE OU SONT STOCKES LES PARAMETRES DE CHOC
! IN NBCHOC : NB DE LIEUX DE CHOC
!     -----------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcncon.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/refdaj.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrem.h"
#include "blas/dnrm2.h"
#include "blas/idamax.h"
!     -----------------------------------------------------------------
!
    aster_logical :: reprise, suite
    character(len=4) :: nomsym(1), kordr
    character(len=8) :: nomres, nomrep, baseno, k8b, nomtab, modrep
    character(len=8) :: typpar(12), typpat(11), vk8
    character(len=14) :: parcho
    character(len=16) :: nomcmd, typres, valk(6), kbif
    character(len=19) :: nompar(12), chamno
    character(len=19) :: nompat(11)
    character(len=24) :: numedd, rigid, masse, valkt(4), matrice(3)
    integer(kind=8) :: iordr, numrep, nbordr0, nbordr, nbpt, vali(4), nbhar, nbsym, isort
    integer(kind=8) :: neq, imat(2), iadd, ieq, ier, ihar, ladpa
    integer(kind=8) :: nbpar, neqv, nmodes, iomega, valit(1), nbpart
    integer(kind=8) :: nbchoc
    integer(kind=8) :: harmaxa, harmax, inspec
    real(kind=8) :: rvide, valr(2), freq, energ, xnorm, valrt(6)
    real(kind=8) :: espec(nbhar+1)
    complex(kind=8) :: cvide
    character(len=8), pointer :: cmp(:) => null()
    real(kind=8), pointer :: orig(:) => null()
    real(kind=8), pointer :: nspec(:) => null()
    integer(kind=8), pointer :: bif(:) => null()
    character(len=8), pointer :: noeu(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    character(len=8), pointer :: type(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    real(kind=8), pointer :: jeu(:) => null()
    real(kind=8), pointer :: reg(:) => null()
    blas_int :: b_incx, b_n
!     -----------------------------------------------------------------
    call jemarq()
    call getres(nomres, typres, nomcmd)
!
    if (reprise) then
        if (nomres .ne. modrep) then
            suite = .true.
        else
            suite = .false.
        end if
    else
        suite = .false.
    end if
!
    rvide = r8vide()
    cvide = dcmplx(rvide, rvide)
!
    nbpar = 12
!
    nompar(1) = 'NUME_ORDRE'
    typpar(1) = 'I'
    nompar(2) = 'FREQUENCE'
    typpar(2) = 'R'
    nompar(3) = 'ENERGIE'
    typpar(3) = 'R'
    nompar(4) = 'NB_COEF_FOURIER'
    typpar(4) = 'I'
    nompar(5) = 'NOM_SD'
    typpar(5) = 'K8'
    nompar(6) = 'TYPE_OBJET'
    typpar(6) = 'K16'
    nompar(7) = 'NOM_OBJET'
    typpar(7) = 'K16'
    nompar(8) = 'CARA_CHOC'
    typpar(8) = 'K8'
    nompar(9) = 'HARM_MAX'
    typpar(9) = 'I'
    nompar(10) = 'BIFURCATION'
    typpar(10) = 'K16'
    nompar(11) = 'STABILITE'
    typpar(11) = 'K16'
    nompar(12) = 'NUME_REPRISE'
    typpar(12) = 'I'
!
    if ((.not. reprise) .or. (suite)) then
!-- INITIALISATION DE LA TABLE_CONTAINER
        call detrsd('TABLE_CONTAINER', nomres)
        call tbcrsd(nomres, 'G')
        call tbajpa(nomres, nbpar, nompar, typpar)
        call gcncon('_', nomtab)
    end if
!
    nbsym = 1
    nomsym(1) = 'DEPL'
!
    rigid = zk24(zi(imat(1)+1))
    masse = zk24(zi(imat(2)+1))
    neqv = zi(imat(1)+2)
    if (neq .ne. neqv) then
        call utmess('F', 'UTILITAI3_21')
    end if
!
    call jeveuo(baseno//'.SORTI', 'L', isort)
    call jeveuo(baseno//'.BIF', 'L', vi=bif)
!
    nmodes = 2*nbhar+1
!
    AS_ALLOCATE(vr=nspec, size=neq)
!
!     BOUCLE SUR LES NUMEROS D ORDRE
    do iordr = 1, nbordr
!   Conversion du numéro d'ordre en chaine de caractère
        ASSERT(nbordr <= 9999)
        write (kordr, '(I4.4)') iordr
!
!-- ATTRIBUTION D UN NOM DE CONCEPT
        call gcncon('_', nomrep)
!
        call rscrsd('G', nomrep, 'MODE_MECA', nmodes)
        matrice(1) = rigid
        matrice(2) = masse
        matrice(3) = ' '
        call refdaj('F', nomrep, nmodes, numedd, 'DYNAMIQUE', &
                    matrice, ier)
        do ihar = 1, nmodes
            call rsexch(' ', nomrep, nomsym(1), ihar, chamno, &
                        ier)
            if (ier .eq. 0) then
            else if (ier .eq. 100) then
                call vtcrem(chamno, rigid, 'G', 'R')
            else
                call utmess('F', 'ALGELINE3_11')
            end if
!
            call jeveuo(chamno//'.VALE', 'E', vr=vale)
            xnorm = 0.d0
            do ieq = 1, neq
                iadd = (iordr-1)*(neq*nmodes+2)+(ihar-1)*neq+ieq
                vale(ieq) = zr(isort-1+iadd)
                xnorm = xnorm+zr(isort-1+iadd)*zr(isort-1+iadd)
            end do
            call rsnoch(nomrep, nomsym(1), ihar)
!
            call rsadpa(nomrep, 'E', 1, 'NUME_MODE', ihar, &
                        0, sjv=ladpa, styp=k8b)
            zi(ladpa) = ihar
!
            call rsadpa(nomrep, 'E', 1, 'FREQ', ihar, &
                        0, sjv=ladpa, styp=k8b)
!     -----------------------------------------------------------------
! reverifier le rangement des contributions harmoniques
! ici : 0 C1 C2  ... CH S1 S2 ... SH
!     -----------------------------------------------------------------
            if (ihar .eq. 1) then
                iomega = 0
            else if (ihar .le. nbhar+1) then
                iomega = ihar-1
            else
                iomega = ihar-(nbhar+1)
            end if
!
            iadd = (iordr-1)*(neq*nmodes+2)+nmodes*neq+1
            zr(ladpa) = iomega*zr(isort-1+iadd)
!
            call rsadpa(nomrep, 'E', 1, 'FACT_PARTICI_DX', ihar, &
                        0, sjv=ladpa, styp=k8b)
            zr(ladpa) = sqrt(xnorm)
        end do
!
        iadd = (iordr-1)*(neq*nmodes+2)+neq*nmodes+1
        freq = zr(isort-1+iadd)
!
        iadd = (iordr-1)*(neq*nmodes+2)+neq*nmodes+2
        energ = zr(isort-1+iadd)
!
        call jelira(nomrep//'           .ORDR', 'LONUTI', nmodes, k8b)
!
!       DETERMINER HARMONIQUE MAX
        harmax = 1
        do ieq = 1, neq
            iadd = (iordr-1)*(neq*nmodes+2)+ieq
            b_n = to_blas_int(2*nbhar+1)
            b_incx = to_blas_int(neq)
            nspec(ieq) = dnrm2(b_n, zr(isort-1+iadd), b_incx)**2
        end do
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        inspec = idamax(b_n, nspec, b_incx)
        iadd = (iordr-1)*(neq*nmodes+2)+inspec
        espec(1) = zr(isort-1+iadd)**2/nspec(inspec)
        do ihar = 1, nbhar
            espec(ihar+1) = ( &
                            zr(isort-1+iadd+ihar*neq)**2+zr(isort-1+iadd+(nbhar+ihar)*neq)**2)/n&
                            &spec(inspec &
                            )
        end do
        b_n = to_blas_int(nbhar+1)
        b_incx = to_blas_int(1)
        harmaxa = idamax(b_n, espec, b_incx)-1
        if (harmaxa .gt. harmax) then
            harmax = harmaxa
        end if
        if (bif(1+int(iordr/(nbpt-1))) .eq. 1) then
            kbif = 'OUI'
        else
            kbif = 'NON'
        end if
!
        if ((reprise) .and. (.not. suite)) then
            vali(1) = iordr+nbordr0
        else
            vali(1) = iordr
        end if
        vali(2) = nmodes
        vali(3) = harmax
        vali(4) = numrep
        valr(1) = freq
        valr(2) = energ
        valk(1) = nomrep
        valk(2) = 'MODE_MECA'
        valk(3) = 'FOURIER'//"_"//kordr
        if (reprise) then
            valk(4) = vk8
        else
            valk(4) = nomtab
        end if
        valk(5) = kbif
        valk(6) = 'NON_EVALUE'
        call tbajli(nomres, nbpar, nompar, vali, valr, &
                    [cvide], valk, 0)
!
    end do
!
    AS_DEALLOCATE(vr=nspec)
!
    if (.not. reprise) then
! TABLE POUR LES CARACTERISTIQUES DE CHOC
        call tbcrsd(nomtab, 'G')
        nbpart = 11
        call jeveuo(parcho//'.RAID', 'L', vr=raid)
        call jeveuo(parcho//'.JEU', 'L', vr=jeu)
        call jeveuo(parcho//'.REG', 'L', vr=reg)
        call jeveuo(parcho//'.TYPE', 'L', vk8=type)
        call jeveuo(parcho//'.NOEU', 'L', vk8=noeu)
        call jeveuo(parcho//'.CMP', 'L', vk8=cmp)
        call jeveuo(parcho//'.ORIG', 'L', vr=orig)
!
        nompat(1) = 'NUME_CHOC'
        typpat(1) = 'I'
        nompat(2) = 'TYPE_CHOC'
        typpat(2) = 'K8'
        nompat(3) = 'NOEUD_CHOC'
        typpat(3) = 'K8'
        nompat(4) = 'NOM_CMP_1'
        typpat(4) = 'K8'
        nompat(5) = 'NOM_CMP_2'
        typpat(5) = 'K8'
        nompat(6) = 'RIGI_NOR'
        typpat(6) = 'R'
        nompat(7) = 'PARA_REGUL'
        typpat(7) = 'R'
        nompat(8) = 'JEU'
        typpat(8) = 'R'
        nompat(9) = 'ORIG_OBST_X'
        typpat(9) = 'R'
        nompat(10) = 'ORIG_OBST_Y'
        typpat(10) = 'R'
        nompat(11) = 'ORIG_OBST_Z'
        typpat(11) = 'R'
        call tbajpa(nomtab, nbpart, nompat, typpat)
!
        do iordr = 1, nbchoc
            valit(1) = iordr
            valkt(1) = type(iordr)
            valkt(2) = noeu(iordr)
            valkt(3) = cmp(2*(iordr-1)+1)
            valkt(4) = cmp(2*(iordr-1)+2)
            valrt(1) = raid(iordr)
            valrt(2) = reg(iordr)
            valrt(3) = jeu(iordr)
            valrt(4) = orig((iordr-1)*3+1)
            valrt(5) = orig((iordr-1)*3+2)
            valrt(6) = orig((iordr-1)*3+3)
            call tbajli(nomtab, nbpart, nompat, valit, valrt, &
                        [cvide], valkt, 0)
        end do
    end if
!
    call jedema()
end subroutine
