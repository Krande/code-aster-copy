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
subroutine op0061()
!
!        MODE_NON_LINE
!        RECHERCHE DE MODES NON-LINEAIRES
!-----------------------------------------------------------------------
!      -POUR LE PROBLEME D'ELASTICITE LINEAIRE AVEC CONTACT LOCALISEE :
!                  ..
!              (M) U  + (K) U + F(U)=0
!              F(U) = [0 ... 0 F  0 ... 0 F   ... 0 ... 0]
!                             I1         I2
!               F  (F  - ALPHA*(U  - 1))(F  - ALPHA*(U  + 1))
!                I   I           I                    I
!          LES MATRICES (M) ET (K) SONT REELLES SYMETRIQUES
!-----------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8depi.h"
#include "asterfort/compno.h"
#include "asterfort/cremnl.h"
#include "asterfort/cresol.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlali.h"
#include "asterfort/mnlbhf.h"
#include "asterfort/mnlbra.h"
#include "asterfort/mnlcdl.h"
#include "asterfort/mnlcho.h"
#include "asterfort/mnlcof.h"
#include "asterfort/mnlcor.h"
#include "asterfort/mnleng.h"
#include "asterfort/mnlgen.h"
#include "asterfort/mnllec.h"
#include "asterfort/mnltan.h"
#include "asterfort/tbacce.h"
#include "asterfort/tbexve.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    character(len=6) :: nompro
    parameter(nompro='OP0061')
    character(len=8) :: baseno
    character(len=24) :: modele
!
!     -----------------------------------------------------------------
    data modele/'                        '/
!                   123456789012345678901234
!     -----------------------------------------------------------------
!
    integer(kind=8) :: imat(2), ordman, nbpt, nchoc, h, hf, itemax, nbran, nbranf, nextr
    real(kind=8) :: epsman, epscor, epsbif
    character(len=24) :: numedd
    character(len=14) :: numdrv
    character(len=19) :: matdrv, nomobj, nomvect, solveu, matas
    character(len=8) :: modini, modrep, vk8, typval, lnm, mailla
    character(len=24) :: masse, grno
!
    integer(kind=8) :: neq, nd, ijmax, iadim, ier, ninc, p, i, k, j, jj, ii, nbordr0, nbordr
    integer(kind=8) :: prodsci, iexi, jrefa
    character(len=14) :: xcdl, parcho, adime
    real(kind=8) :: ampl, amax, ap, epscor2, vr
    aster_logical :: cor, lbif, reprise, lcine
    integer(kind=8) :: iraid, ireg, iorig, ijeu, inddl, ifres, vi, num_ordr, num_lig, nbno, info
    integer(kind=8) :: ivec, iutj, iut1, iups, ius, ifpnl, ieng, isort, icdl, ivect, numrep, ntab
    character(len=14) :: xvect, xut1, xutj, xups, xus, xfpnl, xeng, xsort, xbif
!
    real(kind=8) :: omega, err, prodsc
    complex(kind=8) :: vc
!
    integer(kind=8) :: ibif
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------------------------
!
    call jemarq()
!
    ifres = iunifi('MESSAGE')
    call getvis(' ', 'INFO', scal=info)
    baseno = '&&'//nompro
!
    solveu = '&&OP0061.SOLVEUR'
    call cresol(solveu)
!
    call getvis('ETAT_INIT', 'NUME_ORDRE', iocc=1, scal=num_ordr)
    call getvid('ETAT_INIT', 'MODE_NON_LINE', iocc=1, scal=modrep, nbret=ier)
    if (ier .eq. 1) then
        reprise = .true.
        nomvect = '&&NUME_REPRISE'
        call tbexve(modrep, 'NUME_REPRISE', nomvect, 'V', ntab, &
                    typval)
        call jeveuo(nomvect, 'L', ivect)
        numrep = zi(ivect)
        do k = 2, ntab
            if (zi(ivect-1+k) .gt. numrep) numrep = zi(ivect-1+k)
        end do
        numrep = numrep+1
!
        nomvect = '&&NUME_ORDRE'
        call tbexve(modrep, 'NUME_ORDRE', nomvect, 'V', ntab, &
                    typval)
        call jeveuo(nomvect, 'L', ivect)
        num_lig = 0
        do k = 1, ntab
            if (zi(ivect-1+k) .eq. num_ordr) num_lig = k
        end do
        nbordr0 = zi(ivect-1+ntab)
        if (num_lig .eq. 0) then
            call utmess('F', 'RUPTURE0_37', si=num_ordr)
        end if
        call tbacce(modrep, num_lig, 'NOM_SD', 'L', vi, &
                    vr, vc, modini)
    else
        call getvid('ETAT_INIT', 'MODE_LINE', iocc=1, scal=lnm)
        reprise = .false.
        numrep = 0
        nbordr0 = 0
    end if
! ----------------------------------------------------------------------
! --- LECTURE DES DONNEES DE L'OPERATEUR
! ----------------------------------------------------------------------
    call mnllec(imat, numedd, ordman, epsman, nbpt, &
                epscor, h, hf, itemax, nbran, &
                nextr, epsbif)
! ----------------------------------------------------------------------
! --- TAILLE DE LA MATRICE
! ----------------------------------------------------------------------
    neq = zi(imat(1)+2)
! ----------------------------------------------------------------------
! --- GESTION DES CONDITIONS AUX LIMITES
! ----------------------------------------------------------------------
    xcdl = baseno//'.ICDL'
    call wkvect(xcdl, 'V V I', neq, icdl)
    call mnlcdl(imat, numedd, xcdl, nd, lcine)
! ----------------------------------------------------------------------
! --- GESTION DES PARAMETRES DE CHOC ET ADIMENSIONNEMENT
! ----------------------------------------------------------------------
!
    parcho = baseno//'.PARCH'
    if (reprise) then
        call tbacce(modrep, 1, 'CARA_CHOC', 'L', vi, &
                    vr, vc, vk8)
        nomobj = '&&NUME_CHOC'
        call tbexve(vk8, 'NUME_CHOC', nomobj, 'V', nchoc, &
                    typval)
    else
        call getfac('CHOC', nchoc)
        masse = zk24(zi(imat(2)+1)) (1:19)
        call dismoi('NOM_MAILLA', masse, 'MATR_ASSE', repk=mailla)
        do k = 1, nchoc
            call getvtx('CHOC', 'GROUP_NO', iocc=k, scal=grno, nbret=ier)
            if (ier .eq. 1) then
                call compno(mailla, 1, grno, nbno)
                if (nbno .gt. 1) then
                    call utmess('F', 'MECANONLINE9_68')
                end if
            end if
        end do
        call wkvect(parcho//'.TYPE', 'V V K8', nchoc, ier)
        call wkvect(parcho//'.NOEU', 'V V K8', 2*nchoc, ier)
        call wkvect(parcho//'.RAID', 'V V R', nchoc, iraid)
        call wkvect(parcho//'.REG', 'V V R', nchoc, ireg)
        call wkvect(parcho//'.JEU', 'V V R', nchoc, ijeu)
    end if
!
    call wkvect(parcho//'.JEUMAX', 'V V R', 1, ijmax)
    call wkvect(parcho//'.INDMAX', 'V V I', 1, ier)
    call wkvect(parcho//'.NDDL', 'V V I', 6*nchoc, inddl)
    call wkvect(parcho//'.NEQS', 'V V I', nchoc, ier)
    call wkvect(parcho//'.NCMP', 'V V I', nchoc, ier)
    call wkvect(parcho//'.CMP', 'V V K8', 2*nchoc, ier)
    call wkvect(parcho//'.ORIG', 'V V R', 3*nchoc, iorig)
    adime = baseno//'.ADIME'
    call wkvect(adime, 'V V R', 3, iadim)
!
    call mnlcho(reprise, imat, numedd, xcdl, nd, &
                nchoc, h, hf, parcho, adime, &
                ninc, vk8, lcine, solveu)
! ----------------------------------------------------------------------
! --- NOMBRE D'INCONNUES TOTAL DU SYSTEME
! ----------------------------------------------------------------------
    numdrv = baseno//'.NUDRV'
    matdrv = baseno//'     .MTDRV'
    call mnlgen(numdrv, matdrv, ninc)
! ----------------------------------------------------------------------
! --- AMPLITUDE DE DEMARRAGE DE L'ALGORITHME
! ----------------------------------------------------------------------
    call getvr8('ETAT_INIT', 'COEF_AMPL', iocc=1, scal=ampl)
!    ampl = 0.01d0
! ----------------------------------------------------------------------
! --- INITIALISATION DE L'ALGORITHME
! ----------------------------------------------------------------------
905 format('-- INITIALISATION :     ----------')
    xvect = baseno//'.XVECT'
    call wkvect(xvect, 'V V R', ninc, ivec)
    call mnlali(reprise, modini, imat, xcdl, parcho, &
                adime, ninc, nd, nchoc, h, &
                hf, ampl, xvect, lnm, num_ordr)
    write (ifres, 905)
! ----------------------------------------------------------------------
! --- CORRECTION DU POINT DE DEPART
! ----------------------------------------------------------------------
    cor = .true.
    epscor2 = 1.d-08
    call mnlcor(imat, numdrv, matdrv, xcdl, parcho, &
                adime, ninc, nd, nchoc, h, &
                hf, itemax, epscor, xvect, cor, &
                info)
! ----------------------------------------------------------------------
! --- DIRECTION DE LA CONTINUATION (TANGENTE AU PREMIER POINT)
! ----------------------------------------------------------------------
    xut1 = baseno//'.XUT1'
    xutj = baseno//'.XUTJ'
    call wkvect(xut1, 'V V R', ninc, iut1)
    call wkvect(xutj, 'V V R', ninc, iutj)
    call mnltan(.true._1, imat, numdrv, matdrv, xcdl, &
                parcho, adime, xvect, ninc, nd, &
                nchoc, h, hf, xut1)
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(iut1), b_incx, zr(iutj), b_incy)
    call getvis('ETAT_INIT', 'DIR_EVOLUTION', iocc=1, scal=prodsci)
    prodsc = dble(prodsci)
    if (prodsc .le. 0.d0) then
        b_n = to_blas_int(ninc)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.d0, zr(iut1), b_incx)
    end if
! ----------------------------------------------------------------------
! --- BOUCLE SUR LE NOMBRE DE BRANCHE
! ----------------------------------------------------------------------
    xups = baseno//'.XUPS'
    xus = baseno//'.XUS'
    xfpnl = baseno//'.XFPNL'
    xeng = baseno//'.ENG'
    xsort = baseno//'.SORTI'
    xbif = baseno//'.BIF'
    call wkvect(xfpnl, 'V V R', ninc-1, ifpnl)
    call wkvect(xups, 'V V R', ninc*(ordman+1), iups)
    call wkvect(xus, 'V V R', ninc*nbpt, ius)
    call wkvect(xeng, 'V V R', nbpt-1, ieng)
    call wkvect(xsort, 'V V R', nbran*(nbpt-1)*(neq*(2*h+1)+2), isort)
    call wkvect(xbif, 'V V I', nbran, ibif)
!
900 format('-- BRANCHE NUMERO : ', i8, '   ----------')
903 format('     DETECTION D UNE BIFURCATION')
904 format('             (   ENERGIE   ,  FREQUENCE  )')
901 format('     DEBUT - (', 1pe12.5, ' ,', 1pe12.5, ' )')
902 format('     FIN   - (', 1pe12.5, ' ,', 1pe12.5, ' )')
    k = 1
10  continue
    if (k .le. nbran) then
!    do 10 k = 1, nbran
! ---   CALCUL DES COEFFICIENTS DE LA SERIE ENTIERE
        b_n = to_blas_int((ordman+1)*ninc)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(iups), b_incx)
        call mnlcof(imat, numdrv, matdrv, xcdl, parcho, &
                    adime, xvect, xut1, ninc, nd, &
                    nchoc, h, hf, ordman, xups, &
                    xfpnl, lbif, nextr, epsbif)
        write (ifres, 900) k
        if (lbif) then
            if (info .eq. 2) then
                write (ifres, 903)
            end if
            zi(ibif-1+k) = 1
        end if
! ---   RECONSTITUTION D'UNE BRANCHE
        call mnlbra(xups, xfpnl, ninc, ordman, nbpt, &
                    epsman, amax, xus)
! ---   RECUPERATION DU DERNIER POINT DE LA BRANCHE POUR INITIALISATION
! ---   DE LA PROCHAINE BRANCHE
        b_n = to_blas_int(ninc)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ius+(nbpt-1)*ninc), b_incx, zr(ivec), b_incy)
! ---   CORRECTION DE CE POINT
        cor = .true.
        call mnlcor(imat, numdrv, matdrv, xcdl, parcho, &
                    adime, ninc, nd, nchoc, h, &
                    hf, itemax, epscor, xvect, cor, &
                    info)
        if (cor) then
! ---       DIRECTION DE LA CONTINUATION (TANGENTE DE V(A))
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(iutj), b_incx)
            do p = 1, ordman
                ap = dble(p)*(amax**(p-1))
                b_n = to_blas_int(ninc)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, ap, zr(iups+p*ninc), b_incx, zr(iutj), &
                           b_incy)
            end do
! ---       CALCUL DE LA TANGENTE AU NOUVEAU POINT
            call mnltan(.true._1, imat, numdrv, matdrv, xcdl, &
                        parcho, adime, xvect, ninc, nd, &
                        nchoc, h, hf, xut1)
! ---       SENS DE CONTINUATION
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            prodsc = ddot(b_n, zr(iutj), b_incx, zr(iut1), b_incy)
            if (prodsc .le. 0.d0) then
                b_n = to_blas_int(ninc)
                b_incx = to_blas_int(1)
                call dscal(b_n, -1.d0, zr(iut1), b_incx)
            end if
!
! ---       BON HF ?
!
            call mnlbhf(xvect, parcho, adime, ninc, nd, &
                        nchoc, h, hf, err)
            if (err .gt. 5.d-02) then
                call utmess('A', 'MECANONLINE9_60')
            end if
!
! ---       REDIMENSIONNEMENT DES VECTEURS SOLUTIONS
            do i = 1, nbpt
                do p = 1, nd*(2*h+1)
                    zr(ius-1+(i-1)*ninc+p) = zr(ius-1+(i-1)*ninc+p)*zr(ijmax)
                end do
            end do
            do i = 1, nbpt-1
                zr(ius-1+i*ninc) = zr(ius-1+i*ninc)*zr(iadim-1+3)
            end do
! ---       CALCUL DE L'ENERGIE MECANIQUE
            call mnleng(imat, xcdl, parcho, xus, ninc, &
                        nd, nchoc, h, nbpt, xeng)
            do i = 1, nbpt-1
                omega = zr(ius-1+i*ninc)
            end do
            if (info .eq. 2) then
                write (ifres, 904)
                write (ifres, 901) zr(ieng), zr(ius-1+ninc)/r8depi()
                write (ifres, 902) zr(ieng-1+(nbpt-1)), zr(ius-1+(nbpt-1)*ninc)/r8depi()
            end if
! ---       AJOUTER LES DDLS NON-ACTIFS
            do ii = 1, nbpt-1
                p = (k-1)*(nbpt-1)*(neq*(2*h+1)+2)
                p = p+(ii-1)*(neq*(2*h+1)+2)
                do j = 1, 2*h+1
                    i = 0
                    do jj = 1, neq
                        if (zi(icdl-1+jj) .eq. 0) then
                            i = i+1
                            zr(isort-1+p+(j-1)*neq+jj) = zr(ius-1+(ii-1)*ninc+(j-1)*nd+i)
                        else
                            zr(isort-1+p+(j-1)*neq+jj) = 0.d0
                        end if
                    end do
                end do
                zr(isort-1+p+neq*(2*h+1)+1) = zr(ius-1+ii*ninc)/r8depi()
                zr(isort-1+p+neq*(2*h+1)+2) = zr(ieng-1+ii)
            end do
! ---       REMISE A ZERO DU VECTEUR RESULTAT
            if (k .ne. nbran) then
                b_n = to_blas_int(ninc*nbpt)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(ius), b_incx)
            end if
            nbranf = nbran
            k = k+1
        else
            nbranf = k-1
            k = nbran+1
        end if
        goto 10
    end if
!10  continue
!
    nbordr = nbranf*(nbpt-1)
    call cremnl(reprise, baseno, numrep, nbordr0, nbordr, &
                nbpt, neq, h, imat, numedd, &
                parcho, nchoc, vk8, modrep)
    if (.not. cor) then
        call utmess('S', 'MECANONLINE9_62')
    end if
!
!   -- pour etre OK sdveri, il ne faut pas faire reference au solveur
!      temporaire :
!   -------------------------------------------------------------------
    call getvid(' ', 'MATR_RIGI', scal=matas)
    call jeexin(matas//'.REFA', iexi)
    if (iexi .gt. 0) then
        call jeveuo(matas//'.REFA', 'E', jrefa)
        zk24(jrefa-1+7) = ' '
    end if
    call getvid(' ', 'MATR_MASS', scal=matas)
    call jeexin(matas//'.REFA', iexi)
    if (iexi .gt. 0) then
        call jeveuo(matas//'.REFA', 'E', jrefa)
        zk24(jrefa-1+7) = ' '
    end if
!
    call jedema()
end subroutine
