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

subroutine jjldyn(imode, lmin, ltot)
    use logging_module, only: DEBUG, LOGLEVEL_MEM, is_enabled
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/hpdeallc.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/jermxd.h"
#include "asterfort/jxecro.h"
#include "asterfort/utgtme.h"
#include "asterfort/utmess.h"
#include "asterfort/utptme.h"
#include "asterfort/uttcpu.h"
    integer(kind=8) :: imode, lmin, ltot
! ----------------------------------------------------------------------
! LIBERE LES SEGMENTS DE VALEURS ALLOUES DYNAMIQUEMENT
!
! IN   IMODE :
!              =1 ON NE TRAITE QUE LA BASE VOLATILE
!              =2 ON NE TRAITE QUE LES OBJETS XA
!              =3 ON NE TRAITE QUE LES OBJETS XA DE LA BASE VOLATILE
!              SINON ON EXAMINE TOUTE LA MEMOIRE
! IN   LMIN  : TAILLE MINIMUM EN ENTIERS REQUISE
!              =< 0 ON LIBERE TOUT
!              =-2  ON LIBERE TOUT MAIS ON N'ACTUALISE PAS VMXDYN
! OUT  LTOT  : LONGUEUR CUMULEE EN ENTIERS DES SEGMENTS DESALLOUES
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadmi, iadmoc, iadyn, iadyoc, ibacol
    integer(kind=8) :: ibiadd, ibiadm, iblono, ic, idm
    integer(kind=8) :: isd, isdc, isf, ixdeso, ixiadd, ixiadm
    integer(kind=8) :: ixlono, j, jcara, jdate, jdocu, jgenr, jhcod
    integer(kind=8) :: jiacce, jiadd, jiadm, jindir, jj, jlong, jlono
    integer(kind=8) :: jltyp, jluti, jmarq, jorig, jrnom, jtype, k
    integer(kind=8) :: lonoi, lsv, ltypi, n, nbacce, ncla1, ncla2
    integer(kind=8) :: nmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &               jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &               kitlec(n), kitecr(n), kiadm(n),&
     &               iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &               dn2(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    common/jiacce/jiacce(n), nbacce(2*n)
    common/jindir/jindir(n)
    integer(kind=8) :: isstat
    common/iconje/isstat
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    integer(kind=8) :: icdyn, mxltot
    common/xdynje/icdyn, mxltot
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: datei
    common/iheuje/datei
    integer(kind=8) :: indiq_jjagod, indiq_jjldyn
    common/idagod/indiq_jjagod, indiq_jjldyn
! ----------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idlono
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3, idlono=8)
! ----------------------------------------------------------------------
    character(len=1) :: cgenr
    character(len=8) :: nomk(5)
!    CHARACTER*32   NOM32
    integer(kind=8) :: iaddi(2), lgs, nbioav(2)
    integer(kind=8) :: rang, nbproc, iret, iret2
    real(kind=8) :: valp(5), vx(3), v0
    real(kind=4) :: graine
    mpi_int :: mrank, msize
!
    data nomk/'COUR_JV ', 'RLQ_MEM ', 'VMSIZE  ', 'MEM_TOTA', 'LIMIT_JV'/
!
    indiq_jjldyn = 1
    call uttcpu('CPU.MEMD.1', 'DEBUT', ' ')
!
!     ON LISTE LES OBJETS ALLOUES DYNAMIQUEMENT EN BALAYANT
!     L'ENSEMBLE DES OBJETS, EN COMMENCANT PAR LA BASE VOLATILE
!
    icdyn = icdyn+1
    ltot = 0
    ncla1 = 1
    ncla2 = index(classe, '$')-1
    if (ncla2 .lt. 0) ncla2 = n
    if (imode .eq. 1 .or. imode .eq. 3) then
        ncla2 = index(classe, 'V')
        ncla1 = ncla2
    end if
    do ic = ncla2, ncla1, -1
        if (nreuti(ic) .eq. 0) goto 200
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        if (rang .ne. 0) then
            graine = (rang+1)*datei*1.5d0
            do i = 2, nreuti(ic)
                call random_number(graine)
                k = int(graine*i)+1
                j = indir(jindir(ic)+i)
                indir(jindir(ic)+i) = indir(jindir(ic)+k)
                indir(jindir(ic)+k) = j
            end do
        end if
!
        nbioav(1) = nbacce(2*ic-1)
        nbioav(2) = nbacce(2*ic)
        do jj = 1, nreuti(ic)
            j = indir(jindir(ic)+jj)
            iadmi = iadm(jiadm(ic)+2*j-1)
            if (iadmi .eq. 0) goto 205
            iadyn = iadm(jiadm(ic)+2*j)
            cgenr = genr(jgenr(ic)+j)
!            NOM32 = RNOM(JRNOM(IC)+J)
!
!    ISD DESIGNE LE STATUT DE LA COLLECTION
!        =U ON PASSE PAR LES ROUTINES HABITUELLES (JJALLC, JJLIDE)
!        =X ON TRAITE DIRECTEMENT
!
            isdc = iszon(jiszon+iadmi-1)/isstat
            if (cgenr .eq. 'X' .and. isdc .eq. 2) then
                ibacol = iadmi
                ixiadm = iszon(jiszon+ibacol+idiadm)
                ixiadd = iszon(jiszon+ibacol+idiadd)
                ixdeso = iszon(jiszon+ibacol+iddeso)
                ixlono = iszon(jiszon+ibacol+idlono)
                nmax = iszon(jiszon+ibacol+ivnmax)
                if (ixiadm .gt. 0) then
                    ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
                    ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                    do k = 1, nmax
                        iadmoc = iszon(jiszon+ibiadm-1+2*k-1)
                        iadyoc = iszon(jiszon+ibiadm-1+2*k)
                        if (iadyoc .ne. 0) then
                            idm = iadmoc-4
                            isd = iszon(jiszon+idm+3)/isstat
                            isf = iszon(jiszon+iszon(jiszon+idm)-4)/isstat
                            if (isd .eq. 1) then
!
!     LE SEGMENT DE VALEURS EST MARQUE X A OU X D, ON PEUT LE LIBERER
!
                                if (ixlono .ne. 0) then
                                    iblono = iadm(jiadm(ic)+2*ixlono-1)
                                    lonoi = iszon(jiszon+iblono-1+k)
                                else
                                    lonoi = lono(jlono(ic)+ixdeso)
                                end if
                                ltypi = ltyp(jltyp(ic)+ixdeso)
                                lsv = lonoi*ltypi
                                if (isf .eq. 4) then
                                    if (imode .eq. 2 .or. imode .eq. 3) then
!
!     ON NE TRAITE PAS LE SEGMENT DE VALEURS MARQUE X D
!
                                        goto 210
                                    end if
!
!     LE SEGMENT DE VALEURS EST MARQUE X D, IL FAUT D'ABORD L'ECRIRE
!
                                    iaddi(1) = iszon(jiszon+ibiadd-1+2*k-1)
                                    iaddi(2) = iszon(jiszon+ibiadd-1+2*k)
                                    call jxecro(ic, iadmoc, iaddi, lsv, j, k)
                                    iszon(jiszon+ibiadd-1+2*k-1) = iaddi(1)
                                    iszon(jiszon+ibiadd-1+2*k) = iaddi(2)
                                end if
                                lgs = iszon(jiszon+iadmoc-4)-iadmoc+4
                                mcdyn = mcdyn-lgs
                                mldyn = mldyn+lgs
                                call hpdeallc(iadyoc, nbfree)
                                ltot = ltot+lgs
!       WRITE(6,*) ' OC ',NOM32,' OBJET ',K,' LG =',LSV,LGS,LTOT
                                iszon(jiszon+ibiadm-1+2*k-1) = 0
                                iszon(jiszon+ibiadm-1+2*k) = 0
                                if (lmin .gt. 0) then
                                    if (ltot .ge. lmin) then
                                        lgio(1) = lgio(1)+1024*longbl(ic)*lois* &
                                                  (nbacce(2*ic-1)-nbioav(1))
                                        lgio(2) = lgio(2)+1024*longbl(ic)*lois* &
                                                  (nbacce(2*ic)-nbioav(2))
                                        goto 300
                                    end if
                                end if
                            end if
                        end if
210                     continue
                    end do
                end if
                goto 205
!          ELSE IF ( NOM32(25:32) .EQ. ' ' ) THEN
            else
                if (iadyn .ne. 0) then
                    idm = iadmi-4
                    isd = iszon(jiszon+idm+3)/isstat
                    isf = iszon(jiszon+iszon(jiszon+idm)-4)/isstat
                    if (isd .eq. 1) then
!
!     LE SEGMENT DE VALEURS EST MARQUE X A OU X D, ON PEUT LE LIBERER
!
                        ltypi = ltyp(jltyp(ic)+j)
                        lsv = lono(jlono(ic)+j)*ltypi
                        if (isf .eq. 4) then
                            if (imode .eq. 2 .or. imode .eq. 3) then
!
!     ON NE TRAITE PAS LE SEGMENT DE VALEURS MARQUE X D
!
                                goto 205
                            end if
!
!     LE SEGMENT DE VALEURS EST MARQUE X D, IL FAUT D'ABORD L'ECRIRE
!
                            iaddi(1) = iadd(jiadd(ic)+2*j-1)
                            iaddi(2) = iadd(jiadd(ic)+2*j)
                            call jxecro(ic, iadmi, iaddi, lsv, 0, j)
                            iadd(jiadd(ic)+2*j-1) = iaddi(1)
                            iadd(jiadd(ic)+2*j) = iaddi(2)
                        end if
                        lgs = iszon(jiszon+iadmi-4)-iadmi+4
                        mcdyn = mcdyn-lgs
                        mldyn = mldyn+lgs
                        call hpdeallc(iadyn, nbfree)
                        ltot = ltot+lgs
!            write(6,*) ' OS ',NOM32,' lg =',LSV,LGS,LTOT
                        iadm(jiadm(ic)+2*j-1) = 0
                        iadm(jiadm(ic)+2*j) = 0
                        if (lmin .gt. 0) then
                            if (ltot .ge. lmin) then
                                lgio(1) = lgio(1)+1024*longbl(ic)*lois*(nbacce(2*ic-1)-nbioav(1))
                                lgio(2) = lgio(2)+1024*longbl(ic)*lois*(nbacce(2*ic)-nbioav(2))
                                goto 300
                            end if
                        end if
                    end if
                end if
            end if
205         continue
        end do
!
        lgio(1) = lgio(1)+1024*longbl(ic)*lois*(nbacce(2*ic-1)-nbioav(1))
        lgio(2) = lgio(2)+1024*longbl(ic)*lois*(nbacce(2*ic)-nbioav(2))
200     continue
    end do
300 continue
    mxltot = mxltot+(ltot*lois)/(1024*1024)
!
    if (lmin .ne. -2) then
!
!   ON TESTE LA VALEUR DE VMSIZE PAR RAPPORT AUX ALLOCATIONS JEVEUX ET
!   AU RELIQUAT UNIQUEMENT SI LMIN DIFFERENT DE  -2
!
        call utgtme(5, nomk, valp, iret)
        if (is_enabled(LOGLEVEL_MEM, DEBUG)) then
            do i = 1, 5
                write (6, *) "DEBUG: ", nomk(i), valp(i)
            end do
        end if
        v0 = valp(5)
!   Si VmSize == 0, on ne fait rien
        if (valp(3) .gt. 0 .and. valp(3)-(valp(1)+valp(2)) .gt. 0) then
!
!   ON AJUSTE LA VALEUR DU RELIQUAT ET LA LIMITE DES ALLOCATONS JEVEUX
!
            call utptme('RLQ_MEM ', valp(3)-valp(1), iret)
            vx(1) = valp(4)-(valp(3)-valp(1))
            if (vx(1) .gt. 0) then
                call jermxd(vx(1)*1024*1024, iret)
                if (iret .eq. 0) then
                    call utgtme(5, nomk, valp, iret2)
!
!   ON IMPRIME UN MESSAGE D'INFORMATION SI LA VALEUR DE LIMIT_JV VARIE
!   DE PLUS DE 10 POUR CENT
!
                    if (abs(valp(5)-v0) .gt. v0*0.1d0) then
                        vx(1) = valp(2)
                        vx(2) = valp(5)
                        vx(3) = valp(5)-v0
                        call utmess('I', 'JEVEUX1_74', nr=3, valr=vx)
                    end if
                end if
            end if
        end if
    end if
!
    call uttcpu('CPU.MEMD.1', 'FIN', ' ')
    indiq_jjldyn = 0
!
end subroutine
