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
subroutine vppcom(lcomod, icom1, icom2, resui, resur, &
                  resuk, nbpari, nbparr, nbpark, mxresf, &
                  vectr, nconv, neq, typres)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_split_comm.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/comatr.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/maxint.h"
#include "asterfort/somint.h"
#include "asterfort/vecint.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    aster_logical :: lcomod
    integer(kind=8) :: icom1, icom2, nbpari, nbparr, nbpark, nconv, neq, mxresf
    integer(kind=8) :: resui(*)
    real(kind=8) :: vectr(*), resur(*)
    character(len=16) :: typres
    character(len=*) :: resuk(*)
!     COMMUNICATION DES VECTEURS PROPRES ET SD ASSOCIEES POUR
!     CALC_MODES SUR PLUSIEURS SOUS-BANDES //
!     IN  LCOMOD      :  LOG  : CALC_MODES SUR PLUSIEURS SOUS-BANDES // OR NOT
!     IN  ICOM1/ICOM2 :  IN   : PARAMETRES // ASSOCIES A LCOMOD=.TRUE.
!     IN/OUT RESUI/RESUR/RESUK: VECTEURS I/R/K : VECTEURS ASSOCIES A LA
!                               SD_MODE_MECA A PRODUIRE.
!     IN NBPARI/NBPARR/NBPARK : IN : NBRE DE COLONNES DES CES SDS.
!     IN MXRESF       :  IN   : NBRE DE LIGNES UTILISEES POUR CES SDS.
!     IN/OUT VECTR    : R     : VECTEURS PROPRES REELS
!     IN/OUT NCONV    :  IN   : NBRE DE MODES PROPRES CONVERGES DU PROC
!                (EN IN), NBRE TOTAL DE MODES CV (EN OUT).
!     IN NEQ          :  IN   : TAILLE DU PROBLEME
!     IN TYPRES       :  IN   : TYPE DE RESULTATS (DYNAMIQUE OU FLAMB)
!     ----------------------------------------------------------------
!
!     --- VARIABLES LOCALES
    mpi_int :: rangl, rangll, mpicow, mpicou, mpico0, l1, l2
    integer(kind=8) :: nconvl, nconvg, nconvm
    integer(kind=8) :: ibid, izero, i, idecal, j, i8, jlcom, jlbuff, jlbufs, ifm
    integer(kind=8) :: niv
    real(kind=8) :: rbid
    complex(kind=8) :: cbid
    character(len=24) :: klcom, k24buf, k24bus, k24b
    blas_int :: b_incx, b_incy, b_n
!      LOGICAL      LCPU
!
    call jemarq()
    call infniv(ifm, niv)
!     --- POUR MESURER LE TEMPS CONSOMMEE DS LES ETAPES DE OP0045
!     --- IL FAUT AUSSI DECOMMENTER LES APPELS A SYSTEM_CLOCK (UTILES
!     --- MAIS REFUSES PAR L'AGLA).
!      LCPU=.TRUE.
!      LCPU=.FALSE.
    izero = 0
    k24buf = '&&OP0045.BUFFMPI'
    k24bus = '&&OP0045.BUFFMPI.SAVE'
    klcom = '&&OP0045.LCOMOD1'
    if (lcomod) then
        call asmpi_comm('GET_WORLD', mpicow)
        call asmpi_comm('GET', mpicou)
!       --- ON EST CENSE FONCTIONNER JUSQUE LA EN COM LOCAL
!       --- (POUR SOLVEUR LINEAIRE DU SOLVEUR MODAL)
        if (mpicow .eq. mpicou) then
            ASSERT(.false.)
        end if
        call asmpi_info(comm=mpicou, rank=rangl)
!       ----------------------------------------------------------------
!       --- STEP 0: COMM AU SEIN DU COM_WORLD
!       ----------------------------------------------------------------
!       --- ON REMET LE MPICOW POUR LES COMS DE NCONVG ET DE KLCOM.
!       --- ON GARDE LE COM LOCAL MPICOU POUR LE BCAST ULTERIEUR
!       --- ON PREPARE UN COM LOCAL RESERVE A TOUS LES PROCS DE
!       --- RANGL=0: MPICO0
!       --- DANS CE NEW COM MPICO0, LE RANG EST NOTE RANGLL. IL DOIT
!       --- ETRE IDENTIQUE A ICOM1-1.
        call asmpi_comm('SET', mpicow)
        call asmpi_barrier()
        if (rangl .eq. 0) then
            l1 = 1
        else
            l1 = 2
        end if
        call asmpi_split_comm(mpicow, l1, to_mpi_int(icom1), 'procs0', mpico0)
        ASSERT(mpicow .ne. mpico0)
        call asmpi_info(mpico0, rangll, l2)
        if ((l2 .ne. icom2) .and. (rangl .eq. 0)) then
            ASSERT(.false.)
        end if
        if ((rangll .ne. (icom1-1)) .and. (rangl .eq. 0)) then
            ASSERT(.false.)
        end if
!
!       --- VECTEUR DES DECALAGES (NCONV) PAR PROC POUR VECTEURS PROPRES
!       --- CALCUL ET AFFECTATION DU NBRE TOTAL DE MODES CONVERGES
        call wkvect(klcom, 'V V I', icom2, jlcom)
        call vecint(icom2, izero, zi(jlcom))
!       --- NCONVL: NBRE DE MODES CONVERGES DU CALCUL ASSOCIE AU PROC
!       --- ZI(JLCOM+ICOM1-1)=NCONVL SI RANGL=0, 0 SINON
!       --- NCONV=NCONVG: NBRE DE MODES CONVERGES TOTAL
!       --- NCONVM=MAX(NCONVL) (POUR DIMENSIONNER BUFFER DE COM)
        nconvl = nconv
        if (rangl .eq. 0) zi(jlcom+icom1-1) = nconvl
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=icom2, vi=zi(jlcom))
        nconvg = somint(icom2, zi(jlcom))
        if (nconvg .le. 0) then
            ASSERT(.false.)
        end if
        nconv = nconvg
        nconvm = maxint(icom2, zi(jlcom))
        if (nconvm .le. 0) then
            ASSERT(.false.)
        end if
!
!       ----------------------------------------------------------------
!       --- STEP 1: COMM ENTRE LES MAITRES DE CHAQUE SOUS-BANDES
!       ----------------------------------------------------------------
!       --- ON AFFECTE LE COMCO0 POUR COMMUNIQUER ENTRE SOUS-BANDES
!       --- LES VECTEURS PROPRES
!       --- PUIS COMM DES &&OP0045.RESU_R
!       --- COMM DES &&OP0045.RESU_I ET K
        call asmpi_barrier()
        call asmpi_comm('SET', mpico0)
        if (rangl .eq. 0) then
!         --- 2 BUFFERS: K24BUF POUR LE BCAST, LE K24BUS POUR SAUVE
!         --- GARDER LES VECTEURS PROPRES EN ATTENDANT LEUR COMM AUX
!         --- AUTRES PROCESSUS DE MPICO0.
            call wkvect(k24buf, 'V V R', nconvm*neq, jlbuff)
            call wkvect(k24bus, 'V V R', nconvm*neq, jlbufs)
            b_n = to_blas_int(nconvl*neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vectr, b_incx, zr(jlbufs), b_incy)
            do i = 1, icom2
                idecal = 0
                do j = 1, i-1
                    idecal = idecal+zi(jlcom+j-1)
                end do
                if (idecal .lt. 0) then
                    ASSERT(.false.)
                end if
                i8 = neq*zi(jlcom+i-1)
                if (i .eq. icom1) then
                    b_n = to_blas_int(i8)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(jlbufs), b_incx, zr(jlbuff), b_incy)
                end if
                call asmpi_comm_vect('BCAST', 'R', nbval=i8, bcrank=i-1, vr=zr(jlbuff))
                b_n = to_blas_int(i8)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jlbuff), b_incx, vectr(1+idecal*neq), b_incy)
            end do
            call jedetr(k24bus)
            call jedetr(k24buf)
!
            ibid = 0
            rbid = 0.d0
            cbid = dcmplx(0.d0)
            call comatr('T', 'R', icom2, icom1-1, zi(jlcom), &
                        0, 0, [ibid], mxresf, nbparr, &
                        resur, 0, 0, [cbid])
            call comatr('T', 'I', icom2, icom1-1, zi(jlcom), &
                        mxresf, nbpari, resui, 0, 0, &
                        [rbid], 0, 0, [cbid])
        end if
!
!       ----------------------------------------------------------------
!       --- STEP 2: BARRIERE SUR LE COM_WORLD AU CAS OU
!       ----------------------------------------------------------------
!       -- ON AFFECTE LE COMCOW POUR QUE TOUS LES PROCS S'ATTENDENT
        call asmpi_comm('SET', mpicow)
        call asmpi_barrier()
!        IF (LCPU) THEN
!          CALL SYSTEM_CLOCK(IETFIN)
!          RETFIN=REAL(IETFIN-IETDEB)/REAL(IETRAT)
!          WRITE(IFM,*)'<VPPCOM> COUT COM 1  VP + ZR/ZI: ',RETFIN
!          CALL SYSTEM_CLOCK(IETDEB,IETRAT,IETMAX)
!        ENDIF
!
!       ----------------------------------------------------------------
!       --- STEP 3: COM AU SEIN DE CHAQUE SOUS-BANDES
!       ----------------------------------------------------------------
!       --- ON AFFECTE LE COMCOU POUR COMMUNIQUER AU SEIN
!       --- D'UNE SOUS-BANDE.
        call asmpi_comm('SET', mpicou)
!
!       --- POUR GAGNER DU TEMPS, ON ENVOIE VP PAR VP
!       --- ON COMMUNIQUE LES VECTEURS PROPRES EN SOUS-PAQUETS DE REELS
!       --- DE TAILLE =< SIZBMPI POUR EVITER LES PBS DE CONTENTIONS
!       --- MEMOIRE ET LES LIMITES DES ENTIERS COURTS MPI.
        do i = 1, nconvg
            call asmpi_comm_vect('BCAST', 'R', nbval=neq, bcrank=0, vr=vectr(1+(i-1)*neq))
        end do
!       --- ON COMMUNIQUE LES PETITS OBJETS SUIVANTS
        call asmpi_comm_vect('BCAST', 'R', nbval=nbparr*mxresf, bcrank=0, vr=resur)
        call asmpi_comm_vect('BCAST', 'I', nbval=nbpari*mxresf, bcrank=0, vi=resui)
!
!       ----------------------------------------------------------------
!       --- STEP 4: BARRIERE SUR LE COM_WORLD AU CAS OU
!       ----------------------------------------------------------------
!       --- ON AFFECTE LE COMCOW POUR QUE TOUS LES PROCS S'ATTENDENT
        call asmpi_comm('SET', mpicow)
        call asmpi_barrier()
!        IF (LCPU) THEN
!          CALL SYSTEM_CLOCK(IETFIN)
!          RETFIN=REAL(IETFIN-IETDEB)/REAL(IETRAT)
!          WRITE(IFM,*)'<VPPCOM> COUT COM 2  VP + ZR/ZI: ',RETFIN
!        ENDIF
!
!       --- ON REMET LE COM LOCAL AU CAS OU
        call asmpi_comm('SET', mpicou)
!
!       ----------------------------------------------------------------
!       --- STEP 5: INIT DIVERSES.
!       ----------------------------------------------------------------
        if (typres(1:9) .ne. 'DYNAMIQUE') then
            do i = 1, nconvg
                resui(i) = i
            end do
        end if
        do i = 1, nbpark
            j = 1+(i-1)*mxresf
            k24b = resuk(j)
            resuk(j:j-1+nconvg) = k24b
        end do
!       ----------------------------------------------------------------
!       --- STEP 6: MENAGE.
!       ----------------------------------------------------------------
        call jedetr(klcom)
        call asmpi_comm('FREE', mpico0)
    end if
!
    call jedema()
!
end subroutine
