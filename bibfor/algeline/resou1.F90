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

subroutine resou1(matass, matpre, solveu, chcine, nsecm, &
                  chsecm, chsolu, base, rsolu, csolu, &
                  criter, prepos, istop, iret)
!
!-----------------------------------------------------------------------
! BUT : RESOUDRE UN SYSTEME LINEAIRE D'EQUATIONS (REEL OU COMPLEXE)
!-----------------------------------------------------------------------
!
! ARGUMENTS :
!------------
!
! remarques : on peut appeler resou1 de 2 facons
!   1) avec nsecm = 0 + chsecm, chsolu, base
!   2) avec nsecm > 0 + rsolu (ou csolu) + (chsecm=chsolu=' ')
!
! IN/JXIN  K19 MATASS : MATR_ASSE PREMIER MEMBRE DU SYSTEME LINEAIRE
! IN/JXIN  K19 MATPRE : MATR_ASSE DE PRECONDITIONNEMENT
!                       POUR SOLVEUR ITERATIF GCPC (OU ' ' SINON)
! IN/JXIN  K19 SOLVEU : SD_SOLVEUR (OU ' ')
!                       SI SOLVEU=' ' ON PREND LE SOLVEUR DE MATASS
! IN/JXIN  K*  CHCINE : CHAMP ASSOCIE AUX CHARGES CINEMATIQUES (OU ' ')
! IN       I   NSECM  : / 0 => ON UTILISE CHSECM, CHSOLU, BASE
!                       / N => ON UTILISE RSOLU (OU CSOLU)
!                         N : NOMBRE DE SECONDS MEMBRES
! IN/JXIN  K*  CHSECM : CHAMP SECOND MEMBRE DU SYSTEME LINEAIRE
! IN/JXOUT K*  CHSOLU : CHAMP SOLUTION DU SYSTEME LINEAIRE
! IN       K*  BASE   : BASE SUR LAQUELLE ON CREE CHSOLU
! IN/OUT   R   RSOLU  : TABLEAU (*,NSECM)
!           EN ENTREE : VECTEUR DE REELS CONTENANT LES SECONDS MEMBRES
!           EN SORTIE : VECTEUR DE REELS CONTENANT LES SOLUTIONS
! IN/OUT   C   CSOLU  : TABLEAU (*,NSECM)
!                       IDEM RSOLU POUR LES COMPLEXES
! IN/JXOUT K*  CRITER : SD_CRITER (CRITERES DE CONVERGENCE)
!                       POUR SOLVEUR ITERATIF GCPC (OU ' ' SINON)
! IN       L   PREPOS : / .TRUE.  => ON FAIT LES PRE ET POST-TRAITEMENTS
!                                    DU SMB ET DE LA SOLUTION
!                       / .FALSE. => ON NE FAIT AUCUN TRAITEMENT
!                                    (EN MODAL PAR EXEMPLE)
! IN       I   ISTOP  : COMPORTEMENT EN CAS D'ERREUR (CE PARAMETRE N'A
!                       D'UTILITE QUE POUR UN SOLVEUR ITERATIF)
!                       / 0     : ON S'ARRETE EN <F>
!                       / 2     : ON CONTINUE SANS MESSAGE D'ERREUR
!                       / -9999 : ON PREND LA VALEUR DEFINIE DANS LA
!                                 SD_SOLVEUR POUR STOP_SINGULIER
! OUT      I   IRET   : CODE RETOUR
!                       / 0 : OK (PAR DEFAUT POUR SOLVEURS DIRECTS)
!                       / 1 : ECHEC (NOMBRE MAX. D'ITERATIONS ATTEINT)
!-----------------------------------------------------------------------
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/cheksd.h"
#include "asterfort/amumph.h"
#include "asterfort/apetsc.h"
#include "asterfort/assert.h"
#include "asterfort/dbgobj.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jacopo.h"
#include "asterfort/jedbg2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdscr.h"
#include "asterfort/resgra.h"
#include "asterfort/resldl.h"
#include "asterfort/utmess.h"
#include "asterfort/vtdefs.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMatrix.h"
    character(len=*) :: matass, matpre, solveu, chcine
    integer(kind=8) :: nsecm
    character(len=*) :: chsecm, chsolu, base
    real(kind=8) :: rsolu(*)
    complex(kind=8) :: csolu(*)
    character(len=*) :: criter
    aster_logical :: prepos
    integer(kind=8) :: istop, iret
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, ifm, niv
    character(len=3) :: kmpic, type, typ1, khpc
    character(len=19) :: matr19, mpre19, solv19, cine19
    character(len=19) :: secm19, csol19, crit19
    character(len=24) :: metres
!
    integer(kind=8) :: idbgav, neq, neq1, ier, niter, lmat, jvals
    integer(kind=8) :: jtrav, jval2, imd, istopz
    real(kind=8) :: epsi
    complex(kind=8) :: cbid
    aster_logical :: dbg, l_parallel_matrix
    character(len=1), parameter :: ftype(2) = ['R', 'C']
    character(len=24), pointer :: slvk(:) => null()
    character(len=24), pointer :: refa(:) => null()
    real(kind=8), pointer :: slvr(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
! ----------------------------------------------------------------------
    dbg = .false.
!
    call jemarq()
    call infniv(ifm, niv)
    call jedbg2(idbgav, 0)
!
    matr19 = matass
    mpre19 = matpre
    solv19 = solveu
    cine19 = chcine
    secm19 = chsecm
    csol19 = chsolu
    crit19 = criter
!
    ASSERT(matr19 .ne. ' ')
    call dismoi('MPI_COMPLET', matr19, 'MATR_ASSE', repk=kmpic)
!
    if (solv19 .eq. ' ') call dismoi('SOLVEUR', matr19, 'MATR_ASSE', repk=solv19)
    call jeveuo(solv19//'.SLVK', 'L', vk24=slvk)
    call jeveuo(solv19//'.SLVR', 'L', vr=slvr)
    call jeveuo(solv19//'.SLVI', 'L', vi=slvi)
    metres = slvk(1)
    ASSERT(metres .ne. ' ')
    if (kmpic .eq. 'NON') then
        ASSERT(metres .eq. 'MUMPS' .or. metres .eq. 'PETSC')
    end if
!
!     VERIFICATIONS ET INITIALISATIONS
    ASSERT((istop .eq. 0) .or. (istop .eq. 2) .or. (istop .eq. -9999))
    if (istop .eq. -9999) then
        istopz = slvi(8)
    else
        istopz = istop
    end if
    iret = 0
!
    call mtdscr(matr19)
    call jeveuo(matr19//'.&INT', 'L', lmat)
    neq = zi(lmat+2)
    type = ftype(zi(lmat+3))
!
    ASSERT(nsecm .ge. 0)
    call jeveuo(matr19//'.REFA', 'L', vk24=refa)
    if (refa(11) .eq. 'MATR_DISTR') then
        imd = 1
    else
        imd = 0
    end if
    if (nsecm .eq. 0) then
        ASSERT(secm19 .ne. ' ')
        ASSERT(csol19 .ne. ' ')
        if (csol19 .ne. secm19) then
            call detrsd('CHAMP_GD', csol19)
            call vtdefs(csol19, secm19, base, ' ')
        end if
!
        call jelira(secm19//'.VALE', 'LONMAX', neq1)
        call jelira(secm19//'.VALE', 'TYPE', cval=typ1)
        if ((neq1 .ne. neq) .and. (imd .eq. 0)) then
            call utmess('F', 'FACTOR_67')
        end if
        if (typ1 .ne. type) then
            call utmess('F', 'FACTOR_68')
        end if
!
        call jeveuo(secm19//'.VALE', 'L', jval2)
        if (imd .eq. 0) then
            call wkvect('&&RESOU1.TRAV', 'V V '//type, neq, jtrav)
            call jacopo(neq, type, jval2, jtrav)
        else
            call wkvect('&&RESOU1.TRAV', 'V V '//type, neq1, jtrav)
            call jacopo(neq1, type, jval2, jtrav)
        end if
!
    else
        ASSERT(secm19 .eq. ' ')
        ASSERT(csol19 .eq. ' ')
    end if
!
    if (cine19 .ne. ' ') then
        call jelira(cine19//'.VALE', 'TYPE', cval=typ1)
        ASSERT(typ1 .eq. type)
    end if
!
!
!
!
!
    if (dbg) then
        call cheksd(matr19, 'SD_MATR_ASSE', ier)
        if (nsecm .eq. 0) call dbgobj(secm19//'.VALE', 'OUI', 6, '&&RESOU1 2ND MEMBRE')
        call dbgobj(cine19//'.VALE', 'OUI', 6, '&&RESOU1 CINE19')
        call dbgobj(matr19//'.VALM', 'OUI', 6, '&&RESOU1 MATR.VALM')
        call dbgobj(matr19//'.VALF', 'OUI', 6, '&&RESOU1 MATR.VALF')
        call dbgobj(matr19//'.CONL', 'OUI', 6, '&&RESOU1 MATR.CONL')
        call dbgobj(matr19//'.CCVA', 'OUI', 6, '&&RESOU1 MATR.CCVA')
    end if
!
!
!
!
    if (metres .eq. 'LDLT' .or. metres .eq. 'MULT_FRONT') then
!     ----------------------------------------------------
        l_parallel_matrix = isParallelMatrix(matr19)
        if (l_parallel_matrix) then
            call utmess('F', 'FACTOR_93')
        end if
        if (nsecm .gt. 0) then
            call resldl(solv19, matr19, cine19, nsecm, rsolu, &
                        csolu, prepos)
        else
            if (type .eq. 'R') then
                call resldl(solv19, matr19, cine19, 1, zr(jtrav), &
                            [cbid], prepos)
            else
                call resldl(solv19, matr19, cine19, 1, [0.d0], &
                            zc(jtrav), prepos)
            end if
        end if
!
!
!
    else if (metres .eq. 'MUMPS') then
!     ----------------------------------------------------
        if (nsecm .gt. 0) then
            call amumph('RESOUD', solv19, matr19, rsolu, csolu, &
                        cine19, nsecm, iret, prepos)
        else
            if (type .eq. 'R') then
                call amumph('RESOUD', solv19, matr19, zr(jtrav), [cbid], &
                            cine19, 1, iret, prepos)
            else
                call amumph('RESOUD', solv19, matr19, [0.d0], zc(jtrav), &
                            cine19, 1, iret, prepos)
            end if
        end if
        ASSERT(iret .eq. 0)
!
!
!
    else if (metres .eq. 'GCPC') then
!     ----------------------------------
        niter = slvi(2)
        epsi = slvr(2)
        ASSERT(type .eq. 'R')
        if (nsecm .gt. 0) then
            call resgra(matr19, mpre19, cine19, niter, epsi, &
                        crit19, nsecm, rsolu, solv19, istopz, &
                        iret)
        else
            call resgra(matr19, mpre19, cine19, niter, epsi, &
                        crit19, 1, zr(jtrav), solv19, istopz, &
                        iret)
        end if
!
!
!
    else if (metres .eq. 'PETSC') then
!     ----------------------------------
        ASSERT(type .eq. 'R')
        if (nsecm .gt. 0) then
            call apetsc('RESOUD', solv19, matr19, rsolu, cine19, &
                        nsecm, istopz, iret)
        else
            call apetsc('RESOUD', solv19, matr19, zr(jtrav), cine19, &
                        1, istopz, iret)
        end if
!
    else
        call utmess('F', 'ALGELINE3_44', sk=metres)
    end if
!
!
!     -- RECOPIE DANS LE CHAMP SOLUTION S'IL Y A LIEU :
    if (nsecm .eq. 0) then
        call jeveuo(csol19//'.VALE', 'E', jvals)
        if (imd .eq. 0) then
            call jacopo(neq, type, jtrav, jvals)
        else
            call jacopo(neq1, type, jtrav, jvals)
        end if
    end if
    call jedetr('&&RESOU1.TRAV')
!
!
!
!
    if (dbg .and. (nsecm .eq. 0)) call dbgobj(csol19//'.VALE', 'OUI', 6, '&&RESOU1 SOLU')
!
!
    call jedbg2(ibid, idbgav)
    call jedema()
end subroutine
