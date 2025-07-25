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

subroutine jefini(cond, close)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/asabrt.h"
#include "asterc/xfini.h"
#include "asterfort/assert.h"
#include "asterfort/enlird.h"
#include "asterfort/iunifi.h"
#include "asterfort/jeimpr.h"
#include "asterfort/jelibf.h"
#include "asterfort/jjlidy.h"
#include "asterfort/ulclos.h"
#include "asterfort/utgtme.h"

    character(len=*) :: cond
    aster_logical, intent(in), optional :: close
! ======================================================================
!
    integer(kind=8) :: i, n
!-----------------------------------------------------------------------
    parameter(n=5)
!
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    integer(kind=8) :: nbfic
    common/iparje/nbfic
!     ------------------------------------------------------------------
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
!
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    integer(kind=8) :: icdyn, mxltot
    common/xdynje/icdyn, mxltot
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
!     ==================================================================
    integer(kind=8) :: vali(7), info, ifm, ires, iret
    aster_logical :: close_base
    character(len=8) :: kcond, staou, k8tab(3)
    character(len=24) :: ladate
    real(kind=8) :: rval(3)
!     ------------------------------------------------------------------
!
    if (present(close)) then
        close_base = close
    else
        close_base = ASTER_FALSE
    end if
!

    kcond = cond(1:min(len(cond), len(kcond)))
    ASSERT(kcond .eq. 'NORMAL  ' .or. kcond .eq. 'ERREUR  ' .or. kcond .ne. 'TEST    ')
    if (kcond .eq. 'NORMAL  ' .or. kcond .eq. 'TEST    ') then
        staou = '        '
    else
        staou = 'SAUVE   '
    end if
!
!     EVALUATION DE LA CONSOMMATION MEMOIRE
!
    k8tab(1) = 'VMPEAK'
    k8tab(2) = 'MEM_TOTA'
    call utgtme(2, k8tab, rval, iret)
!
!     -------------  EDITION DES REPERTOIRES ---------------------------
    if (kcond .eq. 'TEST    ') then
        do i = 1, nbfic
            if (classe(i:i) .ne. ' ') then
                call jeimpr(6, classe(i:i), '     JEFINI     '//kcond)
            end if
        end do
    end if
!     -------------  LIBERATION FICHIER --------------------------------
    if (kcond .ne. 'ERREUR  ') then
        if (close_base) then
            info = 0
            do i = 1, nbfic
                if (classe(i:i) .ne. ' ') then
                    call jelibf(staou, classe(i:i), info)
                end if
            end do
!       -----------  DESALLOCATION GESTION DES MARQUES -----------------
            call jjlidy(kdesma(2), kdesma(1))
            call jjlidy(kposma(2), kposma(1))
            kdesma(1) = 0
            kdesma(2) = 0
            kposma(1) = 0
            kposma(2) = 0
!
        end if
    else
        call asabrt(6)
    end if
!
!     --- IMPRESSION DES CONSOMMATIONS MEMOIRE ---
!
    k8tab(1) = 'CMXU_JV'
    k8tab(2) = 'CMAX_JV'
    k8tab(3) = 'VMPEAK'
    call utgtme(3, k8tab, rval, iret)
    ifm = iunifi('MESSAGE')
    ires = iunifi('RESULTAT')
!
    if (ires .gt. 0) then
        write (ires, *) ' '
        write (ires, '(2A,F11.2,A)')&
     &        ' <I> <FIN> MEMOIRE JEVEUX MINIMALE REQUISE POUR ',&
     &        'L''EXECUTION :                ', rval(1), ' Mo'
        write (ires, '(2A,F11.2,A)')&
     &        ' <I> <FIN> MEMOIRE JEVEUX OPTIMALE REQUISE POUR ',&
     &        'L''EXECUTION :                ', rval(2), ' Mo'
        if (rval(3) .gt. 0) then
            write (ires, '(2A,F11.2,A)')&
     &        ' <I> <FIN> MAXIMUM DE MEMOIRE UTILISEE PAR LE PROCESSUS'&
     &        , ' LORS DE L''EXECUTION :', rval(3), ' Mo'
        end if
    end if
!
    if (kcond .ne. 'TEST    ') then
        if (ifm .gt. 0) then
            write (ifm, *) ' '
            write (ifm, *) '<I>       FERMETURE DES BASES EFFECTUEE'
            if (ldyn .eq. 1) then
                vali(1) = nint(mxdyn/(1024*1024))
                vali(2) = liszon*lois/(1024*1024)
                vali(3) = nbdyn
                vali(4) = nbfree
                vali(5) = nint(mldyn/(1024*1024))
                vali(6) = nint(lgio(1)/(1024*1024))
                vali(7) = nint(lgio(2)/(1024*1024))
                write (ifm, *) ' '
                write (ifm, *) '  STATISTIQUES CONCERNANT L'''&
     &                 //'ALLOCATION DYNAMIQUE :'
                write (ifm, *) '    TAILLE CUMULEE MAXIMUM            :',&
     &                 vali(1), ' Mo.'
                write (ifm, *) '    TAILLE CUMULEE LIBEREE            :',&
     &                 vali(5), ' Mo.'
                write (ifm, *) '    NOMBRE TOTAL D''ALLOCATIONS        :',&
     &                 vali(3)
                write (ifm, *) '    NOMBRE TOTAL DE LIBERATIONS       :',&
     &                 vali(4)
                write (ifm, *) '    APPELS AU MECANISME DE LIBERATION :',&
     &                 icdyn
                write (ifm, *) '    TAILLE MEMOIRE CUMULEE RECUPEREE  :',&
     &                 mxltot, ' Mo.'
                write (ifm, *) '    VOLUME DES LECTURES               :',&
     &                  vali(6), ' Mo.'
                write (ifm, *) '    VOLUME DES ECRITURES              :',&
     &                  vali(7), ' Mo.'
                write (ifm, *) ' '
            end if
            write (ifm, '(A,F11.2,A)')&
     &       '   MEMOIRE JEVEUX MINIMALE REQUISE POUR L''EXECUTION :',&
     &       rval(1), ' Mo'
            write (ifm, '(A)') '     - IMPOSE DE NOMBREUX ACCES DISQUE'
            write (ifm, '(A)') '     - RALENTIT LA VITESSE D''EXECUTION'
            write (ifm, '(A,F11.2,A)')&
     &       '   MEMOIRE JEVEUX OPTIMALE REQUISE POUR L''EXECUTION :',&
     &       rval(2), ' Mo'
            write (ifm, '(A)') '     - LIMITE LES ACCES DISQUE'
            write (ifm, '(A)') '     - AMELIORE LA VITESSE D''EXECUTION'
            if (rval(3) .gt. 0) then
                write (ifm, '(A,F11.2,A)')&
     &       '   MAXIMUM DE MEMOIRE UTILISEE PAR LE PROCESSUS     :',&
     &       rval(3), ' Mo'
                write (ifm, '(A)') '     - COMPREND LA MEMOIRE CONSOMMEE PAR '//&
     &       ' JEVEUX, '
                write (ifm, '(A)') '       LE SUPERVISEUR PYTHON, '//&
     &       'LES LIBRAIRIES EXTERNES'
            end if
            write (ifm, *) ' '
!
            call enlird(ladate)
            write (ifm, *) '<I>       FIN D''EXECUTION LE : '//ladate
!
        end if
!
!       --- ON FERME TOUT ---
        call ulclos()
!
        call xfini(19)
    end if
end subroutine
