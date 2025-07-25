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
subroutine cnvesl(lischa, typres, neq, nompar, valpar, &
                  cnvass)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/fointc.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liscpp.h"
#include "asterfort/lisico.h"
#include "asterfort/lislch.h"
#include "asterfort/lislco.h"
#include "asterfort/lislnf.h"
#include "asterfort/lisltc.h"
#include "asterfort/lisltf.h"
#include "asterfort/lisnbg.h"
#include "asterfort/lisnnb.h"
    character(len=19) :: lischa
    character(len=19) :: cnvass
    character(len=1) :: typres
    character(len=8) :: nompar
    integer(kind=8) :: neq
    real(kind=8) :: valpar, tval(1)
!
! ----------------------------------------------------------------------
!
! CALCUL CONTRIBUTION SECOND MEMBRE SI VECT_ASSE
!
! ----------------------------------------------------------------------
!
!
! IN  LISCHA : SD LISTE DES CHARGES
! IN  TYPRES : TYPE DU CHAM_NO RESULTANT 'C'
! IN  NOMPAR : NOM DU PARAMETRE
! IN  VALPAR : VALEUR DU PARAMETRE
! IN  NEQ    : NOMBRE D'EQUATIONS DU SYSTEME
! OUT CNVASS : NOM DU CHAMP
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ichar, nbchar
    integer(kind=8) :: nbveas, nbveag, nbtot, iret, ieq
    integer(kind=8) :: genrec
    integer(kind=8) :: jvale
    character(len=16) :: typfct
    character(len=8) :: nomfct, charge, typech
    real(kind=8) :: valre, valim
    complex(kind=8) :: calpha, calp
    real(kind=8) :: phase, omega
    integer(kind=8) :: npuis
    aster_logical :: lveas, lveag
    character(len=24) :: chamno
    complex(kind=8), pointer :: resu(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ASSERT(typres .eq. 'C')
    omega = r8depi()*valpar
    call jeveuo(cnvass(1:19)//'.VALE', 'E', vc=resu)
    do ieq = 1, neq
        resu(ieq) = dcmplx(0.d0, 0.d0)
    end do
!
! --- NOMBRE DE CHARGEMENTS
!
    call lisnnb(lischa, nbchar)
!
! --- NOMBRE DE CHARGES DE TYPE VECT_ASSE
!
    nbveas = lisnbg(lischa, 'VECT_ASSE')
    nbveag = lisnbg(lischa, 'VECT_ASSE_GENE')
    nbtot = nbveas+nbveag
    if (nbtot .eq. 0) goto 999
!
! --- BOUCLE SUR LES CHARGES
!
    do ichar = 1, nbchar
!
! ----- CODE DU GENRE DE LA CHARGE
!
        call lislco(lischa, ichar, genrec)
!
! ----- FONCTION MULTIPLICATRICE
!
        call lislnf(lischa, ichar, nomfct)
        call lisltf(lischa, ichar, typfct)
!
! ----- MULTIPLICATEUR COMPLEXE
!
        call liscpp(lischa, ichar, phase, npuis)
!
        lveas = lisico('VECT_ASSE', genrec)
        lveag = lisico('VECT_ASSE_GENE', genrec)
        if (lveas .or. lveag) then
            call lislch(lischa, ichar, charge)
            call lisltc(lischa, ichar, typech)
            chamno = charge
!
            valre = 1.d0
            valim = 0.d0
            if (nomfct .ne. ' ') then
                tval(1) = valpar
                if (typfct(7:10) .eq. 'REEL') then
                    call fointe('F', nomfct, 1, nompar, tval, &
                                valre, iret)
                    valim = 0.d0
                else if (typfct(7:10) .eq. 'COMP') then
                    call fointc('F', nomfct, 1, nompar, tval, &
                                valre, valim, iret)
                else
                    ASSERT(.false.)
                end if
            end if
            calp = dcmplx(valre, valim)
            calpha = calp*exp(dcmplx(0.d0, phase*r8dgrd()))
            if (npuis .ne. 0) then
                calpha = calpha*omega**npuis
            end if
            call jeveuo(chamno(1:19)//'.VALE', 'L', jvale)
            if (typech .eq. 'COMP') then
                do ieq = 1, neq
                    resu(ieq) = resu(ieq)+calpha*zc(jvale-1+ieq)
                end do
            else
                do ieq = 1, neq
                    resu(ieq) = resu(ieq)+calpha*zr(jvale-1+ieq)
                end do
            end if
        end if
    end do
!
999 continue
!
    call jedema()
end subroutine
