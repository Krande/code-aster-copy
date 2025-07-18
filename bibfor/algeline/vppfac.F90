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
subroutine vppfac(lmasse, masgen, vect, neq, nbvect, &
                  mxvect, masmod, facpar)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterc/r8vide.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/pteddl.h"
#include "asterfort/wkvect.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsvpar.h"
#include "blas/ddot.h"
#include "asterfort/gettco.h"
!
    integer(kind=8) :: lmasse, neq, nbvect, mxvect
    real(kind=8) :: masgen(*), vect(neq, *)
    real(kind=8) :: masmod(mxvect, *), facpar(mxvect, *)
!     CALCUL DES PARAMETRES MODAUX :
!            FACTEUR DE PARTICIPATION ET MASSE MODALE UNITAIRE
!     ------------------------------------------------------------------
! IN  LMASSE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE MASSE
!     ------------------------------------------------------------------
!     PRESUME L'EXECUTION PREALABLE DE VPPGEN : CALCUL DES PARAMETRES
!     MODAUX.
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: lddl, laux1, laux2, iddl, ia, ieq, ivect, mxddl, iadpar(1), l1, ibid
    parameter(mxddl=6)
    character(len=8) :: nomddl(mxddl), basemo, k8b
    character(len=14) :: nume
    character(len=16) :: nompar(3), typmas, typbas
    character(len=19) :: masse
    character(len=24) :: posddl, vecau1, vecau2
    real(kind=8) :: rmin, rmax, raux, rval
    aster_logical :: gene
    character(len=24), pointer :: refn(:) => null()
    real(kind=8) :: rundef
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
    data nomddl/'DX      ', 'DY      ', 'DZ      ',&
     &              'DRX     ', 'DRY     ', 'DRZ     '/
    data nompar/'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ'/
!
!     ------------------------------------------------------------------
    data posddl/'&&VPPFAC.POSITION.DDL'/
    data vecau1/'&&VPPFAC.VECTEUR.AUX1'/
    data vecau2/'&&VPPFAC.VECTEUR.AUX2'/
!
    rundef = r8vide()
!
!     ------------------------------------------------------------------
!     ----------------- CREATION DE VECTEURS DE TRAVAIL ----------------
!     ------------------------------------------------------------------
!
    call wkvect(vecau1, 'V V R', neq, laux1)
    call wkvect(vecau2, 'V V R', neq, laux2)
    call wkvect(posddl, 'V V I', neq*mxddl, lddl)
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------ RECUPERATION DES POSITIONS DU VECTEUR EXCITATION : --------
!     ---- VECTEUR DE 1 DANS LA DIRECTION iddl DANS LE CAS PHYSIQUE ----
!     --- FACTEURS DE PARTICIPATION DES MODES DANS LE CAS GENERALISE ---
!     ------------------------------------------------------------------
!
    call jemarq()
    masse = zk24(zi(lmasse+1)) (1:19)
    call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=nume)
    call gettco(masse, typmas)
!
    rmin = 100.d0*r8miem()
    rmax = sqrt(r8maem())
!
! DETERMINATION DU CAS : BASE PHYSIQUE OU BASE GENERALISEE
    gene = .false.
    if (typmas(1:14) .eq. 'MATR_ASSE_GENE') then
! SI MATR_ASSE_GENE : BASE GENERALISEE
        call jeveuo(nume(1:14)//'.NUME.REFN', 'L', vk24=refn)
        basemo = refn(1) (1:8)
! SAUF SI NUME.REFN POINTE VERS UN MODELE_GENE ET NON VERS UNE BASE
! ALORS CAS DE LA SSD, TRAITE COMME UN MODE_MECA CLASSIQUE
        call gettco(basemo, typbas)
        if (typbas(1:14) .eq. 'MODE_MECA') then
            gene = .true.
        end if
    end if
    do iddl = 1, 3
        if (gene) then
            do ieq = 1, neq
                call rsvpar(basemo, 1, nompar(iddl), ibid, rundef, &
                            k8b, l1)
                if (l1 .eq. 100) then
                    zr(laux1+ieq-1) = 0.D0
                else
                    call rsadpa(basemo, 'L', 1, nompar(iddl), ieq, &
                                0, tjv=iadpar)
                    zr(laux1+ieq-1) = zr(iadpar(1))
                end if
! SECURITE SI ON EST PASSE PAR DES MODES HETERODOXES AVEC FACTEURS DE PARTICIPATIONS HERETIQUES
            end do
        else
            call pteddl('NUME_DDL', nume, mxddl, nomddl, neq, &
                        tabl_equa=zi(lddl))
            ia = (iddl-1)*neq
            do ieq = 1, neq
                zr(laux1+ieq-1) = zi(lddl+ia+ieq-1)
            end do
        end if
!
!     ------------------------------------------------------------------
!     ----------- CALCUL DE  FREQ * MASSE * UNITAIRE_DIRECTION ---------
!     ------------------------------------------------------------------
        call mrmult('ZERO', lmasse, zr(laux1), zr(laux2), 1, &
                    .false._1)
        do ivect = 1, nbvect
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            rval = ddot(b_n, vect(1, ivect), b_incx, zr(laux2), b_incy)
            raux = masgen(ivect)
            if ((abs(raux) .lt. rmin) .or. (abs(rval) .gt. rmax)) then
                masmod(ivect, iddl) = rmax
                facpar(ivect, iddl) = rmax
            else
                raux = rval/raux
                masmod(ivect, iddl) = rval*raux
                facpar(ivect, iddl) = raux
            end if
        end do
    end do
!
!     ------------------------------------------------------------------
!     ----------------- DESTRUCTION DES VECTEURS DE TRAVAIL ------------
!     ------------------------------------------------------------------
!
    call jedetr(posddl)
    call jedetr(vecau1)
    call jedetr(vecau2)
!
    call jedema()
end subroutine
