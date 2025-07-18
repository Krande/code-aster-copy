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

subroutine xhmini(nomte, nfh, ddld, ddlm, ddlp, nfiss, ddlc, contac)
!
    implicit none
!
#   include "asterfort/elref1.h"
#   include "asterfort/elrefe_info.h"
#   include "asterfort/teattr.h"
#   include "asterfort/tecach.h"
#   include "jeveux.h"
!
    character(len=16) :: nomte
    integer(kind=8) :: nfh, ddlm, ddlp, ddlc
    integer(kind=8) :: nfiss, contac
!
! person_in_charge: daniele.colombo at ifpen.fr
!
!          BUT : INITIALISER LES DIMENSIONS DES DDL DANS UN TE
!                POUR LES ELEMENTS HM-XFEM
!
!
! IN   NOMTE  : NOM DU TYPE ELEMENT
! OUT  NFH    : NOMBRE DE FONCTIONS HEAVISIDES
! OUT  DDLD   : NOMBRE DE DDL (DEPL) A CHAQUE NOEUD SOMMET
! OUT  DDLM   : NOMBRE DE DDL (DEPL) A CHAQUE NOEUD MILIEU
! OUT  DDLP   : NOMBRE DE DDL (PRES) A CHAQUE NOEUD SOMMET
! OUT  DDLC   : NOMBRE DE DDL POUR LE CONTACT (SOMMET UNIQUEMENT)
! OUT  CONTAC : =3 POUR CONTACT P2P1 (PAS DE P1P1 EN HM-XFEM)
! OUT  NFISS  : NOMBRE DE FISSURES
! OUT  NFH    : NOMBRE DE DDL HEAVISIDE PAR NOEUD
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nnop, ier, nnops
    integer(kind=8) :: ddld, jtab(7), iret
    character(len=8) :: elrefp, enr
!
! ----------------------------------------------------------------------
!
    call elref1(elrefp)
    call elrefe_info(elrefe=elrefp, fami='RIGI', ndim=ndim, nno=nnop, nnos=nnops)
!
! --- INITIALISATIONS
!
    nfh = 0
    ddlm = 0
    ddld = 0
    ddlp = 0
    ddlc = 0
    contac = 0
    nfiss = 1
!
    call teattr('S', 'XFEM', enr, ier, typel=nomte)
!
! --- DDL ENRICHISSEMENT : HEAVYSIDE
!
    if (enr(1:2) .eq. 'XH') then
        nfh = 1
        if (enr(1:3) .eq. 'XH2') nfh = 2
        if (enr(1:3) .eq. 'XH3') nfh = 3
!       NOMBRE DE FISSURES
        call tecach('NOO', 'PLST', 'L', iret, nval=7, &
                    itab=jtab)
        nfiss = jtab(7)
    end if

! --- NOMBRE DE DDL POUR LE CONTACT: PRE_FLU, LAG_FLI, LAG_FLS
!     LAG1_HM ET LAG2_HM ET INDICATION DU CONTACT: P2P1
!
    if (enr(3:3) .eq. 'C' .or. enr(4:4) .eq. 'C') then
! --- CONTACT MORTAR
        if (enr(4:4) .eq. '3' .or. enr(5:5) .eq. '3') then
            ddlc = (1+1+1+3*ndim)*nfh
            contac = 2
! --- CONTACT COHESIF
        else
            ddlc = (1+1+1+ndim)*nfh
            contac = 3
        end if
    end if
!
! --- NOMBRE DE DDL AUX NOEUDS SOMMETS (MECANIQUES)
!
    ddld = ndim*(1+nfh)
!
! --- NOMBRE DE DDL AUX NOEUDS MILIEUX (MECANIQUES)
!
    ddlm = ddld
!
! --- NOMBRE DE DDL AUX NOEUDS SOMMETS (HYDRAULIQUES)
!
    ddlp = 1+nfh
end subroutine
