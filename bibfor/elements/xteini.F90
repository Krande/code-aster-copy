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

subroutine xteini(nomte, nfh, nfe, singu, ddlc, &
                  nnom, ddls, nddl, ddlm, nfiss, &
                  contac)
    implicit none
!
#include "jeveux.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/ismali.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
!
    character(len=16) :: nomte
    integer(kind=8) :: nfh, nfe, singu, ddlc, nnom, ddls, nddl, ddlm
    integer(kind=8) :: nfiss, contac
!
! person_in_charge: samuel.geniaut at edf.fr
!
!          BUT : INITIALISER LES DIMENSIONS DES DDL DANS UN TE
!                POUR LES ELEMENTS X-FEM
!
!
! IN   NOMTE  : NOM DU TYPE ELEMENT
! OUT  NFH    : NOMBRE DE FONCTIONS HEAVISIDES
! OUT  NFE    : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! OUT  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! OUT  DDLC   : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! OUT  NNOM   : NB DE NOEUDS MILIEU SERVANT À PORTER DES DDL DE CONTACT
! OUT  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! OUT  NDDL   : NOMBRE DE DDL TOTAL DE L'ÉLÉMENT
! OUT  DDLM   : NOMBRE DE DDL A CHAQUE NOEUD MILIEU
! OUT  NFISS  : NOMBRE DE FISSURES
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, ier, nnos
    integer(kind=8) :: ddld, iadzi, iazk24, jtab(7), iret
    character(len=8) :: elrefp, enr, typma
!
! ----------------------------------------------------------------------
!
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos)
!
! --- INITIALISATIONS
!
    nfh = 0
    nfe = 0
    singu = 0
    ddlc = 0
    ddlm = 0
    nnom = 0
    nfiss = 1
    contac = 0
!
    call teattr('S', 'XFEM', enr, ier, typel=nomte)
!
! --- DDL ENRICHISSEMENT : HEAVYSIDE, ENRICHIS (FOND)
!
    if (enr(1:2) .eq. 'XH') then
        nfh = 1
        if (enr(1:3) .eq. 'XH2') nfh = 2
        if (enr(1:3) .eq. 'XH3') nfh = 3
        if (enr(1:3) .eq. 'XH4') nfh = 4
!       NOMBRE DE FISSURES
        call tecach('NOO', 'PLST', 'L', iret, nval=7, &
                    itab=jtab)
        nfiss = jtab(7)
    end if
!
    if (enr(1:2) .eq. 'XT' .or. enr(3:3) .eq. 'T') then
        nfe = 1
        singu = 1
    end if
!
! --- DDL DE CONTACT
!
    if (enr(1:4) .eq. 'XHC3') then
        ddlc = 3*ndim
    else if (enr(1:3) .eq. 'XHC' .or. enr(1:3) .eq. 'XTC' .or. enr(1:4) .eq. 'XHTC') then
        ddlc = ndim
    end if
    if (enr(1:4) .eq. 'XH2C') ddlc = 2*ndim
    if (enr(1:4) .eq. 'XH3C') ddlc = 3*ndim
    if (enr(1:4) .eq. 'XH4C') ddlc = 4*ndim
!
! --- NOMBRE DE DDL DE DEPLACEMENT
!
    ddld = ndim*(1+nfh+nfe)
!
! --- NOMBRE DE DDL AUX NOEUDS SOMMETS
!
    ddls = ddld+ddlc
!
! --- NOMBRE DE DDL AUX NOEUDS MILIEUX
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
    if (ier .eq. 0) then
        if (ismali(typma)) then
            if (enr(1:4) .eq. 'XHC3') then
                contac = 2
            else
                contac = 1
            end if
            ddlm = 0
        else
            contac = 3
            ddlm = ddld
        end if
    else
        if (.not. iselli(elrefp)) ddlm = ddld
    end if
!
! --- NB DE NOEUDS MILIEUX
!
    nnom = nno-nnos
!
! --- NOMBRE DE DDLS (DEPL+CONTACT) SUR L'ELEMENT
!
    nddl = (nnos*ddls)+(nnom*ddlm)
!
end subroutine
