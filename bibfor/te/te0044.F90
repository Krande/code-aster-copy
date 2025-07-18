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
subroutine te0044(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/ptenci.h"
#include "asterfort/ptenpo.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2pgl.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utmess.h"
#include "asterfort/utppgl.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/vecmap.h"
!
    character(len=*) :: option, nomte
!     CALCUL DE L'ENERGIE DE DEFORMATION, ET CINETIQUE
!     ------------------------------------------------------------------
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'EPOT_ELEM' : CALCUL DE L'ENERGIE DE DEFORMATION
!        'ECIN_ELEM' : CALCUL DE L'ENERGIE CINETIQUE
! IN  NOMTE  : K16 : NOM DU TYPE D'ELEMENT DISCRET :
!         MECA_DIS_T_N
!         MECA_DIS_T_L
!         MECA_DIS_TR_N
!         MECA_DIS_TR_L
!
!      CHARACTER*32 JEXNUM,JEXNOM,JEXR8,JEXATR
    real(kind=8) :: r8bid, ul(12), pgl(3, 3), klc(12, 12), mat(144)
    integer(kind=8) :: infodi, nbterm, nno, nc, ndim, itype, neq, kanl, irep, iiff
    integer(kind=8) :: i, lorie, ldepl, lvite, jende, ldis, jfreq, ibid, iret
    character(len=3) :: stopz
    character(len=8) :: k8bid
!
!     ------------------------------------------------------------------
!
    infodi = 1
!
!     ON VERIFIE QUE LES CARACTERISTIQUES ONT ETE AFFECTEES
!     LE CODE DU DISCRET
    call infdis('CODE', ibid, r8bid, nomte)
!     LE CODE STOKE DANS LA CARTE
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    end if
!
    if (option .eq. 'EPOT_ELEM') then
!        DISCRET DE TYPE RAIDEUR
        call infdis('DISK', infodi, r8bid, k8bid)
        if (infodi .eq. 0) then
            call utmess('A+', 'DISCRETS_27', sk=nomte)
            call infdis('DUMP', ibid, r8bid, 'A+')
        end if
        call infdis('SYMK', infodi, r8bid, k8bid)
    else if (option .eq. 'ECIN_ELEM') then
!        DISCRET DE TYPE MASSE
        call infdis('DISM', infodi, r8bid, k8bid)
        if (infodi .eq. 0) then
            call utmess('A+', 'DISCRETS_26', sk=nomte)
            call infdis('DUMP', ibid, r8bid, 'A+')
        end if
        call infdis('SYMM', infodi, r8bid, k8bid)
    else
        call utmess('F', 'ELEMENTS2_47', sk=option)
    end if
!
! --- INFORMATIONS SUR LES DISCRETS :
!        NBTERM   = NOMBRE DE COEFFICIENTS DANS K
!        NNO      = NOMBRE DE NOEUDS
!        NC       = NOMBRE DE COMPOSANTE PAR NOEUD
!        NDIM     = DIMENSION DE L'ELEMENT
!        ITYPE    = TYPE DE L'ELEMENT
    call infted(nomte, infodi, nbterm, nno, nc, &
                ndim, itype)
    neq = nno*nc
!
!     TYPE DE LA MATRICE DE MASSE
    kanl = 0
!     --- MATRICE DE ROTATION PGL ---
    call jevech('PCAORIE', 'L', lorie)
    call matrot(zr(lorie), pgl)
!     --- RECUPERATION DES DEPLACEMENTS OU DES VITESSES ET PASSAGE
!     --- AU REPERE LOCAL
    if (option .ne. 'ECIN_ELEM') then
        call jevech('PDEPLAR', 'L', ldepl)
        if (ndim .eq. 3) then
            call utpvgl(nno, nc, pgl, zr(ldepl), ul)
        else if (ndim .eq. 2) then
            call ut2vgl(nno, nc, pgl, zr(ldepl), ul)
        end if
    else
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iret, iad=lvite)
! IRET NE PEUT VALOIR QUE 0 (TOUT EST OK) OU 2 (CHAMP NON FOURNI)
        if (iret .eq. 0) then
            if (ndim .eq. 3) then
                call utpvgl(nno, nc, pgl, zr(lvite), ul)
            else if (ndim .eq. 2) then
                call ut2vgl(nno, nc, pgl, zr(lvite), ul)
            end if
        else
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=ldepl)
            if (iret .eq. 0) then
                if (ndim .eq. 3) then
                    call utpvgl(nno, nc, pgl, zr(ldepl), ul)
                else if (ndim .eq. 2) then
                    call ut2vgl(nno, nc, pgl, zr(ldepl), ul)
                end if
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
    end if
!
    if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jende)
!
!        --- MATRICE DE RIGIDITE ---
        call jevech('PCADISK', 'L', ldis)
!        --- GLOBAL VERS LOCAL ? ---
!        --- IREP EQ 1 : MATRICE EN REPERE GLOBAL
!        --- IREP NE 1 : MATRICE EN REPERE LOCAL
        call infdis('REPK', irep, r8bid, k8bid)
        if (irep .eq. 1) then
            if (ndim .eq. 3) then
                if (infodi .eq. 1) then
                    call utpsgl(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call utppgl(nno, nc, pgl, zr(ldis), mat)
                end if
            else if (ndim .eq. 2) then
                if (infodi .eq. 1) then
                    call ut2mgl(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call ut2pgl(nno, nc, pgl, zr(ldis), mat)
                end if
            end if
        else
            do i = 1, nbterm
                mat(i) = zr(ldis+i-1)
            end do
        end if
!
!        ---- MATRICE RIGIDITE LIGNE > MATRICE RIGIDITE CARRE
        if (infodi .eq. 1) then
            call vecma(mat, nbterm, klc, neq)
        else if (infodi .eq. 2) then
            call vecmap(mat, nbterm, klc, neq)
        end if
!
!        --- ENERGIE DE DEFORMATION ---
        iiff = 1
        call ptenpo(neq, ul, klc, zr(jende), itype, &
                    iiff)
!
    else if (option .eq. 'ECIN_ELEM') then
        call jevech('PENERCR', 'E', jende)
!
!        --- MATRICE DE MASSE ---
        call jevech('PCADISM', 'L', ldis)
!        --- GLOBAL VERS LOCAL ? ---
!        --- IREP EQ 1 : MATRICE EN REPERE GLOBAL
!        --- IREP NE 1 : MATRICE EN REPERE LOCAL
        call infdis('REPM', irep, r8bid, k8bid)
        if (irep .eq. 1) then
            if (ndim .eq. 3) then
                if (infodi .eq. 1) then
                    call utpsgl(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call utppgl(nno, nc, pgl, zr(ldis), mat)
                end if
            else if (ndim .eq. 2) then
                if (infodi .eq. 1) then
                    call ut2mgl(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call ut2pgl(nno, nc, pgl, zr(ldis), mat)
                end if
            end if
        else
            do i = 1, nbterm
                mat(i) = zr(ldis+i-1)
            end do
        end if
!
!        ---- MATRICE RIGIDITE LIGNE > MATRICE RIGIDITE CARRE
        if (infodi .eq. 1) then
            call vecma(mat, nbterm, klc, neq)
        else if (infodi .eq. 2) then
            call vecmap(mat, nbterm, klc, neq)
        end if
!
!        --- FREQUENCE ---
        call jevech('POMEGA2', 'L', jfreq)
!
!        --- ENERGIE CINETIQUE  ---
        iiff = 1
        call ptenci(neq, ul, klc, zr(jfreq), zr(jende), &
                    itype, kanl, iiff)
!
    else
        call utmess('F', 'ELEMENTS2_47', sk=option)
    end if
end subroutine
