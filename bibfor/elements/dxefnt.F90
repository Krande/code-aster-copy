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

subroutine dxefnt(nomte, pgl, sigt)
    implicit none
#include "jeveux.h"
#include "asterfort/dxmath.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    real(kind=8) :: pgl(3, 3), sigt(*)
    character(len=16) :: nomte
! --- EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX NOEUDS
! --- POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G DUS :
! ---  .A UN CHAMP DE TEMPERATURES SUR LE PLAN MOYEN DONNANT
! ---        DES EFFORTS DE MEMBRANE
! ---  .A UN GRADIENT DE TEMPERATURES DANS L'EPAISSEUR DE LA COQUE
!     ------------------------------------------------------------------
!     IN  NOMTE        : NOM DU TYPE D'ELEMENT
!     IN  XYZL(3,NNO)  : COORDONNEES DES CONNECTIVITES DE L'ELEMENT
!                        DANS LE REPERE LOCAL DE L'ELEMENT
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     OUT SIGT(1)      : EFFORTS  GENERALISES D'ORIGINE THERMIQUE
!                        AUX NOEUDS
    integer(kind=8) :: icodre(56)
    character(len=10) :: phenom
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3)
    real(kind=8) :: t2iu(4), t2ui(4), t1ve(9)
    real(kind=8) :: tsup(4), tinf(4), tmoy(4), rbid
    integer(kind=8) :: multic, nbcou, jcou, imoy
!     ------------------------------------------------------------------
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: indith, ino, iret1, iret2, iret3, iret4
    integer(kind=8) :: jcara, jmate, nno
    real(kind=8) :: coe1, coe2, epais, somire, tref, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    iret1 = 0
    iret2 = 0
    iret3 = 0
    iret4 = 0
!
    call r8inir(32, 0.d0, sigt, 1)
!
!     -- S'IL N'Y A PAS DE TEMPERATURE, IL N'Y A RIEN A CALCULER :
    call rcvarc(' ', 'TEMP', '+', 'NOEU', 1, 1, rbid, iret4)
    if (iret4 .ne. 0) goto 30
!
!
!
    call jevech('PNBSP_I', 'L', jcou)
    nbcou = zi(jcou)
    call rcvarc(' ', 'TEMP', 'REF', 'NOEU', 1, 1, tref, iret1)
!
!
    if (nomte .eq. 'MEDKTR3 ' .or. nomte .eq. 'MEDSTR3 ' &
        .or. nomte .eq. 'MEDKTG3 ' .or. nomte .eq. 'MET3TR3 ') then
        nno = 3
    else if (nomte .eq. 'MEDKQU4 ' .or. nomte .eq. 'MEDKQG4 ' .or. &
             nomte .eq. 'MEDSQU4 ' .or. nomte .eq. 'MEQ4QU4 ') then
        nno = 4
    else
        call utmess('F', 'ELEMENTS_14', sk=nomte)
    end if
!
!===============================================================
!          -- RECUPERATION DE LA TEMPERATURE  AUX NOEUDS
! COQUE MULTI-COUCHE.
! ON RECUPERE LA TEMPERATURE INFERIEURE, SUPERIEURE ET DANS LA FIBRE
! MOYENNE
    imoy = (3*nbcou+1)/2
    do ino = 1, nno
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, imoy, tmoy(ino), iret4)
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, nbcou*3, tsup(ino), iret3)
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, 1, tinf(ino), iret2)
    end do
!
    call jevech('PMATERC', 'L', jmate)
    call rccoma(zi(jmate), 'ELAS', 1, phenom, icodre(1))
!
    if ((phenom .eq. 'ELAS') .or. (phenom .eq. 'ELAS_COQUE') &
        .or. (phenom .eq. 'ELAS_COQMU') .or. (phenom .eq. 'ELAS_GLRC')) then
!
! --- RECUPERATION DE L'EPAISSEUR DE LA COQUE
!     --------------------------
!
        call jevech('PCACOQU', 'L', jcara)
        epais = zr(jcara)
!
! --- CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE,
! --- MEMBRANE-FLEXION, CISAILLEMENT, CISAILLEMENT INVERSE
!     ----------------------------------------------------
        call dxmath('NOEU', epais, df, dm, dmf, pgl, multic, indith, t2iu, t2ui, t1ve, nno)
        if (indith .ne. -1) then
!
            somire = iret2+iret3+iret4
            if (somire .eq. 0) then
                if (iret1 .eq. 1) then
                    call utmess('F', 'COMPOR5_43')
                else
!
! --- BOUCLE SUR LES NOEUDS
!     ---------------------
                    do ino = 1, nno
!
!  --      LES COEFFICIENTS SUIVANTS RESULTENT DE L'HYPOTHESE SELON
!  --      LAQUELLE LA TEMPERATURE EST PARABOLIQUE DANS L'EPAISSEUR.
!  --      ON NE PREJUGE EN RIEN DE LA NATURE DU MATERIAU.
!  --      CETTE INFORMATION EST CONTENUE DANS LES MATRICES QUI
!  --      SONT LES RESULTATS DE LA ROUTINE DXMATH.
!          ----------------------------------------
                        coe1 = (tsup(ino)+tinf(ino)+4.d0*tmoy(ino))/6.d0-tref
                        coe2 = (tsup(ino)-tinf(ino))/epais
!
                        sigt(1+8*(ino-1)) = coe1*(dm(1, 1)+dm(1, 2))+coe2*(dmf(1, 1)+dmf(1, 2))
                        sigt(2+8*(ino-1)) = coe1*(dm(2, 1)+dm(2, 2))+coe2*(dmf(2, 1)+dmf(2, 2))
                        sigt(3+8*(ino-1)) = coe1*(dm(3, 1)+dm(3, 2))+coe2*(dmf(3, 1)+dmf(3, 2))
                        sigt(4+8*(ino-1)) = coe2*(df(1, 1)+df(1, 2))+coe1*(dmf(1, 1)+dmf(1, 2))
                        sigt(5+8*(ino-1)) = coe2*(df(2, 1)+df(2, 2))+coe1*(dmf(2, 1)+dmf(2, 2))
                        sigt(6+8*(ino-1)) = coe2*(df(3, 1)+df(3, 2))+coe1*(dmf(3, 1)+dmf(3, 2))
                    end do
                end if
            end if
        end if
    end if
30  continue
end subroutine
