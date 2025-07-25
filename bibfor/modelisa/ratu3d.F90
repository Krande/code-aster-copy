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

subroutine ratu3d(iprno, lonlis, klisno, noepou, noma, &
                  ligrel, mod, cara, numddl, &
                  lisrel, coorig, sectio)
    implicit none
#include "jeveux.h"
#include "asterfort/afretu.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/memare.h"
#include "asterfort/mesomm.h"
#include "asterfort/raorfi.h"
#include "asterfort/reajre.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: lonlis, iprno(*)
    character(len=8) :: klisno(lonlis), noepou, noma, cara, mod
    character(len=14) :: numddl
    character(len=19) :: ligrel, lisrel
    real(kind=8) :: coorig(3), rayon
!     RACCORD 3D_TUYAU PAR DES RELATIONS LINEAIRES
!     ECRITURE DES RELATIONS SUR LES DDLS DE FOURIER
!
    integer(kind=8) :: nbcmp, nbmode
    parameter(nbmode=3, nbcmp=6*(nbmode-1))
    character(len=8) :: nocmp(nbcmp), lpain(5), lpaout(6), nomddl(4)
    character(len=24) :: lchin(5), lchout(6), valech
    real(kind=8) :: coef(4), eg1(3), eg2(3), eg3(3), sectio
    integer(kind=8) :: imod, info, ifm
    integer(kind=8) :: nbcoef, idec
    real(kind=8), pointer :: rayon_raccord(:) => null()
!
    call jemarq()
    call infniv(ifm, info)
!
!
!     APPEL DE L OPTION CARA_SECT_POU3R AFIN DE CALCULER LE RAYON
!     MOYEN DE LA SECTION ANNULAIRE DE LA PARTIE 3D
!     L APPEL A MESOMM DONNE L INTEGRALE SURFACIQUE DU RAYON
!
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = '&&RAPO3D.CAORIGE'
    lpaout(1) = 'PRAYONM'
    lchout(1) = '&&RATU3D.PRAYO'
!
    call calcul('S', 'CARA_SECT_POU3R', ligrel, 2, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    AS_ALLOCATE(vr=rayon_raccord, size=1)
    call mesomm(lchout(1), 1, vr=rayon_raccord)
    rayon = rayon_raccord(1)/sectio
!
!
!     CREATION D'UNE CARTE CONTENANT LE POINT P ORIGINE DE PHI
!
    call raorfi(noma, ligrel, noepou, cara, coorig, &
                eg1, eg2, eg3, '&&RATU3D', rayon)
!
! --- DETERMINATION DE 3 LISTES  DE VECTEURS PAR ELEMENT PRENANT
! --- LEURS VALEURS AUX NOEUDS DES ELEMENTS. CES VALEURS SONT :
!     VECT_COEF_UM
! --- SOMME/S_ELEMENT(COS(M.PHI).P(1,J).NI)DS
! --- SOMME/S_ELEMENT(SIN(M.PHI).P(1,J).NI)DS
!     VECT_COEF_VM
! --- SOMME/S_ELEMENT(COS(M.PHI).P(2,J).NI)DS
! --- SOMME/S_ELEMENT(SIN(M.PHI).P(2,J).NI)DS
!     VECT_COEF_WM
! --- SOMME/S_ELEMENT(COS(M.PHI).P(3,J).NI)DS
! --- SOMME/S_ELEMENT(SIN(M.PHI).P(3,J).NI)DS
!
! --- OU P EST LA MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
! --- (E1,E2,E3) DEFINI SUR LE BORD ORIENTE DE LA SURFACE
!     ------------------------------
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PORIGIN'
    lchin(2) = '&&RAPO3D.CAORIGE'
    lpain(3) = 'PCAORIE'
    lchin(3) = '&&RATU3D.CAXE_TUY'
    lpain(4) = 'PORIGFI'
    lchin(4) = '&&RATU3D.CAORIFI'
    lpain(5) = 'PNUMMOD'
    lchin(5) = '&&RATU3D.NUME_MODE'
    lpaout(1) = 'PVECTU1'
    lchout(1) = '&&RATU3D.COSF_UM'
    lpaout(2) = 'PVECTU2'
    lchout(2) = '&&RATU3D.COSF_VM'
    lpaout(3) = 'PVECTU3'
    lchout(3) = '&&RATU3D.COSF_WM'
    lpaout(4) = 'PVECTU4'
    lchout(4) = '&&RATU3D.SINF_UM'
    lpaout(5) = 'PVECTU5'
    lchout(5) = '&&RATU3D.SINF_VM'
    lpaout(6) = 'PVECTU6'
    lchout(6) = '&&RATU3D.SINF_WM'
!
! --- CREATION DES .RERR DES VECTEURS EN SORTIE DE CALCUL
!
    call memare('V', '&&RATU3D', mod, 'CHAR_MECA')
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL WO
!
    imod = 0
    call mecact('V', lchin(5), 'LIGREL', ligrel, 'NUMMOD', &
                ncmp=1, nomcmp='NUM', si=imod)
    call calcul('S', 'CARA_SECT_POUT5', ligrel, 5, lchin, &
                lpain, 6, lchout, lpaout, 'V', &
                'OUI')
    call jedetr('&&RATU3D           .RELR')
    call reajre('&&RATU3D', lchout(3), 'V')
    call assvec('V', 'CH_DEPL_3', 1, '&&RATU3D           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_3          .VALE'
    nbcoef = 1
    idec = 0
    nomddl(1) = 'WO'
    coef(1) = -1.d0*sectio
    call afretu(iprno, lonlis, klisno, noepou, noma, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)
!
!   RELATIONS ENTRE LES NOEUDS DE COQUE ET LE NOEUD POUTRE
!   DDL UIM, UOM, VIM, VOM, WIM, WOM, M VARIANT DE 2 A NBMODE
!
    nocmp(1) = 'UI2'
    nocmp(2) = 'VI2'
    nocmp(3) = 'WI2'
    nocmp(4) = 'UO2'
    nocmp(5) = 'VO2'
    nocmp(6) = 'WO2'
    nocmp(7) = 'UI3'
    nocmp(8) = 'VI3'
    nocmp(9) = 'WI3'
    nocmp(10) = 'UO3'
    nocmp(11) = 'VO3'
    nocmp(12) = 'WO3'
!
    imod = 1
    if (info .eq. 2) then
        write (ifm, *) 'RELATIONS SUR LE MODE ', imod
    end if
    call mecact('V', lchin(5), 'LIGREL', ligrel, 'NUMMOD', &
                ncmp=1, nomcmp='NUM', si=imod)
    call calcul('S', 'CARA_SECT_POUT5', ligrel, 5, lchin, &
                lpain, 6, lchout, lpaout, 'V', &
                'OUI')
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL WIM
!     OU SI IMOD=1 LE DDL DZ DANS REPERE LOCAL DU TUYAU ET WI1
!     IDEC=0 SIGNIFIE QUE ON UTILISE LES TERMES EN COS(M*PHI)
!
    call jedetr('&&RATU3D           .RELR')
    call reajre('&&RATU3D', lchout(3), 'V')
    call assvec('V', 'CH_DEPL_3', 1, '&&RATU3D           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_3          .VALE'
    idec = 0
!
    nbcoef = 4
    nomddl(1) = 'DX'
    nomddl(2) = 'DY'
    nomddl(3) = 'DZ'
    nomddl(4) = 'WI1'
    coef(1) = eg3(1)*sectio/2.d0
    coef(2) = eg3(2)*sectio/2.d0
    coef(3) = eg3(3)*sectio/2.d0
    coef(4) = -sectio/2.d0

    call afretu(iprno, lonlis, klisno, noepou, noma, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL WOM
!     OU SI IMOD=1 LE DDL DY DANS REPERE LOCAL DU TUYAU ET WO1
!
    call jedetr('&&RATU3D           .RELR')
    call reajre('&&RATU3D', lchout(6), 'V')
    call assvec('V', 'CH_DEPL_6', 1, '&&RATU3D           .RELR', [1.d0], numddl)
    valech = 'CH_DEPL_6          .VALE'
    idec = 0
!
    nbcoef = 4
    nomddl(1) = 'DX'
    nomddl(2) = 'DY'
    nomddl(3) = 'DZ'
    nomddl(4) = 'WO1'
    coef(1) = eg2(1)*sectio/2.d0
    coef(2) = eg2(2)*sectio/2.d0
    coef(3) = eg2(3)*sectio/2.d0
    coef(4) = -sectio/2.d0
!
    call afretu(iprno, lonlis, klisno, noepou, noma, &
                valech, nbcoef, idec, coef, nomddl, &
                lisrel)
!
    do imod = 2, nbmode
        if (info .eq. 2) then
            write (ifm, *) 'RELATIONS SUR LE MODE ', imod
        end if
        call mecact('V', lchin(5), 'LIGREL', ligrel, 'NUMMOD', &
                    ncmp=1, nomcmp='NUM', si=imod)
        call calcul('S', 'CARA_SECT_POUT5', ligrel, 5, lchin, &
                    lpain, 6, lchout, lpaout, 'V', &
                    'OUI')
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL UIM
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(1), 'V')
        call assvec('V', 'CH_DEPL_1', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_1          .VALE'
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+1)
        coef(1) = -sectio/2.d0
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL UOM
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(4), 'V')
        call assvec('V', 'CH_DEPL_4', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_4          .VALE'
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+4)
        coef(1) = -sectio/2.d0
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL VOM
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(2), 'V')
        call assvec('V', 'CH_DEPL_2', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_2          .VALE'
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+5)
        coef(1) = -sectio/2.d0
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL VIM
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(5), 'V')
        call assvec('V', 'CH_DEPL_5', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_5          .VALE'
        idec = 0
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+2)
        coef(1) = -sectio/2.d0
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL WIM
!     OU SI IMOD=1 LE DDL DZ DANS REPERE LOCAL DU TUYAU ET WI1
!     IDEC=0 SIGNIFIE QUE ON UTILISE LES TERMES EN COS(M*PHI)
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(3), 'V')
        call assvec('V', 'CH_DEPL_3', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_3          .VALE'
        idec = 0
!
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+3)
        coef(1) = -sectio/2.d0
!
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!         ENDIF
!
!     RELATIONS ENTRE LES NOEUDS DE SURFACE ET LE NOEUD POUTRE DDL WOM
!     OU SI IMOD=1 LE DDL DY DANS REPERE LOCAL DU TUYAU ET WO1
!
        call jedetr('&&RATU3D           .RELR')
        call reajre('&&RATU3D', lchout(6), 'V')
        call assvec('V', 'CH_DEPL_6', 1, '&&RATU3D           .RELR', [1.d0], &
                    numddl)
        valech = 'CH_DEPL_6          .VALE'
        idec = 0
!
        nbcoef = 1
        nomddl(1) = nocmp(6*(imod-2)+6)
        coef(1) = -sectio/2.d0
!
        call afretu(iprno, lonlis, klisno, noepou, noma, &
                    valech, nbcoef, idec, coef, nomddl, &
                    lisrel)
!
!     FIN DE LA BOUCLE SUR LES MODES
!
    end do
!
! --- DESTRUCTION DES OBJETS DE TRAVAIL
!
    AS_DEALLOCATE(vr=rayon_raccord)
    call detrsd('CHAMP_GD', '&&RATU3D.CAXE_TUY')
    call detrsd('CHAMP_GD', '&&RATU3D.CAORIFI')
    call detrsd('CARTE', '&&RATU3D.NUME_MODE')
    call detrsd('CHAMP_GD', '&&RATU3D.PRAYO')
    call detrsd('RESUELEM', '&&RATU3D.COSF_UM')
    call detrsd('RESUELEM', '&&RATU3D.COSF_VM')
    call detrsd('RESUELEM', '&&RATU3D.COSF_WM')
    call detrsd('RESUELEM', '&&RATU3D.SINF_UM')
    call detrsd('RESUELEM', '&&RATU3D.SINF_VM')
    call detrsd('RESUELEM', '&&RATU3D.SINF_WM')
    call jedetr('&&RATU3D           .RELR')
    call jedetr('&&RATU3D           .RERR')
    call detrsd('CHAMP_GD', 'CH_DEPL_1')
    call detrsd('CHAMP_GD', 'CH_DEPL_2')
    call detrsd('CHAMP_GD', 'CH_DEPL_3')
    call detrsd('CHAMP_GD', 'CH_DEPL_4')
    call detrsd('CHAMP_GD', 'CH_DEPL_5')
    call detrsd('CHAMP_GD', 'CH_DEPL_6')
!
    call jedema()
end subroutine
