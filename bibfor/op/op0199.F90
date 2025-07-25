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
subroutine op0199()
    implicit none
! OPERATEUR CALCULANT LA FORCE AJOUTEE : CALC_FORC_AJOU
!
!---------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/cal152.h"
#include "asterfort/calmdg.h"
#include "asterfort/cresol.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mat152.h"
#include "asterfort/phi199.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rigflu.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ibid, nbmo, nbmode(1), ndble, indice, ifm, niv
    integer(kind=8) :: tabad(5), iadesc, iarefe, i, iadirg, imade
    integer(kind=8) :: iphi1, iphi2, iprsto, iret, itxsto
    integer(kind=8) :: itysto, itzsto, ivalk, ivale
    integer(kind=8) :: n1, nbLoad, n3, n4, n5, n6, n7, n8, n9
    real(kind=8) :: rbid, mij, cij, kij
    complex(kind=8) :: cbid
    aster_logical :: vrai
    character(len=2) :: model
    character(len=3) :: nd
    character(len=8) :: nomres, k8bid, modeMeca, phibar, moint, loadName
    character(len=8) :: modelFluid, matrAsse, materField, numgen, modgen
    character(len=14) :: numeDof, num, nugene
    character(len=16) :: typres, nomcom
    character(len=19) :: max, may, maz, chamno, solveu
    character(len=24) :: blanc, nocham, mateco

!
!-----------------------------------------------------------------
!
    call jemarq()
!
    call getres(nomres, typres, nomcom)
!
    call infmaj()
    call infniv(ifm, niv)
!
    nbmo = 0
    ndble = 0
    vrai = .true.
    solveu = '&&OP0199.SOLVEUR'
    nugene = ' '
    materField = ' '
    mateco = ' '
!
! --- RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getvid(' ', 'MODELE_FLUIDE', scal=modelFluid, nbret=n1)
    call getvid(' ', 'CHARGE', scal=loadName, nbret=nbLoad)
    call getvid(' ', 'MODELE_INTERFACE', scal=moint, nbret=n3)
    call getvid(' ', 'CHAM_MATER', scal=materField, nbret=n4)
    call getvid(' ', 'MODE_MECA', scal=modeMeca, nbret=n5)
    call getvid(' ', 'NUME_DDL_GENE', scal=numgen, nbret=n6)
    call getvid(' ', 'MODELE_GENE', scal=modgen, nbret=n7)
    call getvid(' ', 'POTENTIEL', scal=phibar, nbret=n8)
    call getvtx(' ', 'NOEUD_DOUBLE', scal=nd, nbret=n9)
!
! --- LECTURE DES PARAMETRES  SOLVEUR
!
    call cresol(solveu)
!
    if (n4 .ne. 0) call rcmfmc(materField, mateco, l_ther_=ASTER_FALSE)
!
    if (n6 .ne. 0) nugene = numgen
!
    if (n5 .ne. 0) then
        call rsorac(modeMeca, 'LONUTI', 0, rbid, k8bid, &
                    cbid, rbid, 'ABSOLU', nbmode, 1, &
                    ibid)
        nbmo = nbmode(1)
        call rsexch(' ', modeMeca, 'DEPL', 1, nocham, &
                    iret)
    end if
!
    if (n7 .ne. 0) then
        if (nd .eq. 'OUI') ndble = 1
    end if
!
    model = '  '
!
!--------------------------------------------------------------
! --- CALCUL DE LA MATRICE ASSEMBLEE DE RIGIDITE DU FLUIDE
!--------------------------------------------------------------
!
    call rigflu(modelFluid, mateco, &
                nbLoad, loadName, &
                solveu, numeDof, matrAsse)
!
!--------------------------------------------------------------
! CALCUL DES MATR_ELEM AX ET AY DANS L'OPTION FLUX_FLUI_X ET _Y
!---------------SUR LE MODELE INTERFACE(THERMIQUE)-------------
! CALCUL DES MATRICES MODALES BI POUR L OPTION AMOR_AJOU
!--------------------------------------------------------------
!
    call mat152('MASS_AJOU', model, moint, ivalk, &
                nbmo, max, may, maz, num)
!
    call jeexin('&&MAT152.MADE', iret)
    if (iret .gt. 0) call jeveuo('&&MAT152.MADE', 'E', imade)
!
!================================================================
! CALCUL ET STOCKAGE DES POTENTIELS INSTATIONNAIRES PHI1 ET PHI2
! CORRESPONDANT RESPECTIVEMENT AUX EFFETS INERTIELS
! ET AUX EFFETS D'AMORTISSEMENT ET DE RAIDEUR DU FLUIDE
! SUR LA STRUCTURE
!================================================================
!
    call phi199(model, materField, mateco, matrAsse, numeDof, &
                num, nbmo, solveu, indice, tabad)
!
!--------------------------------------------------------------
! VERIFICATION D EXISTENCE DE VECTEUR DE CHAMPS AUX NOEUDS CREES
! DS PHI152 ILS SERONT ENSUITE EXPLOITES DS CAL152 ENTRE AUTRES
! VECTEUR DE NOMS DU POTENTIEL INSTATIONNAIRE PHI1 : MASSE AJOU
! ON Y STOCKE LES NOMS DES POTENTIELS INSTATIONNAIRES POUR
! CHAQUE MODE DE STRUCTURE
!
    call jeexin('&&OP0199.PHI1', iret)
    if (iret .gt. 0) call jeveuo('&&OP0199.PHI1', 'E', iphi1)
    call jeexin('&&OP0199.PHI2', iret)
    if (iret .gt. 0) call jeveuo('&&OP0199.PHI2', 'E', iphi2)
!
!=====================================================================
!---------------------------------------------------------------------
!              CALCUL SUR MODELE GENERALISE
!---------------------------------------------------------------------
!=====================================================================
!
    if (n7 .gt. 0) then
        call calmdg(model, modgen, nugene, num, numeDof, &
                    matrAsse, materField, mateco, moint, ndble, &
                    itxsto, itysto, itzsto, iprsto, nbmo, &
                    iadirg)
    end if
!
!=============================================================
!--------REMPLISSAGE DU  .VALE : CALCUL DU VECTEUR AJOUTE
!=============================================================
!
!---------------------------------------------------------------
    if ((n7 .gt. 0) .or. (indice .eq. 1)) then
!
! CALCUL DU VECTEUR AJOUTE - PRODUITS SCALAIRES SUR MODELE
! GENERALISE - CAS DE LA SOUS-STRUCTURATION DYNAMIQUE
! OU BIEN CAS DE MODES RESTITUES SUR MAILLAGE SQUELETTE
!
        if (indice .eq. 1) then
            itxsto = tabad(1)
            itysto = tabad(2)
            itzsto = tabad(3)
            iprsto = tabad(4)
            iadirg = tabad(5)
            nbmo = nbmode(1)
        end if
    else
!
! --- CREATION DE L OBJET VECT_GENE RESULTAT
!
        call wkvect(nomres//'           .VALE', 'G V R', nbmo, ivale)
        call wkvect(nomres//'           .REFE', 'G V K24', 2, iarefe)
        call wkvect(nomres//'           .DESC', 'G V I', 3, iadesc)
        call jeecra(nomres//'           .DESC', 'DOCU', cval='VGEN')
!
! --- REMPLISSAGE DU .REFE ET .VALE
!
        zk24(iarefe) = modeMeca
        zk24(iarefe+1) = nugene
        zi(iadesc) = 1
        zi(iadesc+1) = nbmo
!
        do i = 1, nbmo
!
            blanc = ' '
            call cal152('MASS_AJOU', max, may, maz, model, &
                        blanc, iphi1, iphi2, imade, modeMeca, &
                        chamno, num, vrai, i, 1, &
                        mij, cij, kij)
!
            zr(ivale+i-1) = mij
!
        end do
    end if
!
!
    call jedetc('G', '&&RIGFLU', 1)
!
    call jedema()
end subroutine
