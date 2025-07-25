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
subroutine op0152()
    implicit none
! AUTEUR : G.ROUSSEAU
! OPERATEUR CALCULANT LA MASSE AJOUTEE, L'AMORTISSEMENT
!  ET LA RIGIDITE AJOUTEE EN THEORIE POTENTIELLE : CALC_MATR_AJOU
!  SUR BASE MODALE DE LA STRUCTURE DANS LE VIDE
!---------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/cal152.h"
#include "asterfort/calmdg.h"
#include "asterfort/chpver.h"
#include "asterfort/cresol.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mag152.h"
#include "asterfort/mamodg.h"
#include "asterfort/mat152.h"
#include "asterfort/phi152.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rigflu.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/ualfva.h"
#include "asterfort/ver152.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: vrai
    integer(kind=8) :: ldblo, hc, ibid
    integer(kind=8) :: nbmo, nbmode, ndble, indice, tabad(5)
    integer(kind=8) :: i, j, jdesc
    integer(kind=8) :: iadirg, iblo, ierd
    integer(kind=8) :: imade
    integer(kind=8) :: iphi1, iphi2, iprsto, iret, itxsto
    integer(kind=8) :: itysto, itzsto, ivalk, kterm
    integer(kind=8) :: jsmdi, jsmde, jsmhc, n1bloc, n2bloc
    integer(kind=8) :: nbid, nbloc, nterm
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, n7, n9, n10, n12, n13, n14
    integer(kind=8) :: ifm, niv, tmod(1), nbLoad
    real(kind=8) :: mij, cij, kij
    real(kind=8) :: bid, ebid
    character(len=2) :: model
    character(len=3) :: nd
    character(len=8) :: nomres, k8bid, modeMeca, phibar
    character(len=8) :: modelFluid, moint, matrAsse, materField
    character(len=8) :: loadName, numgen, modgen
    character(len=9) :: option
    character(len=14) :: numeDof, num, nugene
    character(len=16) :: typres, nomcom
    character(len=19) :: max, may, maz, chamno
    character(len=19) :: stomor, solveu, nomr19
    character(len=24) :: nomcha, nocham
    character(len=24) :: mateco, phib24
    complex(kind=8) :: cbid
!
!-----------------------------------------------------------------
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    vrai = .true.
    solveu = '&&OP0152.SOLVEUR'
!
    call getres(nomres, typres, nomcom)
    nomr19 = nomres
!
!----------RECUPERATION DES ARGUMENTS DE LA COMMANDE--------------
!
!
    materField = ' '
    call getvid(' ', 'MODELE_FLUIDE', scal=modelFluid, nbret=n1)
    call getvid(' ', 'CHARGE', scal=loadName, nbret=nbLoad)
    call getvid(' ', 'MODELE_INTERFACE', scal=moint, nbret=n3)
    call getvid(' ', 'CHAM_MATER', scal=materField, nbret=n4)
    modeMeca = ' '
    call getvid(' ', 'MODE_MECA', scal=modeMeca, nbret=n5)
    call getvid(' ', 'CHAM_NO', nbval=0, nbret=n6)
    call getvid(' ', 'NUME_DDL_GENE', scal=numgen, nbret=n9)
    call getvid(' ', 'MODELE_GENE', scal=modgen, nbret=n10)
    call getvid(' ', 'POTENTIEL', scal=phibar, nbret=n12)
    call getvtx(' ', 'OPTION', scal=option, nbret=n13)
    call getvtx(' ', 'NOEUD_DOUBLE', scal=nd, nbret=n14)
!
! LECTURE DES PARAMETRES DONNES APRES LE MOT CLE FACTEUR SOLVEUR
!
    call cresol(solveu)
!
! VERIFICATIONS SUPPLEMENTAIRES
!
    call ver152(option, modelFluid, moint, n12, model)
!
! EXTRACTION DU POTENTIEL PERMANENT DES VITESSES
!
    if (option .eq. 'AMOR_AJOU' .or. option .eq. 'RIGI_AJOU') then
        call rsexch(' ', phibar, 'TEMP', 0, phib24, &
                    iret)
    end if
!
!
!     CAS NUME_DDL_GENE PRESENT
!     ------------------------------------
    if (n9 .ne. 0) then
        nugene = numgen
        stomor = nugene//'.SMOS'
!
    else
!       CREATION D UN NUME_DDL_GENE
!       ------------------------------------
        nugene = nomres
        stomor = nugene//'.SMOS'
!
        nbmode = -n6
        nbloc = 1
!
        call wkvect(stomor//'.SMDI', 'G V I', nbmode, jsmdi)
        nterm = 0
        do i = 1, nbmode
            nterm = nterm+i
            zi(jsmdi+i-1) = nterm
        end do
!
        call wkvect(stomor//'.SMHC', 'G V S', nterm, jsmhc)
        kterm = 0
        do i = 1, nbmode
            do j = 1, i
                kterm = kterm+1
                zi4(jsmhc-1+kterm) = int(j, 4)
            end do
        end do
!
        call wkvect(stomor//'.SMDE', 'G V I', 6, jsmde)
        zi(jsmde-1+1) = nbmode
        zi(jsmde-1+2) = nterm
        zi(jsmde-1+3) = nbloc
!
    end if
!
!
    if (n6 .ne. 0) then
        n7 = -n6
        vrai = .false.
    else
        n7 = 0
    end if
!
!--------- RECUPERATION DU MATERIAU FLUIDE----------------------------
    if (n4 .ne. 0) then
        call rcmfmc(materField, mateco, l_ther_=ASTER_FALSE)
    else
        mateco = ' '
    end if
!
!--------CALCUL DE LA MATRICE ASSEMBLEE DE RIGIDITE DU FLUIDE---------
!
    call rigflu(modelFluid, mateco, &
                nbLoad, loadName, &
                solveu, numeDof, matrAsse)
!
!=====================================================================
!---------------- ALTERNATIVE CHAMNO OU MODE_MECA OU---------
!-----------------------------MODELE-GENE--------------------
!=====================================================================
!
!----------------------------------------------------------------
    if (n5 .gt. 0) then
        call rsorac(modeMeca, 'LONUTI', ibid, bid, k8bid, &
                    cbid, ebid, 'ABSOLU', tmod, 1, &
                    nbid)
        nbmode = tmod(1)
        nbmo = nbmode
        call rsexch(' ', modeMeca, 'DEPL', 1, nomcha, &
                    iret)
        nocham = nomcha
    else
        if (n7 .gt. 0) then
            nbmo = n7
! 1ERE CREATION DE VECTEUR DE NOMS DES CHAMPS DE DEPL_R
! REPRESENTANT LES MODES
! EN CAS D'UTILISATION DU MOT CLE CHAM-NO, CECI POUR
! MAT152
            call jecreo('&&OP0152.VEC', 'V V K8')
            call jeecra('&&OP0152.VEC', 'LONMAX', n7)
            call jeveuo('&&OP0152.VEC', 'E', ivalk)
            call getvid(' ', 'CHAM_NO', nbval=n7, vect=zk8(ivalk), nbret=n6)
            nocham = zk8(ivalk)
            call chpver('F', nocham, 'NOEU', 'DEPL_R', ierd)
        end if
    end if
!--------------------------------------------------------------
! CALCUL DES MATR_ELEM AX ET AY DANS L'OPTION FLUX_FLUI_X ET _Y
!---------------SUR LE MODELE INTERFACE(THERMIQUE)-------------
! CALCUL DES MATRICES MODALES BI POUR L OPTION AMOR_AJOU
!--------------------------------------------------------------
    call mat152(option, model, moint, ivalk, &
                nbmo, max, may, maz, num)
    call jeexin('&&MAT152.MADE', iret)
    if (iret .gt. 0) call jeveuo('&&MAT152.MADE', 'E', imade)
! DESTRUCTION DU VECTEUR DE NOMS DES DEPL-R POUR RECREATION DS
! PHI152
    call jeexin('&&OP0152.VEC', iret)
    if (iret .gt. 0) call jedetr('&&OP0152.VEC')
!
!================================================================
! CALCUL ET STOCKAGE DES POTENTIELS INSTATIONNAIRES PHI1 ET PHI2
! CORRESPONDANT RESPECTIVEMENT AUX EFFETS INERTIELS
! ET AUX EFFETS D'AMORTISSEMENT ET DE RAIDEUR DU FLUIDE
! SUR LA STRUCTURE
!================================================================
    call phi152(model, option, materField, mateco, phib24, &
                matrAsse, numeDof, num, nbmode, solveu, &
                indice, tabad)
!
! VERIFICATION D EXISTENCE DE VECTEUR DE CHAMPS AUX NOEUDS CREES
! DS PHI152 ILS SERONT ENSUITE EXPLOITES DS CAL152 ENTRE AUTRES
! VECTEUR DE NOMS DU POTENTIEL INSTATIONNAIRE PHI1 : MASSE AJOU
! ON Y STOCKE LES NOMS DES POTENTIELS INSTATIONNAIRES POUR
! CHAQUE MODE DE STRUCTURE
    call jeexin('&&OP0152.PHI1', iret)
    if (iret .gt. 0) call jeveuo('&&OP0152.PHI1', 'E', iphi1)
!
! VECTEUR DE NOMS DU POTENTIEL INSTATIONNAIRE PHI2 : AMOR AJOU
! RIGI AJOU
    call jeexin('&&OP0152.PHI2', iret)
    if (iret .gt. 0) call jeveuo('&&OP0152.PHI2', 'E', iphi2)
!
! VECTEUR DE NOMS DES CHAMPS DE DEPL_R REPRESENTANT LES MODES
! EN CAS D'UTILISATION DU MOT CLE CHAM-NO
    call jeexin('&&OP0152.VEC', iret)
    if (iret .gt. 0) call jeveuo('&&OP0152.VEC', 'E', ivalk)
!
!
!
!================================================================
!----------- CREATION DE LA MATR_ASSE_GENE    -------------------
!----------- CONTENANT LA MASSE AJOUTEE RESULTAT   --------------
!================================================================
!
    call mag152(n9, n10, nomres, nugene, modeMeca, &
                modgen, nbloc, indice)
!
!=====================================================================
!---------------------------------------------------------------------
!              CALCUL SUR MODELE GENERALISE
!---------------------------------------------------------------------
!=====================================================================
!
    if (n10 .gt. 0) then
        if (nd .eq. 'OUI') then
            ndble = 1
        else
            ndble = 0
        end if
        call calmdg(model, modgen, nugene, num, numeDof, &
                    matrAsse, materField, mateco, moint, ndble, &
                    itxsto, itysto, itzsto, iprsto, nbmo, &
                    iadirg)
!
    end if
!
!
!=============================================================
!--------REMPLISSAGE DU  .VALE : CALCUL DE LA MASSE AJOUTEE
!=============================================================
!
!---------------------IMPRESSION DES RESULTATS------------------
!
    if (niv .gt. 1) then
        if (option .eq. 'MASS_AJOU') then
            write (ifm, *) '        '
            write (ifm, *) '          =======MATRICE DE MASSE AJOUTEE======='
            if (n10 .gt. 0) then
                write (ifm, *) '           ========HORS DDL DE LAGRANGE==='
            end if
        end if
        if (option .eq. 'AMOR_AJOU') then
            write (ifm, *) '        '
            write (ifm, *) '         =====MATRICE D AMORTISSEMENT AJOUTE====='
        end if
        if (option .eq. 'RIGI_AJOU') then
            write (ifm, *) '        '
            write (ifm, *) '        =======MATRICE DE RIGIDITE AJOUTEE======='
        end if
    end if
!---------------------------------------------------------------
    if ((n10 .gt. 0) .or. (indice .eq. 1)) then
!
! CALCUL DES MASSES AJOUTEES - PRODUITS SCALAIRES SUR MODELE
! GENERALISE - CAS DE LA SOUS-STRUCTURATION DYNAMIQUE
! OU BIEN CAS DE MODES RESTITUES SUR MAILLAGE SQUELETTE
!
        if (indice .eq. 1) then
            itxsto = tabad(1)
            itysto = tabad(2)
            itzsto = tabad(3)
            iprsto = tabad(4)
            iadirg = tabad(5)
            nbmo = nbmode
        end if
!
        call mamodg(model, stomor, nomres, itxsto, itysto, &
                    itzsto, iprsto, iadirg, nbmo, max, &
                    may, maz, nbloc)
    else
!
! CAS CLASSIQUE
!
        call jeveuo(stomor//'.SMDI', 'L', jsmdi)
        call jeveuo(stomor//'.SMDE', 'L', jsmde)
        call jeveuo(stomor//'.SMHC', 'L', jsmhc)
!
!     BOUCLE SUR LES BLOCS DE LA MATRICE ASSEMBLEE GENE
!
        do iblo = 1, nbloc
            call jecroc(jexnum(nomr19//'.UALF', iblo))
            call jeveuo(jexnum(nomr19//'.UALF', iblo), 'E', ldblo)
!----------------------------------------------------------------
!
!         BOUCLE SUR LES COLONNES DE LA MATRICE ASSEMBLEE
!
            n1bloc = 1
            n2bloc = zi(jsmde)
!
            do i = n1bloc, n2bloc
!
                hc = zi(jsmdi-1+i)
                if (i .gt. 1) hc = hc-zi(jsmdi-1+i-1)
!
                do j = (i-hc+1), i
!
!----------------------------------------------------------------
! ICI ON CALCULE LA MASSE AJOUTEE SUR UN MODELE GENERALISE
!----------------------------------------------------------------
!
!------------------------------------------------------------------
! ICI ON CALCULE LA MATRICE DE MASSE AJOUTEE SOIT SUR UN MODE_MECA
! SOIT SUR UN CHAM_NO
!------------------------------------------------------------------
!
                    if (n7 .gt. 0) then
                        chamno = zk8(ivalk+i-1)
                    end if
!
                    call cal152(option, max, may, maz, model, &
                                phib24, iphi1, iphi2, imade, modeMeca, &
                                chamno, num, vrai, i, j, &
                                mij, cij, kij)
!
!
!
!-----------STOCKAGE DANS LA MATR_ASSE_GENE  ------
!
!        CAS DE LA PROJECTION MODALE OU CHAM_NO
!
!
                    if (option .eq. 'MASS_AJOU') then
                        zr(ldblo+zi(jsmdi+i-1)+j-i-1) = mij
                    end if
                    if (option .eq. 'AMOR_AJOU') then
                        zr(ldblo+zi(jsmdi+i-1)+j-i-1) = cij
                    end if
                    if (option .eq. 'RIGI_AJOU') then
                        zr(ldblo+zi(jsmdi+i-1)+j-i-1) = kij
                    end if
!
!===============================================================
!---------------IMPRESSION DES RESULTATS------------------------
!===============================================================
!
                    if (niv .gt. 1) then
                        if (((n9 .gt. 0) .and. (n5 .ne. 0)) .or. (n6 .ne. 0)) then
                            if (option .eq. 'MASS_AJOU') then
                                write (ifm, 350) i, j, zr(ldblo+j+(i-1)*i/ &
                                                          2-1)
                            end if
                            if (option .eq. 'AMOR_AJOU') then
                                write (ifm, 450) i, j, zr(ldblo+j+(i-1)*i/ &
                                                          2-1)
                            end if
                            if (option .eq. 'RIGI_AJOU') then
                                write (ifm, 550) i, j, zr(ldblo+j+(i-1)*i/ &
                                                          2-1)
                            end if
350                         format(18x, 'M', 2 i 4, 1x, '=', 1x, d 12.5)
450                         format(18x, 'C', 2 i 4, 1x, '=', 1x, d 12.5)
550                         format(18x, 'K', 2 i 4, 1x, '=', 1x, d 12.5)
                        end if
                    end if
                end do
            end do
        end do
    end if
!
    if (niv .gt. 1) then
!
        write (ifm, *) '              ============================'
        write (ifm, *) '              =======FIN IMPRESSION======='
        write (ifm, *) '              ============================'
!
    end if
!
!
!   -- on repasse au stockage morse qui est le stockage normal :
!   ------------------------------------------------------------
    call ualfva(nomres, 'G')
!
!   -- on corrige l'objet .DESC :
!   ------------------------------------------------------------------------------
    call jelira(nomr19//'.CONL', 'LONMAX', n1)
    call jelira(jexnum(nomr19//'.VALM', 1), 'LONMAX', n2)
    call jeveuo(nomr19//'.DESC', 'E', jdesc)
    zi(jdesc) = 2
    zi(jdesc+1) = n1
    if (n2 .eq. n1) then
        zi(jdesc+2) = 1
    else if (n2 .eq. n1*(n1+1)/2) then
        zi(jdesc+2) = 2
    else
        zi(jdesc+2) = 3
    end if
!
!
    call jedetc('G', '&&RIGFLU', 1)
    call jedetc('G', '&&CALMAA', 1)
!
    call jedema()
end subroutine
