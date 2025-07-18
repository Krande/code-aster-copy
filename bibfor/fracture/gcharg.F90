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
subroutine gcharg(modele, lischa, chvolu, ch1d2d, ch2d3d, &
                  chpres, chepsi, chpesa, chrota, lfonc, &
                  time, iord)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/calcG_type.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpver.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcchar.h"
#include "asterfort/gcfonc.h"
#include "asterfort/gcsele.h"
#include "asterfort/getvtx.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lisccc.h"
#include "asterfort/lisdef.h"
#include "asterfort/lislch.h"
#include "asterfort/lislcm.h"
#include "asterfort/lislnf.h"
#include "asterfort/lisltc.h"
#include "asterfort/lisltf.h"
#include "asterfort/lisnnb.h"
#include "asterfort/lisnnl.h"
#include "asterfort/mefor0.h"
#include "asterfort/mepres.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: iord
    character(len=8) :: modele
    character(len=19) :: lischa
    character(len=19) :: chvolu, ch1d2d, ch2d3d, chpres
    character(len=19) :: chepsi, chpesa, chrota
    aster_logical :: lfonc
    real(kind=8) :: time
!
! OPERATEUR CALC_G
!
! CREATION DES CHAMPS PROVENANT DES CHARGEMENTS
!
! ----------------------------------------------------------------------
!
! IN  LISCHA : LISTE DES CHARGES
! IN  MODELE : NOM DU MODELE
! IN  TIME   : INSTANT
! IN  IORD   : NUMERO D'ORDRE CORRESPONDANT A TIME
! OUT CHVOLU : CARTE MODIFIEE POUR CHARGEMENT FORCE_INTERNE
! OUT CH1D2D : CARTE MODIFIEE POUR CHARGEMENT FORCE_CONTOUR
! OUT CH2D3D : CARTE MODIFIEE POUR CHARGEMENT FORCE_FACE
! OUT CHPRES : CARTE MODIFIEE POUR CHARGEMENT PRES_REP
! OUT CHEPSI : CARTE MODIFIEE POUR CHARGEMENT EPSI_INIT
! OUT CHPESA : CARTE MODIFIEE POUR CHARGEMENT PESANTEUR
! OUT CHROTA : CARTE MODIFIEE POUR CHARGEMENT ROTATION
! OUT LFONC  : .TRUE. AU MOINS UNE CHARGE EST DE TYPE FONCTION
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: znbenc
    parameter(znbenc=60)
    integer(kind=8) :: tabaut(znbenc)
!
    character(len=24) :: k24bid
    character(len=24) :: oldfon, cepsi, epselno, modelLigrel
    integer(kind=8) :: jfonci
    integer(kind=8) :: ichar, nbchar, ig, iret, inga, occur, ier, i
    character(len=8) :: charge, typech, nomfct, newfct, ng
    character(len=6) :: nomobj
    character(len=16) :: typfct, motcle, nomcmd, phenom
    character(len=13) :: prefob
    integer(kind=8) :: motclc(2)
    aster_logical :: lfchar, lfmult, lformu, lccomb, lpchar
    integer(kind=8) :: nbauth, nbnaut, mclaut(2), iposit
    integer(kind=8) :: iprec, ibid, itypob(2), ibid2(2)
    character(len=19) :: carteo, cartei
    aster_logical :: lvolu, l1d2d, l2d3d, lpres
    aster_logical :: lepsi, lpesa, lrota
    aster_logical :: lfvolu, lf1d2d, lf2d3d, lfpres
    aster_logical :: lfepsi, lfpesa, lfrota
    integer(kind=8), pointer :: desc(:) => null()
    character(len=8), pointer :: p_vale_epsi(:) => null()
    character(len=8) :: list_option(NB_MAX_OPT)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! - INITIALISATIONS
!
    occur = 0
    lfonc = .false.
    nomcmd = 'CALC_G'
    phenom = 'MECANIQUE'
    cepsi = '&&GCHARG.CEPSI'
    epselno = '&&GCHARG.EPSELNO'
    lvolu = .false.
    l1d2d = .false.
    l2d3d = .false.
    lpres = .false.
    lepsi = .false.
    lpesa = .false.
    lrota = .false.
    lfvolu = .false.
    lf1d2d = .false.
    lf2d3d = .false.
    lfpres = .false.
    lfepsi = .false.
    lfpesa = .false.
    lfrota = .false.
    lformu = .false.
    lpchar = .false.
    lccomb = .false.
    list_option(:) = ' '
    call lisnnb(lischa, nbchar)
!
!   Recuperation du LIGREL
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=modelLigrel)
!
! - STOCKAGE DES TYPES DE CHARGE (FONCTION OU PAS)
!
    oldfon = '&&GCHARG.FONCI'
    if (nbchar .gt. 0) call wkvect(oldfon, 'V V L', nbchar, jfonci)
!
!
! - UNE DES CHARGES EST-ELLE UNE FONCTION ? -> REMPLISSAGE DE LFONC
!
    do ichar = 1, nbchar
! ----- NOM DE LA CHARGE
        call lislch(lischa, ichar, charge)
! ----- TYPE DE LA CHARGE (REELLE OU FONCTION)
        call lisltc(lischa, ichar, typech)
        if (typech(1:4) .eq. 'FONC') then
            lfonc = .true.
            goto 20
        end if
20      continue
    end do
!
! ----- BOUCLE SUR LES CHARGES
!
    do ichar = 1, nbchar
!
! ----- NOM DE LA CHARGE
!
        call lislch(lischa, ichar, charge)
!
! ----- ENTIERS CODES POUR LES MOTS-CLEFS DE LA CHARGE
!
        call lislcm(lischa, ichar, motclc)
!
! ----- TYPE DE LA CHARGE (REELLE OU FONCTION)
!
        call lisltc(lischa, ichar, typech)
        if (typech(1:4) .eq. 'FONC') then
            lfchar = .true.
        else
            lfchar = .false.
        end if
!
! ----- CE CHARGEMENT EST FONCTION
!
        zl(jfonci+ichar-1) = lfchar
!
! ----- FONCTION MULTIPLICATRICE
!
        call lisltf(lischa, ichar, typfct)
        call lislnf(lischa, ichar, nomfct)
        lfmult = .false.
        if (typfct(1:5) .eq. 'FONCT') then
            lfmult = .true.
        end if
!
! ----- CHARGEMENTS UTILISES OU NON DANS CALC_G
!
        call lisccc(nomcmd, motclc, nbauth, nbnaut, mclaut)
!
        if (nbauth .ne. 0) then
!
! --------- DECODAGE
!
            call isdeco(mclaut, tabaut, znbenc)
!
! --------- BOUCLE SUR LES MOTS-CLEFS
!
            do iposit = 1, znbenc
                if (tabaut(iposit) .eq. 1) then
!
! ----------------- MOT-CLEF DE LA CHARGE
!
                    call lisdef('MOTC', k24bid, iposit, motcle, ibid2)
                    if (motcle .eq. 'DIRI_DUAL') goto 12
!
! ----------------- PREFIXE DE L'OBJET DE LA CHARGE
!
                    call lisnnl(phenom, charge, prefob)
!
! ----------------- CARTE D'ENTREE
!
                    call lisdef('CART', motcle, ibid, nomobj, itypob)
                    ASSERT(itypob(1) .eq. 1)
                    cartei = prefob(1:13)//nomobj(1:6)
!
! ----------------- Si le champ PRE_EPSI est de type CHAM_NO ou CARTE, il faut récupérer
!                   le vrai nom du champ correspondant (voir load_neum_spec)
!
                    if (nomobj .eq. '.EPSIN') then
                        call jeveuo(prefob(1:13)//'.EPSIN.DESC', 'L', vi=desc)
                        ig = desc(1)
                        call jenuno(jexnum('&CATA.GD.NOMGD', ig), ng)
! ----------------- recuperation du nom du champ stocké dans la carte "bidon"
                        if (ng .eq. 'NEUT_K8') then
                            call jeveuo(prefob(1:13)//'.EPSIN.VALE', 'L', vk8=p_vale_epsi)
                            cartei = p_vale_epsi(1)
                        end if
! ----------------- transformation si champ ELGA -> ELNO
                        call chpver('C', cartei(1:19), 'ELGA', 'EPSI_R', inga)
!
                        if (inga == 0) then
                            occur = occur+1
                            ASSERT(occur <= 1)
!               traitement du champ pour les elements finis classiques
                            call detrsd('CHAMP', cepsi)
                            call alchml(modelLigrel, 'CALC_G', 'PEPSINR', 'V', cepsi, &
                                        iret, ' ')
                            call chpchd(cartei(1:19), 'ELNO', cepsi, 'OUI', 'V', &
                                        epselno, modele)
                            call chpver('F', epselno(1:19), 'ELNO', 'EPSI_R', iret)
                            cartei(1:19) = epselno(1:19)
                        end if
                    end if
!
! ----------------- SELECTION SUIVANT TYPE
!
                    call gcsele(motcle, chvolu, ch1d2d, ch2d3d, chpres, &
                                chepsi, chpesa, chrota, lvolu, l1d2d, &
                                l2d3d, lpres, lepsi, lpesa, lrota, &
                                lfchar, lfvolu, lf1d2d, lf2d3d, lfpres, &
                                lfepsi, lfpesa, lfrota, carteo, lpchar, &
                                lccomb)
!
!
!------------------ Interdiction d'un chargement PRE_EPSI avec l'option G de CALC_G
!                   (sans re-calcul de contraintes, le résultat est faux)
!
                    if (lepsi .or. lfepsi) then
                        call getvtx(' ', 'OPTION', nbret=ier)
                        if (ier == 1) then
                            call getvtx(' ', 'OPTION', scal=list_option(1))
                            if (list_option(1) .eq. 'G') then
                                call utmess('F', 'RUPTURE0_91')
                            end if
                        else
                            call getvtx(' ', 'OPTION', nbval=-ier, vect=list_option)
                            do i = 1, -ier
                                if (list_option(i) .eq. 'G') then
                                    call utmess('F', 'RUPTURE0_91')
                                end if
                            end do
                        end if
                    end if
!
!
!
! ----------------- PREPARATION NOM DE LA FONCTION RESULTANTE
!
                    call gcfonc(ichar, iord, cartei, lfchar, lfmult, &
                                newfct, lformu)
!
! ----------------- CAS PARTICULIER (CARACTERE GENERIQUE A D'AUTRE CAS A DETERMINER):
! ----------------- LE CHARGEMENT DE TYPE PRES_REP DOIT ETRE UNE FONCTION S'IL EXITE UN
! ----------------- AUTRE CAHRGEMENT FONCTION
! ----------------- DANS LE CAS CONTRAIRE, ON CREE UN CHAMPS FONCTION UNIFORME
                    if (motcle .eq. 'PRES_REP') then
                        if (lfonc .and. (.not. lfpres)) then
                            call jeveuo(cartei//'.VALE', 'L', iret)
                            call mepres(modele, cartei, lfonc, zr(iret), 0.d0)
                            lfpres = .true.
                            lfchar = .true.
                        end if
                    end if
!
! ----------------- CONSTRUIT LA CARTE A PARTIR DU CHARGEMENT
!
                    call gcchar(ichar, iprec, time, carteo, lfchar, &
                                lpchar, lformu, lfmult, lccomb, cartei, &
                                nomfct, newfct, oldfon)
!
!
12                  continue
                end if
            end do
        end if
    end do
!
! - SI ABSENCE D'UN CHAMP DE FORCES, CREATION D'UN CHAMP NUL
!
    if (.not. lvolu) then
        call mefor0(modele, chvolu, lfonc)
    end if
    if (.not. l1d2d) then
        call mefor0(modele, ch1d2d, lfonc)
    end if
    if (.not. l2d3d) then
        call mefor0(modele, ch2d3d, lfonc)
    end if
    if (.not. lpres) then
        call mepres(modele, chpres, lfonc, 0.d0, 0.d0)
    end if
!
!
    call jedetr(oldfon)
    call jedema()
end subroutine
