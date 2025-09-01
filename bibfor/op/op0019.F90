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

subroutine op0019()
!
!
! --------------------------------------------------------------------------------------------------
!
!                O P E R A T E U R    AFFE_CARA_ELEM
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterf_types.h"
#include "asterfort/ace_crea_carte.h"
#include "asterfort/ace_masse_repartie.h"
#include "asterfort/aceaba.h"
#include "asterfort/aceaca.h"
#include "asterfort/aceaco.h"
#include "asterfort/aceadi.h"
#include "asterfort/aceagb.h"
#include "asterfort/aceama.h"
#include "asterfort/aceamb.h"
#include "asterfort/aceamr.h"
#include "asterfort/aceaor.h"
#include "asterfort/aceapc.h"
#include "asterfort/aceapf.h"
#include "asterfort/aceapo.h"
#include "asterfort/acearm.h"
#include "asterfort/acearp.h"
#include "asterfort/acecel.h"
#include "asterfort/aceinc.h"
#include "asterfort/acevba.h"
#include "asterfort/acevco.h"
#include "asterfort/acevdi.h"
#include "asterfort/acevgb.h"
#include "asterfort/acevma.h"
#include "asterfort/acevmb.h"
#include "asterfort/acevmr.h"
#include "asterfort/acevor.h"
#include "asterfort/acevpf.h"
#include "asterfort/acevpo.h"
#include "asterfort/acevrm.h"
#include "asterfort/acevrp.h"
#include "asterfort/alcart.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/checkCaraElem.h"
#include "asterfort/coqucf.h"
#include "asterfort/detrsd_vide.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/pmfd00.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
#include "asterfort/verif_affe.h"
#include "asterfort/verima.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)             :: elem_supp_num(ACE_NB_TYPE_ELEM)
    character(len=16)   :: elem_supp_nom(ACE_NB_TYPE_ELEM)
    integer(kind=8)             :: elem_supp_typ(ACE_NB_TYPE_ELEM)
    integer(kind=8)             :: nb_type_elem(ACE_NB_ELEMENT)
! --------------------------------------------------------------------------------------------------
!   Pour les cartes :
    type(cara_elem_carte)   :: info_carte(ACE_NB_CARTE)
!   Infomation sur le concept : maillage, modele, nb noeuds, nb mailles, ...
    type(cara_elem_info) :: info_concept
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbocc(ACE_NB_MCLEF)
    integer(kind=8) :: nbtout(ACE_NB_MCLEF)
    character(len=8) :: mclef_type
!
    integer(kind=8) :: ivr(4), nbcart, iret, jadr, ii, nbtel, ireponse, nbelemdi
    integer(kind=8) :: nbver, nlm, nlg, lxc, lxo, nln, nlj
    integer(kind=8) :: lxb, lxm, lxpf, lxgb, lxmb, lmax, ifm, niv, lxp, nbvm
    integer(kind=8) :: lxd, nboccd, lxrp, noemaf, lxrm, noemf2, nbmail, nbnoeu
    integer(kind=8) :: lxmr, noemf3
    integer(kind=8) :: npoutr, ncable, nbarre, nbdisc
    integer(kind=8) :: iclf, ioc, icle, ng, nocc, nocctout, nocctot
    integer(kind=8) :: depart
    aster_logical :: locaco, locagb, locamb, l_pmesh
    character(len=8) :: ver(3), nomu, nomo, noma
    character(len=16) :: concep, cmd, mclef, k16bid
    character(len=19) :: cartcf
    character(len=24) :: mlgnma, modnom, tmpncf, mlgnno
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), pointer            :: affe_mail(:) => null()
    character(len=24), pointer  :: grp_lmax(:) => null()
! --------------------------------------------------------------------------------------------------
    call jemarq()
!   CALL ONERRF('ABORT', K16BID, IRET)
    iret = 0
! --------------------------------------------------------------------------------------------------
!   Récupération des arguments de la commande
    call getres(nomu, concep, cmd)
! --------------------------------------------------------------------------------------------------
!   Modèle
    call getvid(' ', 'MODELE', scal=nomo, nbret=nbvm)
!   Enregistre le nom du modèle dans la SD de AFFE_CARA_ELEM
    call wkvect(nomu//'.MODELE', 'G V K8', 1, jadr)
    zk8(jadr) = nomo
!   Récupération du nom du maillage associé
    call dismoi('NOM_MAILLA', nomo, 'MODELE', repk=noma)
    l_pmesh = isParallelMesh(noma)
!   Construction des noms jeveux du concept maillage associé
    mlgnma = noma//'.TYPMAIL'
    mlgnno = noma//'.COORDO    .VALE'
!   Nombre de mailles du maillage
    call jelira(mlgnma, 'LONMAX', nbmail)
!   Nombre de noeuds du maillage
    call jelira(mlgnno, 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
!   Récupération de la dimension géométrique du modèle
    call dismoi('DIM_GEOM', nomo, 'MODELE', repi=ireponse)
    info_concept%dimmod = ireponse
    if (ireponse .ge. 100) then
        ireponse = ireponse-100
        info_concept%dimmod = 1
    end if
    if (ireponse .ge. 20) then
        ireponse = ireponse-20
        info_concept%dimmod = 2
    end if
    if (ireponse .eq. 3) info_concept%dimmod = 3
!   Mémorisation des informations pour ne plus le refaire
    info_concept%nomu = nomu
    info_concept%concept = concep
    info_concept%commande = cmd
    info_concept%modele = nomo
    info_concept%maillage = noma
    info_concept%nbmail = nbmail
    info_concept%nbnoeu = nbnoeu
!
! --------------------------------------------------------------------------------------------------
!   Pour faire des vérifications de cohérence d'affectation
    call getvtx(' ', 'VERIF', nbval=2, vect=ver, nbret=nbver)
    ivr(:) = 0
!   ivr(1)=1    : vérification MAILLES
!   ivr(2)      : libre
!   ivr(3)=niv  : niveau d'impression
!   ivr(4)=ifm  : unité d'impression

    if (nbver .gt. 0) then
        do ii = 1, nbver
            if (ver(ii) .eq. 'MAILLE  ') ivr(1) = 1
        end do
    end if
!   Récupération du niveau d'impression
    call infmaj()
    call infniv(ifm, niv)
    ivr(3) = niv
    ivr(4) = ifm
    info_concept%ivr(:) = ivr(:)
!
! --------------------------------------------------------------------------------------------------
!   Occurence des mots clefs facteur
    do ii = 1, ACE_NB_MCLEF
        call getfac(ACE_MCLEF(ii), nbocc(ii))
    end do
!
! --------------------------------------------------------------------------------------------------
!   Initialisation des éléments pouvant être affetés :
!       elem_supp_nom  elem_supp_typ
    nbtel = 0
    do ii = nbtel+1, nbtel+ACE_NB_POUTRE
        elem_supp_nom(ii) = ACE_EL_POUTRE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_POUTRE
    end do
    nbtel = nbtel+ACE_NB_POUTRE
!
    do ii = nbtel+1, nbtel+ACE_NB_DISCRET
        elem_supp_nom(ii) = ACE_EL_DISCRET(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_DISCRET
    end do
    nbtel = nbtel+ACE_NB_DISCRET
!
    do ii = nbtel+1, nbtel+ACE_NB_COQUE
        elem_supp_nom(ii) = ACE_EL_COQUE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_COQUE
    end do
    nbtel = nbtel+ACE_NB_COQUE
!
    do ii = nbtel+1, nbtel+ACE_NB_CABLE
        elem_supp_nom(ii) = ACE_EL_CABLE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_CABLE
    end do
    nbtel = nbtel+ACE_NB_CABLE
!
    do ii = nbtel+1, nbtel+ACE_NB_BARRE
        elem_supp_nom(ii) = ACE_EL_BARRE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_BARRE
    end do
    nbtel = nbtel+ACE_NB_BARRE
!
    do ii = nbtel+1, nbtel+ACE_NB_MASSIF
        elem_supp_nom(ii) = ACE_EL_MASSIF(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_MASSIF
    end do
    nbtel = nbtel+ACE_NB_MASSIF
!
    do ii = nbtel+1, nbtel+ACE_NB_GRILLE
        elem_supp_nom(ii) = ACE_EL_GRILLE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_GRILLE
    end do
    nbtel = nbtel+ACE_NB_GRILLE
!
    do ii = nbtel+1, nbtel+ACE_NB_MEMBRANE
        elem_supp_nom(ii) = ACE_EL_MEMBRANE(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_MEMBRANE
    end do
    nbtel = nbtel+ACE_NB_MEMBRANE
!
    do ii = nbtel+1, nbtel+ACE_NB_THHMM
        elem_supp_nom(ii) = ACE_EL_THHMM(ii-nbtel)
        elem_supp_typ(ii) = ACE_NU_THHMM
    end do
    nbtel = nbtel+ACE_NB_THHMM
    ASSERT(nbtel .eq. ACE_NB_TYPE_ELEM)
! --------------------------------------------------------------------------------------------------
!   Récupération des numéros des types éléments
    do ii = 1, ACE_NB_TYPE_ELEM
        call jenonu(jexnom('&CATA.TE.NOMTE', elem_supp_nom(ii)), elem_supp_num(ii))
    end do
!
! --------------------------------------------------------------------------------------------------
!   Vérification de l'existence des GROUP_MA et MAILLE déclarés
!       Après cette vérification il n'est plus nécessaire d'utiliser les routines
!           verima   getvem(getvtx+verima)
!   On ne fait pas la vérification si on a TOUT (= 'OUI' seule valeur possible)
!   sous le mot-clé facteur.
!
!   Comptage des GROUP_MA et MAILLE. Pour ne pas faire des ALLOCATE dans la boucle.
    lmax = 10
    nocc = 0
    nocctout = 0
    nocctot = 0
    do iclf = 1, ACE_NB_MCLEF
        do ioc = 1, nbocc(iclf)
!   Le superviseur doit avoir déjà vérifié qu'on ne peut avoir qu'une seule occurrence
!   du mot-clé facteur si on trouve TOUT='OUI' sous celui-ci (fonction compat_syntax).
!   On met quand même des assert au cas où.
            call getvtx(ACE_MCLEF(iclf), 'TOUT', iocc=ioc, nbval=1, nbret=nbtout(iclf))
            ASSERT((nbtout(iclf) .eq. 1) .or. (nbtout(iclf) .eq. 0))
            if (nbtout(iclf) .eq. 1) then
                ASSERT(nbocc(iclf) .eq. 1)
                nocctout = nocctout+1
            else
                do icle = 1, ACE_NB_GRMA_MA
                    ii = MCLEF_GRP_MA(icle+(iclf-1)*ACE_NB_GRMA_MA)
                    if (ii .ne. ACE_NOTHING) then
                        mclef = ACE_GRMA_MA(ii)
                        call getvtx(ACE_MCLEF(iclf), mclef, iocc=ioc, nbval=0, nbret=ng)
                        lmax = max(lmax, -ng)
                        nocc = max(nocc, -ng)
                    end if
                end do
            end if
            nocctot = nocctot+nbocc(iclf)
        end do
    end do
!   Si on a moins d'occurrence de TOUT = 'OUI' que de mot-clés facteur, et qu'on a pas
!   GROUP_MA défini, alors il manque des affections sous au moins
!   un des mot-clés facteur donné par l'utilisateur.
    if ((nocc .le. 0) .and. (nocctout .lt. nocctot)) then
        call utmess('F', 'AFFECARAELEM_2')
    end if
!   Vérification
    AS_ALLOCATE(vk24=grp_lmax, size=lmax)
    do iclf = 1, ACE_NB_MCLEF
        do ioc = 1, nbocc(iclf)
            if (nbtout(iclf) .eq. 0) then
                do icle = 1, ACE_NB_GRMA_MA
                    ii = MCLEF_GRP_MA(icle+(iclf-1)*ACE_NB_GRMA_MA)
                    if (ii .ne. ACE_NOTHING) then
                        mclef = ACE_GRMA_MA(ii)
                        mclef_type = ACE_GRMA_TY(ii)
                        call getvtx(ACE_MCLEF(iclf), mclef, iocc=ioc, nbval=lmax, vect=grp_lmax, &
                                    nbret=ng)
                        call verima(noma, grp_lmax, ng, mclef_type)
                    end if
                end do
            end if
        end do
    end do
!
! --------------------------------------------------------------------------------------------------
!   Pour mémoriser les mailles affectées.
!           Si affe_mail(i)= -elem_supp_num     la maille est affectée et a été traitée
!                          = elem_supp_num      la maille est affectée et pas traitée
!                          = 0                  la maille n'est pas affectée
    AS_ALLOCATE(vi=affe_mail, size=nbmail)
! --------------------------------------------------------------------------------------------------
!   Création des cartes utilisées
    call ace_crea_carte(info_concept, info_carte)
!
! --------------------------------------------------------------------------------------------------
!   Vérification de la syntaxe pour :
!       ACE_CABLE   : éléments CABLE CABLE_POULIE
!                           call acevca(nbocc(ACE_CABLE), nlm, nlg, iret)
!
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS POUTRE
    lxp = 0
    if (nbocc(ACE_POUTRE) .ne. 0) then
        call acevpo(nbocc(ACE_POUTRE), nlm, nlg, iret)
        lxp = max(nlm, nlg)
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS COQUE
    lxc = 0
    if (nbocc(ACE_COQUE) .ne. 0) then
        call acevco(nbocc(ACE_COQUE), nlg, iret)
        lxc = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ORIENTATIONS DES ELEMENTS
    lxo = 0
    if (nbocc(ACE_ORIENTATION) .ne. 0) then
        call acevor(nbocc(ACE_ORIENTATION), nlg, iret)
        lxo = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS BARRE
    lxb = 0
    if (nbocc(ACE_BARRE) .ne. 0) then
        call acevba(nbocc(ACE_BARRE), nlg, iret)
        lxb = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS MASSIF :
    lxm = 0
    if (nbocc(ACE_MASSIF) .ne. 0) then
        call acevma(nbocc(ACE_MASSIF), nlg)
        lxm = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS POUTRE_FLUI
    lxpf = 0
    if (nbocc(ACE_POUTRE_FLUI) .ne. 0) then
        if (nbocc(ACE_POUTRE) .eq. 0) then
            call utmess('F', 'MODELISA5_56')
        end if
        call acevpf(nbocc(ACE_POUTRE_FLUI), nlg)
        lxpf = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS "GRILLE"
    lxgb = 0
    if (nbocc(ACE_GRILLE) .ne. 0) then
        call acevgb(nbocc(ACE_GRILLE), nlm, nlg)
        lxgb = max(nlm, nlg)
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS "MEMBRANE"
    lxmb = 0
    if (nbocc(ACE_MEMBRANE) .ne. 0) then
        call acevmb(nbocc(ACE_MEMBRANE), nlg)
        lxmb = nlg
    end if
! --------------------------------------------------------------------------------------------------
!   LONGUEUR MAXIMUM D UNE LISTE DE MAILLE/NOEUD/GROUP_MA/GROUP_NO
    lmax = max(lmax, lxp, lxc, lxo, lxb, lxm, lxpf, lxgb, lxmb)
!
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA SYNTAXE DES ELEMENTS DISCRET
    lxd = 0
    if (nbocc(ACE_DISCRET) .ne. 0 .or. nbocc(ACE_DISCRET_2D) .ne. 0) then
        nboccd = nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)
        if (nbocc(ACE_DISCRET) .ne. 0) k16bid = ACE_MCLEF(ACE_DISCRET)
        if (nbocc(ACE_DISCRET_2D) .ne. 0) k16bid = ACE_MCLEF(ACE_DISCRET_2D)
        call acevdi(nboccd, noma, nomo, k16bid, nlm, nlg, nln, nlj, iret)
        lxd = max(nlm, nln, nlg, nlj)
        lmax = max(lmax, lxd)
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA DIMENSION DES RAIDEURS REPARTIES
    lxrp = 0
    if (nbocc(ACE_RIGI_PARASOL) .ne. 0) then
        call acevrp(nbocc(ACE_RIGI_PARASOL), noma, lxrp, noemaf)
        lmax = max(lmax, lxrp)
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA DIMENSION DES RAIDEURS MISS
    lxrm = 0
    if (nbocc(ACE_RIGI_MISS_3D) .ne. 0) then
        call acevrm(nbocc(ACE_RIGI_MISS_3D), noma, lxrm, noemf2)
        lmax = max(lmax, lxrm)
    end if
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA DIMENSION DES MASSES REPARTIES
    lxmr = 0
    if (nbocc(ACE_MASS_AJOU) .ne. 0) then
        call acevmr(nbocc(ACE_MASS_AJOU), noma, lxmr, noemf3)
        lmax = max(lmax, lxmr)
    end if
!
! --------------------------------------------------------------------------------------------------
!   COMPTEUR D'ELEMENTS ET VERIFICATION COHERENCE DES AFFECTATIONS
    call acecel(noma, nomo, nbocc, elem_supp_num, elem_supp_typ, &
                nb_type_elem, affe_mail, iret)
!
    npoutr = nb_type_elem(ACE_NU_POUTRE)
    ncable = nb_type_elem(ACE_NU_CABLE)
    nbarre = nb_type_elem(ACE_NU_BARRE)
    nbdisc = nb_type_elem(ACE_NU_DISCRET)
!
    if (iret .ne. 0) then
        call utmess('F', 'MODELISA5_57')
    end if
!
! --------------------------------------------------------------------------------------------------
!   VERIFICATION DE LA BONNE  AFFECTATION  DES  CARACTERISTIQUES
!     POUR TOUTES LES MAILLES ET NOEUDS AFFECTES , IMPR SI DEMANDE
!     INCREMENTATION DES COMPTEURS D APPELS A NOCART (DISCRET,COQUE,
!     DEFI_ARC,CABLE,POUTRE,BARRE)
    iret = 0
    call aceinc(noma, nomo, elem_supp_num, nbocc, ivr, &
                locaco, locagb, locamb, affe_mail, lmax, iret)
    if (iret .ne. 0) then
        call utmess('F', 'MODELISA5_59')
    end if
!     FABRICATION DE LA CARTE COMMUNE A TOUS LES ELEMENTS LINEIQUE
!     S'IL Y EN A D'AFFECTE
    nbcart = 0
    if (nbocc(ACE_POUTRE)+nbocc(ACE_BARRE)+nbocc(ACE_CABLE) .ne. 0) then
        nbcart = npoutr+nbarre+ncable
        if (nbcart .gt. 0) then
            cartcf = nomu//'.CVENTCXF'
            call alcart('G', cartcf, noma, 'VENTCX_F')
        end if
    end if
    if ((nbocc(ACE_POUTRE) .eq. 0) .and. (npoutr .ne. 0)) then
        call utmess('A', 'MODELISA5_60')
    end if
    if ((nbocc(ACE_BARRE) .eq. 0) .and. (nbarre .ne. 0)) then
        call utmess('A', 'MODELISA5_61')
    end if
    if ((nbocc(ACE_CABLE) .eq. 0) .and. (ncable .ne. 0)) then
        call utmess('A', 'MODELISA5_62')
    end if
!
! --------------------------------------------------------------------------------------------------
!   Traitement des masses réparties
    call ace_masse_repartie(nbocc(ACE_MASS_REP), info_concept, grp_lmax, lmax, info_carte, &
                            nbdisc, affe_mail)
!
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES ORIENTATIONS AUX ELEMENTS POUTRES ET DISCRETS  ET
!     BARRES ET AFFECTATION DE LA CARTE ORIENTATION
    if (nbocc(ACE_POUTRE) .ne. 0 .or. nbocc(ACE_DISCRET) .ne. 0 .or. nbocc(ACE_DISCRET_2D) .ne. 0 &
        .or. nbocc(ACE_BARRE) .ne. 0 .or. nbocc(ACE_RIGI_PARASOL) .ne. 0) then
        call aceaor(noma, nomo, lmax, ACE_NB_POUTRE, &
                    elem_supp_num, elem_supp_nom, ivr, nbocc)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES AUX ELEMENTS POUTRES
    if (nbocc(ACE_POUTRE) .ne. 0) then
!        NBEPO + NBEDI + NBECO + NBECA + NBEBA + NBEMA + NBEGB
        depart = 1
        call aceapo(noma, nomo, lmax, npoutr, nbocc(ACE_POUTRE), &
                    ACE_MCLEF(ACE_POUTRE), ACE_NB_POUTRE, &
                    elem_supp_num(depart), ivr, affe_mail)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES EPAISSEURS/COURBURES/ANGLES AUX ELEMENTS COQUES
    if (nbocc(ACE_COQUE) .ne. 0) then
        call aceaco(nomu, noma, lmax, locagb, locamb, nbocc(ACE_COQUE))
    end if
! --------------------------------------------------------------------------------------------------
!   Affectation des coefficients de correction pour les COUDES
    if (nbocc(ACE_POUTRE) .ne. 0) then
        call aceapc(nomu, noma, lmax, nbocc(ACE_POUTRE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES SECTIONS AUX ELEMENTS CABLE :
    if (nbocc(ACE_CABLE) .ne. 0) then
        call aceaca(nomu, noma, lmax, nbocc(ACE_CABLE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES AUX ELEMENTS BARRE
    if (nbocc(ACE_BARRE) .ne. 0) then
!        NBEPO + NBEDI + NBECO + NBECA + NBEBA + NBEMA + NBEGB
        depart = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ACE_NB_CABLE+1
        call aceaba(noma, nomo, lmax, nbarre, nbocc(ACE_BARRE), &
                    ACE_MCLEF(ACE_BARRE), ACE_NB_BARRE, &
                    elem_supp_num(depart), ivr, affe_mail)
    end if

! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES REPERES AUX ELEMENTS THERMIQUES ET MECANIQUES
    if (nbocc(ACE_MASSIF) .ne. 0) then
        call aceama(nomu, noma, lmax, nbocc(ACE_MASSIF))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES REPERES AUX ELEMENTS POUTRE_FLUI
    if (nbocc(ACE_POUTRE_FLUI) .ne. 0) then
        call aceapf(nomu, noma, lmax, nbocc(ACE_POUTRE_FLUI))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX RAIDEURS REPARTIES
    if (nbocc(ACE_RIGI_PARASOL) .ne. 0) then
        call acearp(info_concept, lmax, noemaf, nbocc(ACE_RIGI_PARASOL), info_carte, ivr, &
                    affe_mail)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX ELEMENTS DISCRETS
    if (nbocc(ACE_DISCRET) .ne. 0 .or. nbocc(ACE_DISCRET_2D) .ne. 0) then
        nboccd = nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)
        if (nbocc(ACE_DISCRET) .ne. 0) k16bid = ACE_MCLEF(ACE_DISCRET)
        if (nbocc(ACE_DISCRET_2D) .ne. 0) k16bid = ACE_MCLEF(ACE_DISCRET_2D)
        call aceadi(noma, nomo, k16bid, lmax, nboccd, info_carte, ivr)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT "GRILLE"
    if (nbocc(ACE_GRILLE) .ne. 0) then
        call aceagb(nomu, noma, locamb, nbocc(ACE_GRILLE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX RAIDEURS MISS
    if (nbocc(ACE_RIGI_MISS_3D) .ne. 0) then
        call acearm(info_concept, lmax, nbocc(ACE_RIGI_MISS_3D), info_carte, ivr)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT "MEMBRANE"
    if (nbocc(ACE_MEMBRANE) .ne. 0) then
        call aceamb(nomu, noma, lmax, nbocc(ACE_MEMBRANE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX MASSES REPARTIES
    if (nbocc(ACE_MASS_AJOU) .ne. 0) then
        call aceamr(info_concept, lmax, nbocc(ACE_MASS_AJOU), info_carte, ivr)
    end if
! --------------------------------------------------------------------------------------------------
!   COMPACTAGE DE LA CARTE : '.CVENTCXF'
    if (nbcart .gt. 0) then
!        PAS APPELE POUR UNE SURCHARGE "FINE" MAIS POUR LE COMPACTAGE
        call tecart(cartcf)
!        DESTRUCTION DES CHAMPS
        tmpncf = cartcf//'.NCMP'
        call jedetr(tmpncf)
        tmpncf = cartcf//'.VALV'
        call jedetr(tmpncf)
    end if
!
!     POUR LES COQUES, GRILLES IL PEUT EXISTER UNE CARTE FONCTION
!     IL FAUT L'EVALUER ET METTRE LE RESULTAT DANS LA CARTE DES REELS
    if ((nbocc(ACE_COQUE) .ne. 0) .or. (nbocc(ACE_GRILLE) .ne. 0)) then
        call coqucf(nomu)
    end if
!
! --------------------------------------------------------------------------------------------------
!   TRAITEMENT DES MOTS CLES
!           MULTIFIBRE  /  GEOM_FIBRE
!           COQUE       /  COQUE_NCOU
!           GRILLE      /  COQUE_NCOU
!           MEMBRANE    /  COQUE_NCOU
!           POUTRE      /  TUYAU_NCOU  TUYAU_NSEC
    call pmfd00()
! --------------------------------------------------------------------------------------------------

! - Some checks by element (calcul)
    call checkCaraElem(nomo, nomu)

!   Certaines cartes peuvent etre vides : il faut les detruire.
    do ii = 1, ACE_NB_CARTE
        call detrsd_vide('CARTE', info_carte(ii)%nom_carte)
    end do
!
!   Destruction sélective des CARTES, si elles n'ont pas lieu d'être
    nbelemdi = nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)+nbocc(ACE_RIGI_PARASOL)+ &
               nbocc(ACE_RIGI_MISS_3D)+nbocc(ACE_MASS_AJOU)+nbocc(ACE_MASS_REP)
    if (nbelemdi .eq. 0) then
        do ii = 1, ACE_NB_CARTE
            if (ACE_CARTE(3+(ii-1)*ACE_NB_CARTE_CMP) .eq. 'DISCRET') then
                call detrsd('CHAMP', info_carte(ii)%nom_carte)
            end if
        end do
    end if
! --------------------------------------------------------------------------------------------------
!   Audit assignments
    call verif_affe(modele=nomo, sd=nomu)
!
    AS_DEALLOCATE(vi=affe_mail)
    AS_DEALLOCATE(vk24=grp_lmax)
!
    call jedema()
end subroutine
