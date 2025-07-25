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
subroutine avgrma(vwork, tdisp, vnbpg, nbpgt, nbordr, &
                  nmaini, nbmap, numpaq, tspaq, nomcri, &
                  nomfor, grdvie, forvie, fordef, proaxe, &
                  cesr)
! person_in_charge: van-xuan.tran at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/loisem.h"
#include "asterc/lor8em.h"
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/avplcr.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedisp.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rnomat.h"
#include "asterfort/utmess.h"
#include "asterfort/vecnuv.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: tdisp, nbmap, vnbpg(nbmap), nbpgt, nbordr, nmaini
    integer(kind=8) :: numpaq, tspaq
    real(kind=8) :: vwork(tdisp)
    character(len=16) :: nomcri, proaxe, nomfor, forvie, grdvie
    character(len=19) :: cesr
    aster_logical :: fordef, post
! ---------------------------------------------------------------------
! BUT: DETERMINER LE PLAN DANS LEQUEL LE DOMMAGE EST MAXIMAL
! ---------------------------------------------------------------------
! ARGUMENTS:
! VWORK     IN    R  : VECTEUR DE TRAVAIL CONTENANT
!                      L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                      ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                      DU <<PAQUET>> DE MAILLES.
! TDISP     IN    I  : DIMENSION DU VECTEUR VWORK
! VNBPG     IN    I  : VECTEUR CONTENANT LE NOMBRE DE
!                      POINT DE GAUSS DE CHAQUE MAILLE DU MAILLAGE.
! NBPGT     IN    I  : NOMBRE TOTAL DE POINTS DE GAUSS A TRAITER.
! NBORDR    IN    I  : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                      STRUCTURE DE DONNEES RESULTAT.
! NMAINI    IN    I  : NUMERO DE LA 1ERE MAILLE DU <<PAQUET>> DE
!                      MAILLES COURANT.
! NBMAP     IN    I  : NOMBRE DE MAILLES DANS LE <<PAQUET>> DE
!                      MAILLES COURANT.
! NUMPAQ    IN    I  : NUMERO DU PAQUET DE MAILLES COURANT.
! TSPAQ     IN    I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                      COURANT.
! NOMCRI    IN    K16: NOM DU CRITERE AVEC PLANS CRITIQUES.
! PROAXE    IN    K16: TYPE DE PROJECTION (UN OU DEUX AXES).
! CESR      IN    K19: NOM DU CHAMP SIMPLE DESTINE A RECEVOIR LES
!                      RESULTATS.
!
! REMARQUE :
!  - LA TAILLE DU SOUS-PAQUET EST EGALE A LA TAILLE DU <<PAQUET>> DE
!    MAILLES DIVISEE PAR LE NOMBRE DE NUMERO D'ORDRE (NBORDR).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: nbvecm
    integer(kind=8) :: jcerd, jcerl, iret, icesd, icesl, icesv, ibid
    integer(kind=8) :: tneces, tdisp2(1), jvecpg, n, k
    integer(kind=8) :: ideb, dim, j, ngam, tab2(18), ifin
    integer(kind=8) :: l, nbpg, nbpgp, kwork, sompgw, imap
    integer(kind=8) :: ipg
    integer(kind=8) :: icmp, jad
    integer(kind=8) :: vali(2)
!
    real(kind=8) :: fatsoc, dgam, gamma, pi, dphi, tab1(18), phi0
    real(kind=8) :: vala, valb, coefpa, cudomx
    real(kind=8) :: nxm(2), nym(2), nzm(2)
    real(kind=8) :: vresu(24)
!
    character(len=8) :: chmat1, nommat
    character(len=10) :: optio
    character(len=19) :: chmat, cesmat
    real(kind=8), pointer :: vect_norma(:) => null()
    real(kind=8), pointer :: vect_tangu(:) => null()
    real(kind=8), pointer :: vect_tangv(:) => null()
    real(kind=8), pointer :: cerv(:) => null()
!
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
    data tab1/180.0d0, 60.0d0, 30.0d0, 20.0d0, 15.0d0, 12.857d0,&
     &             11.25d0, 10.588d0, 10.0d0, 10.0d0, 10.0d0, 10.588d0,&
     &             11.25d0, 12.857d0, 15.0d0, 20.0d0, 30.0d0, 60.0d0/
!
    data tab2/1, 3, 6, 9, 12, 14, 16, 17, 18, 18, 18, 17, 16, 14,&
     &           12, 9, 6, 3/
!
    pi = r8pi()
!-----------------------------------------------------------------------
!
    call jemarq()
!
! CONSTRUCTION DU VECTEUR NORMAL SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR U DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
! CONSTRUCTION DU VECTEUR V DANS LE PLAN TANGENT, SUR UNE DEMI SPHERE
!
    AS_ALLOCATE(vr=vect_norma, size=627)
    AS_ALLOCATE(vr=vect_tangu, size=627)
    AS_ALLOCATE(vr=vect_tangv, size=627)
!
!
! OBTENTION DES ADRESSES '.CESD', '.CESL' ET '.CESV' DU CHAMP SIMPLE
! DESTINE A RECEVOIR LES RESULTATS : DOMMAGE_MAX, COORDONNEES VECTEUR
! NORMAL CORRESPONDANT
!
    call jeveuo(cesr//'.CESD', 'L', jcerd)
    call jeveuo(cesr//'.CESL', 'E', jcerl)
    call jeveuo(cesr//'.CESV', 'E', vr=cerv)
!
! RECUPERATION MAILLE PAR MAILLE DU MATERIAU DONNE PAR L'UTILISATEUR
!
    call getvid(' ', 'CHAM_MATER', scal=chmat1, nbret=iret)
    chmat = chmat1//'.CHAMP_MAT'
    cesmat = '&&AVGRMA.CESMAT'
    call carces(chmat, 'ELEM', ' ', 'V', cesmat, &
                'A', iret)
    call jeveuo(cesmat//'.CESD', 'L', icesd)
    call jeveuo(cesmat//'.CESL', 'L', icesl)
    call jeveuo(cesmat//'.CESV', 'L', icesv)
!
! DEFINITION DU VECTEUR CONTENANT LES VALEURS DU CISAILLEMENT POUR TOUS
! LES INSTANTS ET TOUS LES PLANS
!
    tneces = 209*nbordr*2
    call jedisp(1, tdisp2)
    tdisp2(1) = (tdisp2(1)*loisem())/lor8em()
    if (tdisp2(1) .lt. tneces) then
        vali(1) = tdisp2(1)
        vali(2) = tneces
        call utmess('F', 'PREPOST5_8', ni=2, vali=vali)
    else
        call wkvect('&&AVGRMA.VECTPG', 'V V R', tneces, jvecpg)
    end if
!
! COEFFICIENT PERMETTANT D'UTILISER LES MEMES ROUTINES POUR LES
! CONTRAINTES ET LES DEFORMATIONS
!
    if ((nomcri(1:16) .eq. 'FATESOCI_MODI_AV') .or. fordef) then
        fatsoc = 1.0d4
    else
        fatsoc = 1.0d0
    end if
!
! CONSTRUCTION DES VECTEURS N, U ET V
!
    dgam = 10.0d0
!
    n = 0
    k = 1
    ideb = 1
    dim = 627
    do j = 1, 18
        gamma = (j-1)*dgam*r8dgrd()
        dphi = tab1(j)*r8dgrd()
        ngam = tab2(j)
        ifin = ngam
        phi0 = dphi/2.0d0
!
        call vecnuv(ideb, ifin, gamma, phi0, dphi, &
                    n, k, dim, vect_norma, vect_tangu, &
                    vect_tangv)
!
    end do
!
! CONSTRUCTION DU VECTEUR : CISAILLEMENT = F(NUMERO D'ORDRE) EN CHAQUE
! POINT DE GAUSS DU PAQUET DE MAILLES.
    l = 1
    nbpg = 0
    nbpgp = 0
    kwork = 0
    sompgw = 0
!
    do imap = nmaini, nmaini+(nbmap-1)
        if (imap .gt. nmaini) then
            kwork = 1
            sompgw = sompgw+vnbpg(imap-1)
        end if
        nbpg = vnbpg(imap)
! SI LA MAILLE COURANTE N'A PAS DE POINTS DE GAUSS, LE PROGRAMME
! PASSE DIRECTEMENT A LA MAILLE SUIVANTE.
        if (nbpg .eq. 0) then
            goto 400
        end if
!
        nbpgp = nbpgp+nbpg
        if ((l*int(nbpgt/10.0d0)) .lt. nbpgp) then
            write (6, *) numpaq, '   ', (nbpgp-nbpg)
            l = l+1
        end if
!
! RECUPERATION DU NOM DU MATERIAU AFFECTE A LA MAILLE COURANTE
! ET DES PARAMETRES ASSOCIES AU CRITERE CHOISI POUR LA MAILLE COURANTE.
!
        optio = 'DOMA_ELGA'
        call rnomat(icesd, icesl, icesv, imap, nomcri, &
                    ibid, ibid, ibid, optio, vala, &
                    valb, coefpa, nommat)
!
!
!
        nbvecm = 209
!
        post = .false.
!
        do ipg = 1, nbpg
!
!
! REMPLACER PAR AVPLCR
            call avplcr(nbvecm, vect_norma, vect_tangu, vect_tangv, nbordr, &
                        kwork, sompgw, vwork, tdisp, tspaq, &
                        ipg, nomcri, nomfor, grdvie, forvie, &
                        fordef, fatsoc, proaxe, nommat, vala, &
                        coefpa, post, cudomx, nxm, nym, &
                        nzm)
!
! RECUPERER LES RESULTATS
            do icmp = 1, 24
                vresu(icmp) = 0.0d0
            end do
!
            vresu(2) = nxm(1)
            vresu(3) = nym(1)
            vresu(4) = nzm(1)
            vresu(11) = cudomx
            vresu(13) = nxm(2)
            vresu(14) = nym(2)
            vresu(15) = nzm(2)
!
! 12. AFFECTATION DES RESULTATS DANS UN CHAM_ELEM SIMPLE
!
            do icmp = 1, 24
                call cesexi('C', jcerd, jcerl, imap, ipg, &
                            1, icmp, jad)
!
                ASSERT(jad .ne. 0)
                jad = abs(jad)
                zl(jcerl-1+jad) = .true.
                cerv(jad) = vresu(icmp)
!
            end do
!
        end do
400     continue
    end do
!
! MENAGE
!
    call detrsd('CHAM_ELEM_S', cesmat)
!
    AS_DEALLOCATE(vr=vect_norma)
    AS_DEALLOCATE(vr=vect_tangu)
    AS_DEALLOCATE(vr=vect_tangv)
!
    call jedema()
end subroutine
