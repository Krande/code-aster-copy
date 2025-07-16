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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine dichoc_endo_ldc(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
!        COMPORTEMENT GRILLE ASSEMBLAGE COMBUSTIBLE
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type, only: te0047_dscr
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/ldc_dichoc_endo.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvala.h"
#include "asterfort/rk5adp.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imatri, ivarim, irep, ifono, icontp, ivarip, iadzi, iazk24, iiter, iterat
    integer(kind=8) :: icarcr, idf, ipi, imate, jmater, nbmater
    integer(kind=8) :: ii, kk
    character(len=24) :: messak(6)
!
    integer(kind=8) :: imater, igeom, icontm, jdc, ivitp, idepen, iviten, jtm, jtp
    integer(kind=8) :: iretlc
    real(kind=8) :: klc(for_discret%neq*for_discret%neq), klv(for_discret%nbt)
    real(kind=8) :: dvl(for_discret%neq), dpe(for_discret%neq), dve(for_discret%neq)
    real(kind=8) :: fl(for_discret%neq)
    real(kind=8) :: raide(6), force(1)
    real(kind=8) :: r8bid
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: k8bid
    aster_logical :: rigi, resi, Prediction, Dynamique
! --------------------------------------------------------------------------------------------------
!   Pour le matériau
    integer(kind=8), parameter :: nbre1 = 3
    real(kind=8) :: valre1(nbre1)
    integer(kind=8) :: codre1(nbre1)
    character(len=16) :: materiau
! --------------------------------------------------------------------------------------------------
!   Pour l'intégration de la loi de comportement
    real(kind=8) :: temps0, temps1, dtemps
!   Paramètres de la loi :     jeu
    integer(kind=8), parameter :: ijeu = 1
    integer(kind=8), parameter :: nbpara = 1, nbpain = 2*3+3
    real(kind=8) :: ldcpar(nbpara)
    integer(kind=8) :: ldcpai(nbpain)
    character(len=8) :: ldccar(1)
!   Équations du système
    integer(kind=8), parameter :: nbequa = 6
    real(kind=8) :: y0(nbequa), dy0(nbequa), resu(nbequa*2), errmax, ynorme(nbequa)
    integer(kind=8) :: nbdecp
!   Variables internes
    integer(kind=8), parameter :: nbvari = 5, nbcorr = 4, idebut = nbvari
    integer(kind=8) :: Correspond(nbcorr)
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: xl(6), xd(3), rignor, deplac, evoljeu0, evoljeu1, xjeu
    real(kind=8) :: LongDist, Dist12, forceref, rigidref, deplaref
! --------------------------------------------------------------------------------------------------
!   Paramètres associés au matériau codé
    integer(kind=8), parameter :: lmat = 9, lfct = 10
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!   RIGI_MECA_TANG ->        DSIDEP        -->  RIGI
!   FULL_MECA      ->  SIGP  DSIDEP  VARP  -->  RIGI  RESI
!   RAPH_MECA      ->  SIGP          VARP  -->        RESI
    rigi = (for_discret%option(1:4) .eq. 'RIGI' .or. for_discret%option(1:4) .eq. 'FULL')
    resi = (for_discret%option(1:4) .eq. 'RAPH' .or. for_discret%option(1:4) .eq. 'FULL')
    iret = 0
! --------------------------------------------------------------------------------------------------
    call jevech('PCOMPOR', 'L', vk16=compor)
!   Seulement en 3D, sur un segment, avec seulement de la translation
    if ((for_discret%nomte(1:12) .ne. 'MECA_DIS_T_L') .or. (for_discret%ndim .ne. 3) .or. &
        (for_discret%nno .ne. 2) .or. (for_discret%nc .ne. 3)) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_22', nk=5, valk=messak)
    end if
! --------------------------------------------------------------------------------------------------
!   Paramètres en entrée
    call jevech('PMATERC', 'L', imater)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PCADISK', 'L', jdc)
!   on recupere le no de l'iteration de newton
    call jevech('PITERAT', 'L', iiter)
    iterat = zi(iiter)
!   Seulement en repère local : irep = 2
    call infdis('REPK', irep, r8bid, k8bid)
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_5', nk=5, valk=messak)
    end if
!   les caractéristiques sont toujours dans le repère local. on fait seulement une copie
    b_n = to_blas_int(for_discret%nbt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
!   Récupère les termes diagonaux de la matrice de raideur
    call diraidklv(for_discret%nomte, raide, klv)
! --------------------------------------------------------------------------------------------------
!   Champ de vitesse
    Dynamique = ASTER_FALSE
    call tecach('ONO', 'PVITPLU', 'L', iretlc, iad=ivitp)
    if (iretlc .eq. 0) then
        call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        Dynamique = ASTER_TRUE
    else
        dvl(:) = 0.0d0
    end if
!   Champ de déplacement d'entrainement
    call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
    if (iretlc .eq. 0) then
        call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
    else
        dpe(:) = 0.0d0
    end if
!   Champ de vitesse d'entrainement
    call tecach('ONO', 'PVITENT', 'L', iretlc, iad=iviten)
    if (iretlc .eq. 0) then
        call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
    else
        dve(:) = 0.d0
    end if
! --------------------------------------------------------------------------------------------------
!   Variables internes
    call jevech('PVARIMR', 'L', ivarim)
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
        varpl(ii) = varmo(ii)
    end do
! --------------------------------------------------------------------------------------------------
!   Coordonnees du discret dans le repère local
    xl(:) = 0.0
    call utpvgl(for_discret%nno, 3, for_discret%pgl, zr(igeom), xl)
! --------------------------------------------------------------------------------------------------
!   Adresse de la SD mater
    jmater = zi(imater)
!   Nombre de matériau sur la maille : 1 seul autorisé
    nbmater = zi(jmater)
    ASSERT(nbmater .eq. 1)
!   Adresse du matériau codé
    imate = jmater+zi(jmater+nbmater+1)
!   Recherche du matériau dans la SD compor
    materiau = 'DIS_CHOC_ENDO'
    ipi = 0
    do kk = 1, zi(imate+1)
        if (zk32(zi(imate)+kk-1) (1:16) .eq. materiau) then
            ipi = zi(imate+2+kk-1)
            goto 10
        end if
    end do
!   Le matériau n'est pas trouvé
    messak(1) = for_discret%nomte
    messak(2) = for_discret%option
    messak(3) = compor(INCRELAS)
    messak(4) = materiau
    call tecael(iadzi, iazk24)
    messak(5) = zk24(iazk24-1+3)
    call utmess('F', 'DISCRETS_7', nk=5, valk=messak)
10  continue
! --------------------------------------------------------------------------------------------------
!   Le bloc d'instruction précédent permet de déterminer : ipi
    ldcpai(:) = -1
    idf = zi(ipi)+zi(ipi+1)
!   Pour les fonctions :
!       Nombre de point     : zi(ii)
!       Adresse des valeurs : zi(ii+2)
    do kk = 1, zi(ipi+2)
        if ('FXP   ' .eq. zk16(zi(ipi+3)+idf+kk-1)) then
            ii = ipi+lmat-1+lfct*(kk-1)
            ldcpai(1) = zi(ii)
            ldcpai(2) = zi(ii+2)
        else if ('RIGIP ' .eq. zk16(zi(ipi+3)+idf+kk-1)) then
            ii = ipi+lmat-1+lfct*(kk-1)
            ldcpai(3) = zi(ii)
            ldcpai(4) = zi(ii+2)
        else if ('AMORP ' .eq. zk16(zi(ipi+3)+idf+kk-1)) then
            ii = ipi+lmat-1+lfct*(kk-1)
            ldcpai(5) = zi(ii)
            ldcpai(6) = zi(ii+2)
        end if
    end do
    ASSERT((ldcpai(1) .eq. ldcpai(3)) .and. (ldcpai(1) .eq. ldcpai(5)))
    if (ldcpai(1) .lt. 5) then
        messak(1) = materiau
        messak(2) = '[FXP|RIGIP|AMORP]'
        call utmess('F', 'DISCRETS_64', nk=2, valk=messak)
    end if
! --------------------------------------------------------------------------------------------------
!   Le premier point de FXP     : forceref
!                       RIGIP   : rigidref
!   Le déplacement de référence : deplaref = forceref/rigidref
    forceref = zr(ldcpai(1)+ldcpai(2))
    rigidref = zr(ldcpai(3)+ldcpai(4))
    deplaref = forceref/rigidref
! --------------------------------------------------------------------------------------------------
!   loi de comportement non-linéaire : récupération du temps + et - , calcul de dt
    call jevech('PINSTPR', 'L', jtp)
    call jevech('PINSTMR', 'L', jtm)
    temps0 = zr(jtm)
    temps1 = zr(jtp)
    dtemps = temps1-temps0
!   contrôle de rk5 : découpage successif, erreur maximale
    call jevech('PCARCRI', 'L', icarcr)
!   nombre d'itérations maxi (ITER_INTE_MAXI=-20 par défaut)
    nbdecp = abs(nint(zr(icarcr)))
!   tolérance de convergence (RESI_INTE=1.0E-06 par défaut)
    errmax = zr(icarcr+2)
!
! --------------------------------------------------------------------------------------------------
!   Paramètres de la loi de comportement
    valre1(:) = 0.0
    call rcvala(jmater, ' ', 'DIS_CHOC_ENDO', 0, ' ', &
                [0.0d0], 3, ['DIST_1   ', 'DIST_2   ', 'CRIT_AMOR'], valre1, codre1, &
                1)
!   Calcul du jeu final
!   Longueur du discret
    xd(1:3) = xl(1+for_discret%ndim:2*for_discret%ndim)-xl(1:for_discret%ndim)
    LongDist = xd(1)
    Dist12 = -valre1(1)-valre1(2)
!   Amortissement dans le critère ou pas : 1 ou 2
    ldcpai(9) = nint(valre1(3))
!
    ldcpar(ijeu) = LongDist+Dist12
!   Loi complète
    ldcpai(7) = 1
!
!   Traitement de l'évolution du jeu
!   Pour l'instant ce n'est pas utile
!       ldcpai(8) = 0 Pas de fonction d'évolution du jeu
!                   1 fonction d'évolution du jeu
!       ldcpai(7) = 1 Loi complète
!                 = 2 Loi sans frottement, ni amortissement
    ldcpai(8) = 0
    evoljeu0 = 1.0
    evoljeu1 = 1.0
! --------------------------------------------------------------------------------------------------
!
!   équations du système :
!              1   2   3   4     5    6
!       y0   : ux  fx  vx  pcum  uxan jeu
!       vari :         3   1     2    4
!
    Correspond(:) = [4, 5, 3, 6]
    y0(:) = 0.0; dy0(:) = 0.0
    do ii = 1, nbcorr
        y0(Correspond(ii)) = varmo(ii)
    end do
!   Les grandeurs et leurs dérivées
    y0(1) = (for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1)+dpe(1+for_discret%nc)-dpe(1))
    y0(2) = zr(icontm)
    dy0(1) = (for_discret%dul(1+for_discret%nc)-for_discret%dul(1))/dtemps
    dy0(3) = (dvl(1+for_discret%nc)-dvl(1)+dve(1+for_discret%nc)-dve(1)-y0(3))/dtemps
!   calcul de la vitesse d'évolution du jeu
    if (nint(varmo(idebut)) .eq. 0) then
        y0(6) = LongDist
    end if
    dy0(6) = Dist12*(evoljeu1-evoljeu0)/dtemps
!   RIGIP est obligatoire Kp( y0(4) ) = raideur
!       La fonction raideur et sa dérivée
!   Traitement de l'évolution du jeu
    valre1(:) = 0.0
    call rcvala(jmater, ' ', 'DIS_CHOC_ENDO', 1, 'PCUM', &
                [y0(4)], 1, ['RIGIP'], valre1, codre1, &
                0, nan='NON')
    rignor = valre1(1)
!
    force(:) = 0.0
!   Prédiction en dynamique, on retourne les efforts précédents
    Prediction = ((iterat .eq. 1) .and. resi .and. Dynamique)
!
!   Soit on intègre le jeu soit on prend sa valeur
!       ldcpai(8) = 1 : intégration du jeu
!       ldcpai(8) = 0 : valeur finale
    xjeu = y0(6)*ldcpai(8)+ldcpar(ijeu)*(1.0-ldcpai(8))
!
    if (Prediction .and. (ldcpai(8) .eq. 0)) then
        r8bid = y0(1)+dy0(1)*dtemps+xjeu-y0(5)
        raide(1) = 0.0
        if (r8bid <= 0.0) then
            raide(1) = rignor
            if (for_discret%option .eq. 'RAPH_MECA') then
                force(1) = raide(1)*r8bid
            end if
        end if
        goto 888
    end if
!
!   Norme pour le critère d'erreur
    ynorme(1) = deplaref
    ynorme(2) = forceref
    ynorme(3) = deplaref/dtemps
    ynorme(4) = deplaref
    ynorme(5) = deplaref
    ynorme(6) = deplaref
!
    call rk5adp(nbequa, ldcpar, ldcpai, ldccar, temps0, &
                dtemps, nbdecp, errmax, y0, dy0, &
                ldc_dichoc_endo, resu, iret, ynorme)
!   resu(1:nbeq)            : variables intégrées
!   resu(nbeq+1:2*nbeq)     : d(resu)/d(t) a t+dt
    if (iret .ne. 0) goto 999
!
!   Les variables internes
    do ii = 1, nbcorr
        varpl(ii) = resu(Correspond(ii))
    end do
    varpl(idebut) = 1.0
!
    call rcvala(jmater, ' ', 'DIS_CHOC_ENDO', 1, 'PCUM', &
                [resu(4)], 1, ['RIGIP'], valre1, codre1, &
                0, nan='NON')
    raide(1) = valre1(1)
!
    deplac = resu(1)-y0(1)
    force(1) = min(0.0, resu(2))
    if (abs(deplac) > r8prem()) then
        raide(1) = min(raide(1), abs((force(1)-y0(2))/deplac))
    end if
!
888 continue
! --------------------------------------------------------------------------------------------------
!   Actualisation de la matrice tangente : klv(i,i) = raide(i)
    call diklvraid(for_discret%nomte, klv, raide)
!   Actualisation de la matrice quasi-tangente
    if (rigi) then
        call jevech('PMATUUR', 'E', imatri)
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatri))
    end if
!
!   Efforts généralisés, Forces nodales
!       Sortie : efforts généralisés
!       Calcul : Forces nodales  (Mise à jour après)
    if (resi) then
        call jevech('PCONTPR', 'E', icontp)
!       demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, for_discret%neq)
!       calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', for_discret%neq, klc, for_discret%dul, fl)
!       efforts généralisés aux noeuds 1 et 2 (repère local)
!       on change le signe des efforts sur le premier noeud pour les MECA_DIS_TR_L et MECA_DIS_T_L
        do ii = 1, for_discret%nc
            zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
            zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc&
                                             &)
            fl(ii) = fl(ii)-zr(icontm-1+ii)
            fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
        end do
        zr(icontp-1+1) = force(1)
        zr(icontp-1+1+for_discret%nc) = force(1)
        fl(1) = -force(1)
        fl(1+for_discret%nc) = force(1)
! Sortie : forces nodales
! forces nodales aux noeuds 1 et 2 (repère global)
        call jevech('PVECTUR', 'E', ifono)
        call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
! Sortie : variables internes
! mise à jour des variables internes
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, nbvari
            zr(ivarip+ii-1) = varpl(ii)
            zr(ivarip+ii-1+nbvari) = varpl(ii)
        end do
    end if
!
999 continue
end subroutine
