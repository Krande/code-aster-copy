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
subroutine coefmo(typflu, zrigi, nbm, nmode, indic, &
                  x, pulsc, vgap, xsi0, veci1, &
                  vecr1, vecr2, vecr3, vecr4, vecr5, &
                  xmf, xkf, xcf)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/bijmoc.h"
#include "asterfort/bmocca.h"
#include "asterfort/cajgr2.h"
#include "asterfort/cfrott.h"
#include "asterfort/coefam.h"
#include "asterfort/coefra.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/kajgr2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=8) :: typflu
    real(kind=8) :: xsi0, vgap, x(2), pulsc, xmf, xcf, rkf
    real(kind=8) :: vecr1(*), vecr2(*), vecr3(*), vecr4(*), vecr5(*)
    integer(kind=8) :: nbm, nmode, indic, veci1(*)
    aster_logical :: zrigi
    complex(kind=8) :: xkf
!
! ----------------------------------------------------------------------
    real(kind=8) :: depi, vrmin, vrmax, vred
    real(kind=8) :: mcf0, kaj1, kaj2
    integer(kind=8) :: jvired, jcompt, jextr, jvrzo, jalarm, iret
    integer(kind=8) :: izone, nzone, n1, n2
    complex(kind=8) :: bii, biie, poids, z
    character(len=24) :: fsic, fsvi, fsvr, nom1, nom2, nom3, nom4, nom5
! ----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, imasse, imatb, imod, imode, ipas, ires
    integer(kind=8) :: itypfl, ivecc, j, jmod, lfsic, lfsvi, lfsvr
    integer(kind=8) :: nbp
    real(kind=8) :: aire, caj1, caj2, cd, cf0, ck, cocaj1
    real(kind=8) :: cocaj2, cokaj1, cokaj2, de, hmoy, p1, p2
    real(kind=8) :: phie, rug, val1, val2, visc, vr
    integer(kind=8), pointer :: tempo(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    depi = r8depi()
!
    fsic = typflu//'           .FSIC'
    call jeveuo(fsic, 'L', lfsic)
    itypfl = zi(lfsic)
    nom1 = '&&COEFAM.CDR2'
    nom2 = '&&COEFMO.COMPT'
    nom3 = '&&COEFMO.EXTR'
    nom4 = '&&COEFMO.VRZO'
    nom5 = '&&COEFMO.ALARM'
!
!
! --- 1.CAS D'UN FAISCEAU_TRANS
!
    if (itypfl .eq. 1) then
!
! ---   VECR1(I) EST LA VITESSE DU POINT I
! ---   VECR5(I) EST L'ABSCISSE CURVILIGNE POINT I.
! ---   VECR1(I) EST LA VITESSE REELLE DU POINT I.
! ---   VECR1(I+NBP) EST LA VITESSE MOYENNE DE LA ZONE D'EXCITATION
! ---   CONTENANT LE POINT I.
! ---   VECR1(1+2*NBP) EST LA VITESSE MOYENNE DE LA REUNION DE
!       TOUTES LES ZONES D'EXCITATION CONSTITUANT LE TUBE.
! ---   VECR2(I) EST LA MASSE VOLUMIQUE EXTERIEURE DU FLUIDE AU POINT I
! ---   VECR3(IMODE+I) EST LA DEFORMEE MODALE AU POINT I, POUR LA
!       COMPOSANTE DEFINIE DANS DEFI_FLUI_STRU
! ---   VECR4(1) EST LE DIAMETRE EXTERIEUR
! ---   VECI1(I) EST LE TYPE DE RESEAU DE LA ZONE CONTENANT LE POINT I
!
        nbp = indic
        de = vecr4(1)
        vr = vgap*depi/(pulsc*de*vecr1(2*nbp+1))
!
        imode = nbp*(nmode-1)
!
!   --- CALCUL DU COEFFICIENT D'AMORTISSEMENT AJOUTE ---
!
        fsvi = typflu//'           .FSVI'
        call jeveuo(fsvi, 'L', lfsvi)
        ipas = zi(lfsvi)
!
        call jeveuo('&&MDCONF.TEMPO', 'L', vi=tempo)
        nzone = tempo(1)
!
! ---   REINITIALISATION DES TABLEAUX :
!       NOM2 : NB POINTS HORS PLAGE (MIN, TOT, MAX), PAR ZONE
!       NOM3 : VREDMIN ET VREDMAX DES NOEUDS HORS PLAGE, PAR ZONE
!       NOM4 : VRMIN ET VRMAX ASSOCIES A CHAQUE ZONE
        call jedetr(nom2)
        call jedetr(nom3)
        call jedetr(nom4)
        call wkvect(nom2, 'V V I', 3*nzone, jcompt)
        call wkvect(nom3, 'V V R', 2*nzone, jextr)
        call wkvect(nom4, 'V V R', 2*nzone, jvrzo)
        do i = 1, 3*nzone
            zi(jcompt+i-1) = 0
        end do
!
        aire = 0.d0
        do i = 1, nbp
            ires = veci1(i)
!
            if ((ipas .eq. 1) .and. (ires .eq. 1003) .and. (vecr1(i) .ne. 0.d0)) then
                vred = vr*vecr1(nbp+i)
            else
                vred = vr*vecr1(i)
            end if
!
            call coefam(ipas, ires, vred, xsi0, cd)
!
            if (i .eq. 1) then
                val1 = cd*vecr2(1)*vecr1(1)*vecr3(imode+1)*vecr3(imode+1)
            else
                val2 = cd*vecr2(i)*vecr1(i)*vecr3(imode+i)*vecr3(imode+i)
                j = i-1
                aire = aire+(vecr5(i)-vecr5(j))*(val2+val1)/2.d0
                val1 = val2
            end if
!
! ---     DETERMINATION S'IL Y EN A UNE DE LA ZONE CONTENANT LE POINT I
            izone = 0
            do j = 1, nzone
                n1 = tempo(1+2*(j-1)+1)
                n2 = tempo(1+2*(j-1)+2)
                if ((i .ge. n1) .and. (i .le. n2)) then
                    izone = j
                    goto 12
                end if
            end do
12          continue
!
            if (izone .ne. 0) then
!
! ---       DETERMINATION DES VALEURS VRMIN ET VRMAX ASSOCIEES A LA ZONE
                if (i .eq. n1) then
                    call jeveuo(nom1, 'L', jvired)
                    vrmin = zr(jvired+0)
                    vrmax = zr(jvired+1)
                    zr(jvrzo+2*(izone-1)+0) = vrmin
                    zr(jvrzo+2*(izone-1)+1) = vrmax
                    zr(jextr+2*(izone-1)+0) = vrmin
                    zr(jextr+2*(izone-1)+1) = vrmax
                end if
!
! ---       DETERMINATION S'IL S'AGIT D'UN NOEUD HORS PLAGE
                if ((vred .lt. vrmin) .or. (vred .gt. vrmax)) then
                    call jeexin(nom5, iret)
                    if (iret .eq. 0) then
                        call wkvect(nom5, 'V V I', 1, jalarm)
                        call utmess('A', 'ALGELINE_27')
                    end if
!
! ---         INCREMENTATION DES NOMBRES DE NOEUDS HORS PLAGE :
!             ...+0 : NOMBRE DE NOEUDS AVEC VRED INFERIEURE A VRMIN
!             ...+1 : NOMBRE DE NOEUDS TOTAL HORS PLAGE
!             ...+2 : NOMBRE DE NOEUDS AVEC VRED SUPERIEURE A VRMAX
                    if (vred .lt. vrmin) then
                        zi(jcompt+3*(izone-1)+0) = zi(jcompt+3*(izone-1)+0)+1
                    else
                        zi(jcompt+3*(izone-1)+2) = zi(jcompt+3*(izone-1)+2)+1
                    end if
                    zi(jcompt+3*(izone-1)+1) = zi(jcompt+3*(izone-1)+1)+1
!
! ---         DETERMINATION DE LA PLUS PETITE ET DE LA PLUS GRANDE
!             VITESSE REDUITE HORS PLAGE
                    zr(jextr+2*(izone-1)+0) = min(vred, zr(jextr+2*( &
                                                           izone-1)+0))
                    zr(jextr+2*(izone-1)+1) = max(vred, zr(jextr+2*( &
                                                           izone-1)+1))
!
                end if
            end if
        end do
!
        call jelibe(nom1)
!
        xcf = de*vgap*aire/(2.d0*vecr1(2*nbp+1))
!
!   --- CALCUL DU COEFFICIENT DE RAIDEUR AJOUTE ---
!
        if (zrigi) then
            ires = veci1(1)
            if (ipas .eq. 1 .and. ires .eq. 1003 .and. vecr1(1) .ne. 0.d0) then
                vred = vr*vecr1(nbp+1)
            else
                vred = vr*vecr1(1)
            end if
            call coefra(ipas, ires, vred, xsi0, ck)
            val1 = ck*vecr2(1)*vecr1(1)*vecr1(1)*vecr3(imode+1)*vecr3(imode+1)
            aire = 0.d0
!
            do i = 2, nbp
                ires = veci1(i)
                if (ipas .eq. 1 .and. ires .eq. 1003 .and. vecr1(i) .ne. 0.d0) then
                    vred = vr*vecr1(nbp+i)
                else
                    vred = vr*vecr1(i)
                end if
                call coefra(ipas, ires, vred, xsi0, ck)
                val2 = ck*vecr2(i)*vecr1(i)*vecr1(i)*vecr3(imode+i)*vecr3(imode+i)
                aire = aire+(vecr5(i)-vecr5(i-1))*(val2+val1)/2.d0
                val1 = val2
            end do
            rkf = vgap*vgap*aire/(2.d0*vecr1(2*nbp+1)*vecr1(2*nbp+1))
            xkf = dcmplx(rkf, 0.d0)
        else
            xkf = dcmplx(0.d0, 0.d0)
        end if
!
!   --- CALCUL DU COEFFICIENT DE MASSE AJOUTE ---
!
        xmf = 0.d0
!
! --- 2.CAS D'UNE GRAPPE
!
    else if (itypfl .eq. 2) then
!
        phie = vecr4(1)
        vr = vgap*depi/(pulsc*phie)
        p1 = vecr3(2*(nmode-1)+1)
        p2 = vecr3(2*(nmode-1)+2)
!
!   --- CALCUL DU COEFFICIENT D'AMORTISSEMENT AJOUTE ---
!
        caj1 = vecr2(1)*vgap
        caj2 = vecr2(2)*vgap
        call cajgr2(indic, vr, cocaj1, cocaj2)
        xcf = caj1*cocaj1*p1+caj2*cocaj2*p2
!
!   --- CALCUL DU COEFFICIENT DE RAIDEUR AJOUTE ---
!
        if (zrigi) then
            kaj1 = vecr2(3)*vgap*vgap
            kaj2 = vecr2(4)*vgap*vgap
            call kajgr2(indic, vr, cokaj1, cokaj2)
            rkf = kaj1*cokaj1*p1+kaj2*cokaj2*p2
            xkf = dcmplx(rkf, 0.d0)
        else
            xkf = dcmplx(0.d0, 0.d0)
        end if
!
!   --- CALCUL DU COEFFICIENT DE MASSE AJOUTE ---
!
        xmf = 0.d0
!
! --- 3.CAS DE COQUE_COAX
!
    else if (itypfl .eq. 4) then
!
        imasse = indic
!
        hmoy = vecr4(1)
!
        fsvr = typflu//'           .FSVR'
        call jeveuo(fsvr, 'L', lfsvr)
        visc = zr(lfsvr+1)
        rug = zr(lfsvr+2)
!
        if (vgap .lt. 1.d-5) then
            cf0 = 0.d0
            mcf0 = 1.d0
            z = dcmplx(x(1), x(2))
            poids = z*z
        else
            call cfrott(visc, rug, hmoy, vgap, cf0, &
                        mcf0)
            poids = dcmplx(1.d0, 0.d0)
        end if
!
        if (imasse .eq. 0) then
!
            call bijmoc(vgap, vecr4, cf0, mcf0, zr(lfsvr), &
                        nmode, nmode, nbm, veci1, vecr2, &
                        vecr3, x(1), x(2), bii)
!
            biie = bii*poids
!
        else
!
            call wkvect('&&COEFMO.TEMP.MATB', 'V V C', nbm*nbm, imatb)
            call wkvect('&&COEFMO.TEMP.VECC', 'V V C', nbm, ivecc)
!
            call bmocca(vgap, vecr4, cf0, mcf0, zr(lfsvr), &
                        nbm, veci1, vecr2, vecr3, x(1), &
                        x(2), zc(imatb))
!
            do imod = 1, nbm
                do jmod = 1, nbm
                    zc(ivecc+imod-1) = zc(ivecc+imod-1)+zc(imatb+nbm*(jmod-1)+imod-1)*poids &
                                      &*vecr5(nbm*(nmode-1)+jmod)
                end do
            end do
!
            biie = dcmplx(0.d0, 0.d0)
            do imod = 1, nbm
                biie = biie+vecr5(nbm*(nmode-1)+imod)*zc(ivecc+imod-1)
            end do
!
        end if
!
!   --- CALCUL DU COEFFICIENT D'AMORTISSEMENT AJOUTE ---
!
        xcf = 0.d0
!
!   --- CALCUL DU COEFFICIENT DE RAIDEUR AJOUTE ---
!
        if (zrigi) then
            xkf = -biie
        else
            xkf = dcmplx(0.d0, 0.d0)
        end if
!
!   --- CALCUL DU COEFFICIENT DE MASSE AJOUTE ---
!
        if (zrigi) then
            xmf = -vecr1(nmode)
        else
            xmf = 0.d0
        end if
!
    end if
!
    call jedetr('&&COEFMO.TEMP.MATB')
    call jedetr('&&COEFMO.TEMP.VECC')
    call jedema()
end subroutine
