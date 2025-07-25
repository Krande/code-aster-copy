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
subroutine tfveri(nommcf, nbocc, itypfl)
    implicit none
!-----------------------------------------------------------------------
!     VERIFICATIONS DE PREMIER NIVEAU
!     APPELANT : OP0143 , OPERATEUR DEFI_FLUI_STRU
!-----------------------------------------------------------------------
!  IN   : NOMMCF : NOM DU MOT-CLE FACTEUR UTILISE
!  IN   : NBOCC  : NOMBRE D'OCCURENCES DU MOT-CLE FACTEUR UTILISE
!  IN   : ITYPFL : INDICE CARACTERISTIQUE DE LA CONFIGURATION ETUDIEE
!-----------------------------------------------------------------------
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/tfvegr.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: itypfl
    character(len=16) :: nommcf
! ----------------------------------------------------------------------
    integer(kind=8) :: count1, count2, count3, count4, count5
    integer(kind=8) :: ocgril
    character(len=2) :: carapa(4)
    character(len=3) :: ouinon
!      CHARACTER*9   TYPAS(2)
    real(kind=8) :: vect(3), valepa(4)
!
!      DATA TYPAS   /'CARRE_LIGN ','TRIA_LIGN'/
! ----------------------------------------------------------------------
!
!
!
! ----1.CAS D'UN FAISCEAU_TRANS
!       -----------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iangl, ibid, icapa, icara, icm, icmp, icoup
    integer(kind=8) :: ier2, igra2, ihy, ihz, iocc, ipas, ipesan
    integer(kind=8) :: ir, irayon, irho, irhoe, irhoi, irugo, itpas
    integer(kind=8) :: itres, ivapa, ivect, ivisc, jcoup, nbcoor, nbhy
    integer(kind=8) :: nbhz, nbocc, nbr, nbtub, nbtub2, nbyc, nbzc
    integer(kind=8) :: ncara, ncm, ncmp, ncoup, npas, nrhoe, nrhoi
    integer(kind=8) :: ntpas, ntres, ntypg
!-----------------------------------------------------------------------
    if (itypfl .eq. 1) then
! ---    VERIFICATION DE LA PRESENCE D AU MOINS UNE OCCURENCE DU
!        MOT-CLE COUPLAGE
        ncoup = 0
        do iocc = 1, nbocc
            call getvtx(nommcf, 'COUPLAGE', iocc=iocc, nbval=0, nbret=icoup)
            if (icoup .ne. 0) then
                ncoup = ncoup+1
                jcoup = iocc
            end if
        end do
        if (ncoup .eq. 0) then
            call utmess('E', 'MODELISA7_19')
            goto 999
        end if
        ncara = 0
        nrhoi = 0
        nrhoe = 0
        ncmp = 0
        do iocc = 1, nbocc
            call getvid(nommcf, 'CARA_ELEM', iocc=iocc, nbval=0, nbret=icara)
            if (icara .ne. 0) ncara = ncara+1
            call getvid(nommcf, 'PROF_RHO_F_INT', iocc=iocc, nbval=0, nbret=irhoi)
            if (irhoi .ne. 0) nrhoi = nrhoi+1
            call getvid(nommcf, 'PROF_RHO_F_EXT', iocc=iocc, nbval=0, nbret=irhoe)
            if (irhoe .ne. 0) nrhoe = nrhoe+1
            call getvtx(nommcf, 'NOM_CMP       ', iocc=iocc, nbval=0, nbret=icmp)
            if (icmp .ne. 0) ncmp = ncmp+1
        end do
!
        call getvtx(nommcf, 'COUPLAGE', iocc=jcoup, scal=ouinon, nbret=ibid)
!
! -------1.1.SI PRISE EN COMPTE DU COUPLAGE
!
        if (ouinon .eq. 'OUI') then
            ntpas = 0
            ntres = 0
            npas = 0
            ncm = 0
            do iocc = 1, nbocc
                call getvtx(nommcf, 'TYPE_PAS', iocc=iocc, nbval=0, nbret=itpas)
                if (itpas .ne. 0) then
                    ntpas = ntpas+1
                end if
                call getvis(nommcf, 'TYPE_RESEAU', iocc=iocc, nbval=0, nbret=itres)
                if (itres .ne. 0) then
                    ntres = ntres+1
                end if
                call getvr8(nommcf, 'PAS', iocc=iocc, nbval=0, nbret=ipas)
                if (ipas .ne. 0) then
!                  JPAS = IOCC
                    npas = npas+1
                end if
            end do
            if (ntpas .eq. 0 .or. ntres .ne. nbocc .or. npas .eq. 0) then
                call utmess('E', 'MODELISA7_20')
            end if
!
! ------1.2.SI NON PRISE EN COMPTE DU COUPLAGE
!
        else
            ncm = 0
            do iocc = 1, nbocc
                call getvr8(nommcf, 'COEF_MASS_AJOU', iocc=iocc, nbval=0, nbret=icm)
                if (icm .ne. 0) then
                    ncm = ncm+1
                end if
            end do
            if (ncm .eq. 0) then
                call utmess('E', 'MODELISA7_21')
            end if
        end if
!
! ------1.3.VERIFICATION DE LA PRESENCE  DES MOT-CLE DEVANT APPARAITRE
!       AU MOINS UNE FOIS DANS L UNE DES OCCURENCES DU MOT-CLE FACTEUR
!
        if (ncara .eq. 0) then
            call utmess('E', 'MODELISA7_22')
        end if
        if (nrhoi .eq. 0) then
            call utmess('E', 'MODELISA7_23')
        end if
        if (nrhoe .eq. 0) then
            call utmess('E', 'MODELISA7_24')
        end if
        if (ncara .eq. 0) then
            call utmess('E', 'MODELISA7_25')
        end if
!
! ----2.CAS D'UNE GRAPPE
!       ----------------
!
    else if (itypfl .eq. 2) then
!
        call getvtx(nommcf, 'COUPLAGE', iocc=1, scal=ouinon, nbret=ibid)
        if (ouinon .eq. 'OUI') then
            call getvtx(nommcf, 'GRAPPE_2', iocc=1, nbval=0, nbret=igra2)
            if (igra2 .eq. 0) then
                call utmess('E', 'MODELISA7_26')
            end if
        end if
!
!
! ----3.CAS D'UN FAISCEAU_AXIAL
!       -----------------------
!
    else if (itypfl .eq. 3) then
!
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        ocgril = 0
!
        do iocc = 1, nbocc
!
! --------3.1.SI PLUSIEURS OCCURENCES <RAYON_TUBE> ET <COOR_TUBE>
! --------    OBLIGATOIRES A CHAQUE OCCURENCE
! --------    VERIFICATION DES DONNEES POUR <COOR_TUBE>
!
            call getvr8(nommcf, 'RAYON_TUBE', iocc=iocc, nbval=0, nbret=irayon)
            if (irayon .eq. 0) then
                if (nbocc .gt. 1) then
                    call utmess('E', 'MODELISA7_27')
                end if
            else
                call getvr8(nommcf, 'COOR_TUBE', iocc=iocc, nbval=0, nbret=nbcoor)
                nbcoor = abs(nbcoor)
                nbtub = int(nbcoor/2)
                nbtub2 = 2*nbtub
                if (nbtub2 .ne. nbcoor) then
                    call utmess('E', 'MODELISA7_28')
                end if
            end if
!
! --------3.2.INCREMENTATION DU COMPTEUR POUR <VECT_X> ET VERIFICATION
! --------    DES DONNEES SI PRESENCE
!
            call getvr8(nommcf, 'VECT_X', iocc=iocc, nbval=0, nbret=ivect)
            if (ivect .ne. 0) then
                count1 = count1+1
                if (abs(ivect) .ne. 3) then
                    call utmess('E', 'MODELISA7_29')
                else
                    ier2 = 0
                    call getvr8(nommcf, 'VECT_X', iocc=iocc, nbval=3, vect=vect(1), &
                                nbret=ibid)
                    if (vect(1) .eq. 1.d0) then
                        if (vect(2) .ne. 0.d0 .or. vect(3) .ne. 0.d0) ier2 = 1
                    else if (vect(2) .eq. 1.d0) then
                        if (vect(1) .ne. 0.d0 .or. vect(3) .ne. 0.d0) ier2 = 1
                    else if (vect(3) .eq. 1.d0) then
                        if (vect(1) .ne. 0.d0 .or. vect(2) .ne. 0.d0) ier2 = 1
                    else
                        ier2 = 1
                    end if
                    if (ier2 .eq. 1) then
                        call utmess('E', 'MODELISA7_30')
                    end if
                end if
            end if
!
! --------3.3.INCREMENTATION DES COMPTEURS POUR <PROF_RHO_FLUI>,
! --------    <PROF_VISC_CINE> ET <RUGO_TUBE>
!
            call getvid(nommcf, 'PROF_RHO_FLUI', iocc=iocc, nbval=0, nbret=irho)
            if (irho .ne. 0) count2 = count2+1
!
            call getvid(nommcf, 'PROF_VISC_CINE', iocc=iocc, nbval=0, nbret=ivisc)
            if (ivisc .ne. 0) count3 = count3+1
!
            call getvr8(nommcf, 'RUGO_TUBE', iocc=iocc, nbval=0, nbret=irugo)
            if (irugo .ne. 0) count4 = count4+1
!
! --------3.4.VERIFICATION DES DONNEES POUR <PESANTEUR> SI PRESENCE
!
            call getvr8(nommcf, 'PESANTEUR', iocc=iocc, nbval=0, nbret=ipesan)
            ipesan = abs(ipesan)
            if (ipesan .ne. 0 .and. ipesan .ne. 4) then
                call utmess('E', 'MODELISA7_31')
            end if
!
! --------3.5.INCREMENTATION DU COMPTEUR POUR <CARA_PAROI>
! --------    VERIFICATION DES DONNEES POUR <CARA_PAROI>, <VALE_PAROI>
! --------    ET <ANGL_VRIL> SI PRESENCE
!
            call getvtx(nommcf, 'CARA_PAROI', iocc=iocc, nbval=0, nbret=icapa)
            icapa = abs(icapa)
            if (icapa .ne. 0) then
                count5 = count5+1
                if (icapa .ne. 3 .and. icapa .ne. 4) then
                    call utmess('E', 'MODELISA7_32')
                else
                    call getvr8(nommcf, 'VALE_PAROI', iocc=iocc, nbval=0, nbret=ivapa)
                    ivapa = abs(ivapa)
                    if (ivapa .ne. icapa) then
                        call utmess('E', 'MODELISA7_33')
                    else
                        call getvtx(nommcf, 'CARA_PAROI', iocc=iocc, nbval=icapa, vect=carapa(1), &
                                    nbret=ibid)
                        call getvr8(nommcf, 'VALE_PAROI', iocc=iocc, nbval=ivapa, vect=valepa(1), &
                                    nbret=ibid)
                        nbyc = 0
                        nbzc = 0
                        nbr = 0
                        nbhy = 0
                        nbhz = 0
                        if (icapa .eq. 3) then
                            do icara = 1, icapa
                                if (carapa(icara) .eq. 'YC') nbyc = nbyc+1
                                if (carapa(icara) .eq. 'ZC') nbzc = nbzc+1
                                if (carapa(icara) (1:1) .eq. 'R') then
                                    nbr = nbr+1
                                    ir = icara
                                end if
                            end do
                            if (nbyc .ne. 1 .or. nbzc .ne. 1 .or. nbr .ne. 1) then
                                call utmess('E', 'MODELISA7_34')
                            else if (valepa(ir) .le. 0.d0) then
                                call utmess('E', 'MODELISA7_35')
                            end if
                        else
                            do icara = 1, icapa
                                if (carapa(icara) .eq. 'YC') nbyc = nbyc+1
                                if (carapa(icara) .eq. 'ZC') nbzc = nbzc+1
                                if (carapa(icara) .eq. 'HY') then
                                    nbhy = nbhy+1
                                    ihy = icara
                                end if
                                if (carapa(icara) .eq. 'HZ') then
                                    nbhz = nbhz+1
                                    ihz = icara
                                end if
                            end do
                            if (nbyc .ne. 1 .or. nbzc .ne. 1 .or. nbhy .ne. 1 .or. nbhz &
                                .ne. 1) then
                                call utmess('E', 'MODELISA7_36')
                            else if (valepa(ihy) .le. 0.d0 .or. valepa( &
                                     ihz) .le. 0.d0) then
                                call utmess('E', 'MODELISA7_37')
                            else
                                call getvr8(nommcf, 'ANGL_VRIL', iocc=iocc, nbval=0, nbret=iangl)
                                if (iangl .eq. 0) then
                                    call utmess('E', 'MODELISA7_38')
                                end if
                            end if
                        end if
                    end if
                end if
            end if
!
! --------3.6.DETECTION DE LA DERNIERE OCCURENCE POUR LAQUELLE LES
! --------    OPERANDES ASSOCIEES AUX CARACTERISTIQUES DES GRILLES
! --------    SONT PRESENTES
!
            call getvr8(nommcf, 'LONG_TYPG', iocc=iocc, nbval=0, nbret=ntypg)
            if (ntypg .ne. 0) then
                ocgril = iocc
            end if
!
        end do
!
! ------3.7.VERIFICATION DES COMPTEURS
!
        if (count1 .eq. 0) then
            call utmess('E', 'MODELISA7_39')
        else if (count2 .eq. 0) then
            call utmess('E', 'MODELISA7_40')
        else if (count3 .eq. 0) then
            call utmess('E', 'MODELISA7_41')
        else if (count4 .eq. 0) then
            call utmess('E', 'MODELISA7_42')
        else if (count5 .eq. 0) then
            call utmess('E', 'MODELISA7_43')
        end if
!
! ------3.8.VERIFICATION DES DONNEES CARACTERISTIQUES DES GRILLES
!
        if (ocgril .ne. 0) then
            call tfvegr(nommcf, ocgril)
        end if
!
!
! ----4.CAS DE COQUE_COAX
!       -----------------
!
    else
!
        call getvr8(nommcf, 'VECT_X', iocc=1, nbval=0, nbret=ivect)
        if (abs(ivect) .ne. 3) then
            call utmess('E', 'MODELISA7_44')
        else
            ier2 = 0
            call getvr8(nommcf, 'VECT_X', iocc=1, nbval=3, vect=vect(1), &
                        nbret=ibid)
            if (vect(1) .eq. 1.d0) then
                if (vect(2) .ne. 0.d0 .or. vect(3) .ne. 0.d0) ier2 = 1
            else if (vect(2) .eq. 1.d0) then
                if (vect(1) .ne. 0.d0 .or. vect(3) .ne. 0.d0) ier2 = 1
            else if (vect(3) .eq. 1.d0) then
                if (vect(1) .ne. 0.d0 .or. vect(2) .ne. 0.d0) ier2 = 1
            else
                ier2 = 1
            end if
            if (ier2 .eq. 1) then
                call utmess('E', 'MODELISA7_45')
            end if
        end if
!
    end if
!
999 continue
!
end subroutine
