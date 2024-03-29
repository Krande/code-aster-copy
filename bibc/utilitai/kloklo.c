/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "aster.h"

/***************************************************************
 *
 *  Sous programme       : KLOKLO
 *
 *   Donne la date ( jour/mois/annee heure:minute:seconde )
 *   au moment d'appel de cette fonction.
 *
 *       Ce sous programme remplace le sous programme KLOAK de la
 *       bibliotheque IBM.
 *
 *   Cette fonction donne les memes informations que la commande date
 *       de UNIX.
 *
 *  Parametres d'entree  :
 *  Parametres de sortie :
 *
 *   date : tableau de 6 entiers dont les significations (FORTRAN)
 *          sont les suivantes
 *
 *       date(1) : Numero du jour (0-6) de la semaine compte
 *                         a partir du Lundi.
 *                         Ex si date(1) est 1 alors c'est le Mardi.
 *       date(2) : jour (1-28,29,30,31 suivant les mois)
 *       date(3) : mois (1-12)
 *       date(4) : annee ( qu'il faut additionner avec
 *             1900 pour avoir la vraie valeur )
 *       date(5) : heures (0-23)
 *       date(6) : minutes (0-59)
 *       date(7) : secondes (0-59)
 *       date(8) : Numero du jour dans l'annee (1-366)
 *       date(9) : Numero de la semaine dans l'annee (1-52)
 *
 *
 *  Variables globales   :
 *
 *     - les variables simples
 *
 *     - les tableaux
 *
 *  Sous programmes      :
 *
 *       time
 *       localtime
 *
 *  Exemple d'utilisation
 *
 *      program tkloklo
 *      character*8 jour(7)
 *      character*4 mois(12)
 *      integer tt(9)
 *      data jour/'Lundi', 'Mardi', 'Mercredi', 'Jeudi',
 *     &          'Vendredi', 'Samedi', 'Dimanche'/
 *      data mois/'Jan','Fev','Mars','Avr','Mai','Juin','Juil','Aout',
 *     &          'Sept', 'Oct','Nov','Dec'/
 *      call kloklo(tt)
 *      print *,jour(tt(1)+1), ', le ', tt(2), ' ', mois(tt(3)),' 19',tt(4),
 *     &                     ' a ', tt(5),':',tt(6),':',tt(7)
 *      print *,'Nous sommes au ', tt(8),
 *     &        'eme jour de l''annee, a la semaine ', tt(9)
 *      end
 *
 *  Auteur      : Albert Y
 *
 ****************************************************************/

void DEFP( KLOKLO, kloklo, ASTERINTEGER *date ) {
    time_t timval;
    struct tm *timeptr;

    time( &timval );
    timeptr = localtime( &timval );
    date[6] = ( ASTERINTEGER )( timeptr->tm_sec );
    date[5] = ( ASTERINTEGER )( timeptr->tm_min );
    date[4] = ( ASTERINTEGER )( timeptr->tm_hour );
    date[1] = ( ASTERINTEGER )( timeptr->tm_mday );
    date[2] = ( ASTERINTEGER )( timeptr->tm_mon + 1 );
    date[3] = ( ASTERINTEGER )( timeptr->tm_year + 1900 );
    date[0] = ( ASTERINTEGER )( timeptr->tm_wday == 0 ? 6 : timeptr->tm_wday - 1 );
    date[7] = ( ASTERINTEGER )( timeptr->tm_yday + 1 );
    date[8] = ( ASTERINTEGER )( ( date[7] + 6 ) / 7 );
}
