# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------
import numpy


def loi_de_kelvin_hr(pc, Temp, rho_liquide, R=8.314, Mmolaire=18.01528e-3):
    hr = numpy.exp(pc * Mmolaire / (rho_liquide * R * Temp))
    return hr


def diffusion_vapeur_air(Temp):
    """
    Diffusion de la vapeur d'eau dans l'air libre'
    """
    D0 = 0.217 * 1e-4 * ((Temp / 273.15)) ** 1.88
    return D0


def modele_de_Rankine_pvs(Temp, alpha_rankine=13.7, beta_rankine=5120, Pa=101325):
    """
    Equation de Rankine pour le calcul de la pression de vapeur saturante
    """
    Pvs = Pa * numpy.exp(alpha_rankine - (beta_rankine / Temp))
    return Pvs


def tension_superficielle(Temp):
    """
    Tension superficielle de l'eau'
    """
    gamma = 0.1558 * (1 - (Temp / 647.1)) ** 1.26
    return gamma


def viscosite_liquide(Temp):
    """
    Viscosite liquide
    """
    nul = 0.6612 * (Temp - 229) ** (-1.562)
    return nul


def densite_liquide(Temp):
    """
    Masse volumique de l'eau liquide
    """
    rho_liquide = 314.4 + 685.6 * (1 - ((Temp - 273.15) / 374.14) ** (1 / 0.55)) ** 0.55
    return rho_liquide


def permeabilite_liquide(K0, Temp, Ea_R, T0):
    """
    Permeabilite_liquide
    """
    Kl = K0 * numpy.exp(
        numpy.exp((Temp - T0) / Ea_R) - 1.0
    )  # (K0*Temp/T0)*numpy.exp(-Ea_R*((1/Temp)-(1/T0))) # [m2]
    return Kl


def isotherme_hr_T(C, alpha, beta, Ad, T0, poro, Temp, R=8.314, Mmolaire=18.01528e-3):
    """
    Isotherme Leverett(C,T)
    C : Concentration en eau
    Temp : Température en Kelvin
    """
    pc = pression_capillaire(C, alpha, beta, Ad, T0, poro, Temp, R, Mmolaire)

    rho_liquide = densite_liquide(Temp)
    HR = loi_de_kelvin_hr(pc, Temp, rho_liquide)
    return HR


def pression_capillaire(C, alpha, beta, Ad, T0, poro, Temp, R=8.314, Mmolaire=18.01528e-3):
    """
    Pression capillaire(C,T)
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    gamma0 = tension_superficielle(T0)
    gamma = tension_superficielle(Temp)
    a = densite_liquide(T0) * R * T0 / (alpha * Mmolaire)
    K0_KT = (10 ** (Ad * (2 * 10**-3 * (Temp - T0) - 1e-6 * (Temp - T0) ** 2))) ** (-1)
    pc = (
        -a
        * (((C / (poro * 1000)) ** (-1.0 / beta) - 1.0) ** (1.0 - beta))
        * (gamma0 / gamma)
        * numpy.sqrt(K0_KT)
    )
    return pc


def permeabilte_relative_VGM(C, p, beta, poro):
    """
    Permeabilite relative de Mualem-Van-guenuchten
    """
    krl = ((C / (poro * 1000)) ** p) * (
        1.0 - (1.0 - (C / (poro * 1000)) ** (1 / beta)) ** beta
    ) ** 2
    return krl


def derivee_isotherme_pc(C, alpha, beta, Ad, T0, poro, Temp, R=8.314, Mmolaire=18.01528e-3):
    """
    Derivee isotherme de desorption de Van-guenuchten
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    gamma0 = tension_superficielle(T0)
    gamma = tension_superficielle(Temp)
    a = densite_liquide(T0) * R * T0 / (alpha * Mmolaire)
    KT_K0 = 10 ** (Ad * (2 * 10**-3 * (Temp - T0) - 1e-6 * (Temp - T0) ** 2))
    dpc = (
        ((gamma / gamma0) * ((KT_K0) ** 0.5) * a * (-beta + 1.0) / (beta))
        * ((C / (poro * 1000)) ** (-(1.0 + beta) / beta))
        * (-1.0 + (C / (poro * 1000)) ** (-1.0 / beta)) ** (-beta)
    )
    return dpc


def facteur_de_resistance(C, poro, a, b):
    """
    Facteur de resistance de Millington
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    R = ((poro) ** a) * (1.0 - C / (poro * 1000)) ** b
    return R


def coefficient_de_diffusion_fick(C, poro, a, b, Temp):
    """
    Coefficient de diffusion effectif du béton
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    D0 = diffusion_vapeur_air(Temp)
    R = facteur_de_resistance(C, poro, a, b)
    D = D0 * R
    return D


def coefficient_de_diffusion_vapeur(
    C, poro, alpha, beta, a, b, Ad, T0, Temp, Mmolaire=18.01528e-3, R=8.314
):
    """
    Coefficient de diffusion de la vapeur d'eau
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    Temp = Temp + 273.15
    rho_liquide = densite_liquide(Temp)
    D = coefficient_de_diffusion_fick(C, poro, a, b, Temp)
    dpc = derivee_isotherme_pc(C, alpha, beta, Ad, T0, poro, Temp, R, Mmolaire)
    HR = isotherme_hr_T(C, alpha, beta, Ad, T0, poro, Temp, R, Mmolaire)
    Pv = HR * modele_de_Rankine_pvs(Temp)
    Dv = (D * Pv * dpc * (Mmolaire / (R * Temp)) ** 2) / (poro * rho_liquide**2)
    return Dv


def coefficient_de_diffusion(
    C, K0, poro, alpha, beta, p, Ad, T0, Temp, Ea_R, Mmolaire=18.01528e-3, R=8.314
):
    """
    Coefficient de diffusion de Richards
    """
    C = numpy.minimum(C, 0.999 * poro * 1000)
    Temp = Temp + 273.15
    Kl = permeabilite_liquide(K0, Temp, Ea_R, T0)
    nul = viscosite_liquide(Temp)
    krl = permeabilte_relative_VGM(C, p, beta, poro)
    dpc = derivee_isotherme_pc(C, alpha, beta, Ad, T0, poro, Temp, R, Mmolaire)
    Dl = (Kl / (poro * nul)) * krl * dpc
    return Dl
