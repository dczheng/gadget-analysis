#!/usr/bin/env python3

from astropy.cosmology import LambdaCDM
import astropy.constants  as ac

m_p      =                       ac.m_p.cgs.value
m_e      =                       ac.m_e.cgs.value
k_b      =                       ac.k_B.cgs.value
c        =                       ac.c.cgs.value
sigma_t  =                       ac.sigma_T.cgs.value
G        =                       ac.G.cgs.value
pc       =                       ac.pc.cgs.value
Kpc      =                       pc * 13
Mpc      =                       pc * 1e6
c2       =                       c * c
m_ec     =                       m_e * c
m_ec2    =                       m_e * c2
h        =                       ac.h.cgs.value
hbar     =                       ac.hbar.cgs.value
alpha    =                       ac.alpha.value
Msun     =                       ac.M_sun.cgs.value

H0       =                       3.2407789e-18  # /* in h/sec */
Bcmb     =                       (3.24e-6) # // gauss
Xh       =                       0.760
elec_frac =                      (1+Xh) / (2*Xh)
Gamma    =                       5.0 / 3.0
e        =                       4.8032e-10
e2       =                       e**2
e3       =                       e**3
Jy       =                       1e-23 # erg s^-1 cm^-2 Hz^-1
mJy       =                      Jy * 1e-3

yr       =                       3.155e7
Myr      =                       yr * 1e6
Gyr      =                       Myr * 1e6

Omega0      =                    0.302
OmegaLambda =                    0.698
OmegaBaryon =                    0.04751
HubbleParam =                    0.68

gadget_length_in_cm            =              Kpc  / HubbleParam
gadget_mass_in_g               =              Msun / HubbleParam
gadget_velocity_in_cm_per_s    =              1e5  / HubbleParam
gadget_time_in_s               =              gadget_length_in_cm / gadget_velocity_in_cm_per_s
gadget_energy_in_erg           =              gadget_mass_in_g * gadget_length_in_cm**2 / (gadget_time_in_s**2)

cosmology = LambdaCDM( H0 = HubbleParam*100, Om0=Omega0, Ode0 = OmegaLambda )

D_c = lambda z: cosmology.comoving_distance( z ).cgs.value
D_a = lambda z: cosmology.angular_diameter_distance( z ).cgs.value
D_l = lambda z: cosmology.luminosity_distance( z ).cgs.value

rho_crit = lambda z: cosmology.critical_density( z ).cgs.value
rho_crit0  = cosmology.critical_density0.cgs.value

rho_bar_crit = lambda z: rho_crit( z ) * OmegaBaryon
rho_bar_crit0 = rho_crit0 * OmegaBaryon

n_p_crit = lambda z: rho_bar_crit(z) / m_p * Xh
n_p_crit0 = rho_bar_crit0 / m_p * Xh


n_e_crit =  lambda z: n_p_crit( z ) * elec_frac
n_e_crit0 = n_p_crit0 * elec_frac
