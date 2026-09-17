# Original code written by Taras Verkholyak,
# Institute for Condensed Matter Physics,
# National Academy of Sciences of Ukraine
#
# Adapted and modified by Vishnu Prasad Kurupath
# (https://vishnu-prasad-kurupath.github.io/)
#
# Merged from the formerly-duplicated mpore_gcmc_serial.py /
# mpore_gcmc_serial_c.py into a single library used by run_gcmc.py
# (plain voltage loop, or --auto-resume for fine-grained checkpointing).
# Ion confinement (volume_free/box_bound/random_position) uses the
# accessible radius (Rex), not the wall-atom/carbon-centre radius (Rc) --
# one of the two pre-merge files got this wrong; see analyze_gcmc.py for
# the (correctly Rex-based) charge/capacitance/energy normalization.

import numpy as np
import math
import mplib_ctypes as mplib
from scipy import interpolate
import pickle

mplib_cyl_u2_vec = np.vectorize(mplib.cyl_u2, otypes=[float])
mplib_slit_u2_vec = np.vectorize(mplib.slit_u2, otypes=[float])


def resolve_pore_width(pore_width_accessible=None, pore_width_nominal=None,
                        wall_atom_radius=0.0, tol=1e-6):
    """Resolve (pore_width_accessible, pore_width_nominal) from user input.

    Accepts either or both of pore_width_accessible (diameter for a
    cylindrical pore, gap width for a slit pore -- the width available to
    ion centres) and pore_width_nominal (to the wall-atom/carbon centres),
    related by pore_width_accessible = pore_width_nominal -
    2*wall_atom_radius for both pore types (for a slit this is two
    physical walls each contributing wall_atom_radius; for a cylinder it
    is a single wall but pore_width is a diameter, so the single-radius
    offset doubles when expressed as a width). If both are given, they
    must be consistent -- this is exactly the nominal-vs-accessible mixup
    that caused a real confinement bug in the pre-merge mpore_gcmc_serial.py
    (fixed above). If only one is given, the other is derived.
    """
    if pore_width_accessible is None and pore_width_nominal is None:
        raise ValueError(
            "must provide at least one of pore_width_accessible or "
            "pore_width_nominal")
    if pore_width_accessible is not None and pore_width_nominal is not None:
        expected = pore_width_nominal - 2.0 * wall_atom_radius
        if abs(expected - pore_width_accessible) > tol:
            raise ValueError(
                f"inconsistent pore geometry: pore_width_nominal="
                f"{pore_width_nominal} - 2*wall_atom_radius="
                f"{wall_atom_radius} = {expected}, but "
                f"pore_width_accessible={pore_width_accessible} was also "
                f"given (differs by {expected - pore_width_accessible:+.6f} A)")
        return pore_width_accessible, pore_width_nominal
    if pore_width_accessible is not None:
        return pore_width_accessible, pore_width_accessible + 2.0 * wall_atom_radius
    return pore_width_nominal - 2.0 * wall_atom_radius, pore_width_nominal


class pore(object):

    def __init__(self, excl_width=0.0, length=0.0, eshift=0, wall_atom_radius=0.):

        if excl_width == 0.0:
            raise NameError("Pore width cannot be zero")

        self.width_ex = excl_width
        self.width_c = excl_width + 2. * wall_atom_radius
        self.width_el = self.width_c - 2. * eshift
        self.eshift = eshift

        self.Rc = self.width_c / 2.
        self.Rel = self.width_el / 2.
        self.Rex = self.width_ex / 2.

        if length == 0.0:
            raise NameError("Pore latteral size cannot be zero")

        self.length = length
        self.length_2 = length / 2.
        self.length_x = length
        self.length_y = length
        self.length_x_2 = self.length_x / 2.
        self.length_y_2 = self.length_y / 2.


class cylinder(pore):

    def __init__(self, excl_width=0, length=0, eshift=0, wall_atom_radius=0,
                 aion=None, k_B=None, d_ion=0.0, e_charge=0.0, temp=0.0):
        super().__init__(excl_width, length, eshift, wall_atom_radius)

        self.rho_factor = d_ion / self.length  # for 1D?? 3D vol_ion/vol_tube?

        r_data, U1_data = mplib.cyl_u1_calc_array(self.Rel, self.Rex - np.min(aion), 200)
        self.u1_interpol = interpolate.interp1d(r_data, U1_data)
        self.volume_free = self.length * math.pi * (self.Rex - aion)**2
        self.N_half_packing = int(self.length / (2.0 * d_ion))

    def u1(self, r):
        return self.u1_interpol(np.sqrt(r[0]**2 + r[1]**2))

    def u2(self, r1, r2):
        dz = np.abs(r1[:, 2] - r2[2])
        dz_wrap_mask = dz > self.length_2
        dz[dz_wrap_mask] = self.length - dz[dz_wrap_mask]
        rxy1 = np.sqrt(r1[:, 0]**2 + r1[:, 1]**2)
        rxy2 = np.sqrt(r2[0]**2 + r2[1]**2)
        phi = np.arctan2(r1[:, 1], r1[:, 0]) - np.arctan2(r2[1], r2[0])
        return mplib_cyl_u2_vec(rxy1, rxy2, phi, dz, self.Rel)

    def overlap(self, r_pos, n_comp, ion_coord, ion_type, hc_dist2):
        ''' checks if the hard-core condition is satisfied for inserted ion '''
        dz = np.abs(ion_coord[:, 2] - r_pos[2])
        dz_wrap_mask = dz > self.length_2
        dz[dz_wrap_mask] = self.length - dz[dz_wrap_mask]
        ion_rad_sum = hc_dist2[ion_type, n_comp]
        overlap_mask = (ion_coord[:, 0] - r_pos[0])**2 + (ion_coord[:, 1] - r_pos[1])**2 + dz**2 < ion_rad_sum
        if True in overlap_mask:
            return True
        return False

    def box_bound(self, r, i_comp, aion):
        if (r[0]**2 + r[1]**2 < ((self.Rex - aion)[i_comp])**2):
            return True
        else:
            return False

    def random_displacement(self, ion_coord, ion_type, max_step, i_ion, aion):
        change = np.array([max_step[0, ion_type[i_ion]] * (np.random.rand() - 0.5),
                            max_step[1, ion_type[i_ion]] * (np.random.rand() - 0.5),
                            max_step[2, ion_type[i_ion]] * (np.random.rand() - 0.5)])
        new_coord = ion_coord[i_ion] + change
        while self.box_bound(new_coord, ion_type[i_ion], aion) == False:
            change[0] = max_step[0, ion_type[i_ion]] * (np.random.rand() - 0.5)
            change[1] = max_step[1, ion_type[i_ion]] * (np.random.rand() - 0.5)
            new_coord = ion_coord[i_ion] + change
        new_coord[2] = new_coord[2] % self.length
        return new_coord

    def random_position(self, i_comp, aion):
        lims = self.Rex - aion
        new_pos = np.array([2 * lims[i_comp] * (np.random.rand() - 0.5),
                             2 * lims[i_comp] * (np.random.rand() - 0.5),
                             self.length * (np.random.rand())])
        while new_pos[0]**2 + new_pos[1]**2 > lims[i_comp]**2:
            new_pos[0] = 2 * lims[i_comp] * (np.random.rand() - 0.5)
            new_pos[1] = 2 * lims[i_comp] * (np.random.rand() - 0.5)
        return new_pos


class slit(pore):

    def __init__(self, excl_width=0, length=0, eshift=0, wall_atom_radius=0,
                 aion=None, k_B=None, d_ion=0.0, e_charge=0.0, temp=0.0):
        super().__init__(excl_width, length, eshift, wall_atom_radius)

        self.rho_factor = (d_ion / self.length)**2

        r_data, U1_data = mplib.slit_u1_calc_array(self.width_el, min(aion), 200)
        r_data = r_data - self.Rel
        self.u1_interpol = interpolate.interp1d(r_data, U1_data)
        self.volume_free = self.length_x * self.length_y * (self.Rex - aion) * 2
        self.N_half_packing = (int(self.length / (2 * d_ion)))**2

    def u1(self, r):
        return self.u1_interpol(r[2])

    def u2(self, r1, r2):
        dx = np.abs(r1[:, 0] - r2[0])
        dx_wrap_mask = dx > self.length_x_2
        dx[dx_wrap_mask] = self.length_x - dx[dx_wrap_mask]
        dy = np.abs(r1[:, 1] - r2[1])
        dy_wrap_mask = dy > self.length_y_2
        dy[dy_wrap_mask] = self.length_y - dy[dy_wrap_mask]
        R = np.sqrt(dx**2 + dy**2)
        z1 = (r1[:, 2] + self.Rel)
        z2 = (r2[2] + self.Rel)
        return mplib_slit_u2_vec(z1, z2, R, self.width_el)

    def overlap(self, r_pos, n_comp, ion_coord, ion_type, hc_dist2):
        ''' checks if the hard-core condition is satisfied for inserted ion '''
        if ion_type.size > 0:
            dx = np.abs(ion_coord[:, 0] - r_pos[0])
            dy = np.abs(ion_coord[:, 1] - r_pos[1])
            dz = np.abs(ion_coord[:, 2] - r_pos[2])
            dx_wrap_mask = dx > self.length_x_2
            dx[dx_wrap_mask] = self.length_x - dx[dx_wrap_mask]
            dy_wrap_mask = dy > self.length_y_2
            dy[dy_wrap_mask] = self.length_y - dy[dy_wrap_mask]
            ion_rad_sum = hc_dist2[ion_type, n_comp]
            overlap_mask = dx**2 + dy**2 + dz**2 < ion_rad_sum
            if True in overlap_mask:
                return True
        return False

    def box_bound(self, r, i_comp, aion):
        if (np.abs(r[2]) < (self.Rex - aion)[i_comp]):
            return True
        else:
            return False

    def random_displacement(self, ion_coord, ion_type, max_step, i_ion, aion):
        ''' n_ion is the number of ion to be moved '''
        change = np.array([max_step[0, ion_type[i_ion]] * (np.random.rand() - 0.5),
                            max_step[1, ion_type[i_ion]] * (np.random.rand() - 0.5),
                            max_step[2, ion_type[i_ion]] * (np.random.rand() - 0.5)])
        new_coord = ion_coord[i_ion] + change
        while self.box_bound(new_coord, ion_type[i_ion], aion) == False:
            change[2] = max_step[2, ion_type[i_ion]] * (np.random.rand() - 0.5)
            new_coord = ion_coord[i_ion] + change
        new_coord[0] = new_coord[0] % self.length_x
        new_coord[1] = new_coord[1] % self.length_y
        return new_coord

    def random_position(self, i_comp, aion):
        ''' generates a random position of the ion in the slit '''
        new_pos = np.array([self.length_x * np.random.rand(),
                             self.length_y * np.random.rand(),
                             (self.Rex - aion)[i_comp] * 2. * (np.random.rand() - 0.5)])
        return new_pos


class mc_settings(object):

    def __init__(self, Ntherm=1000000, Nsim=1000000):
        self.n_therm = Ntherm
        self.n_sim = Nsim
        # Continuation/resume bookkeeping -- only meaningfully driven when
        # run_gcmc.py is run with --auto-resume; harmless no-ops (stay at
        # 0) otherwise.
        self.c_therm = 0
        self.c_sim = 0
        self.complete = 0

        self.p_trans = 0.3
        self.p_widom = 0.6
        self.p_MTSW = 1.0

    def set_prob(self, trans, widom):
        self.p_trans = trans
        self.p_widom = widom
        self.p_MTSW = 1.

    def print(self, filename):
        f = open(filename, 'a')
        f.write("# Number of thermalisation steps {0}\n".format(self.n_therm))
        f.write("# Number of simulation steps {0}\n".format(self.n_sim))
        f.write("# Translational move: {0}\n".format(self.p_trans))
        f.write("# Widom insertion/deletion move: {0}\n".format(self.p_widom))
        f.write("# MTSW move: {0}\n".format(self.p_MTSW))
        f.close()


class output(object):

    def __init__(self) -> None:
        pass

    def dump_count(self, ion_type, n_comp, prefix, voltage):
        counts, _ = np.histogram(ion_type, bins=np.arange(n_comp + 1) - 0.5)
        with open(prefix + '_' + str(voltage) + '.count', 'a+') as f:
            np.savetxt(f, counts.reshape(-1, n_comp))

    def dump_coords(self, ion_coords, ion_types, i_step, prefix, voltage):
        with open(prefix + '_' + str(voltage) + '.coords', 'a+') as f:
            f.write('step:' + str(i_step) + '\n')
            np.savetxt(f, np.c_[ion_types, ion_coords])
        f.close()

    def print(self, step, ion_type, n_comp):
        if ion_type.size > 0:
            counts, _ = np.histogram(ion_type, bins=np.arange(n_comp + 1) - 0.5)
            print('{}, {}'.format(step, counts))
        else:
            print('{}, None'.format(step))


class state(object):

    def __init__(self, T=300., epsr=2.5, aion=np.array([2.5, 2.5]),
                 q=np.array([-1, 1]), ptype=None,
                 pwidth=0.0, plength=0.0, eshift=0.0, wall_atom_radius=0.0,
                 n_therm=1000000, n_sim=10000000, p_trans=0.5, p_widom=1.0):
        self.k_B = 1.38064852e-23  # J/K
        self.e_charge = 1.6e-19  # C
        self.k_e = 8.9875517923e9  # N*m^2/C^2,  Coulomb constant
        self.factor_J_to_eV = 6.242e18
        self.factor_K_to_eV = 8.61732814974056E-05

        self.ptype = ptype
        if ptype is None:
            raise NameError("Pore cannot be None")

        self.temp = T
        self.eps_r = epsr
        self.factor_kBT_to_eV = T * self.k_B * self.factor_J_to_eV
        self.U_conv_coeff = 332.0636 * 503.2166 / self.eps_r / self.temp

        if len(aion) != len(q):
            raise NameError("Array sizes of the radii and charges are not equal")

        self.aion = aion
        self.n_comp = len(aion)
        self.d_ion = 2. * np.max(aion)

        self.q_comp = q

        self.hc_dist2 = (self.aion[:, np.newaxis] + self.aion[np.newaxis, :])**2

        if self.ptype == 'cyl':
            self.pore = cylinder(pwidth, plength, eshift, wall_atom_radius,
                                  self.aion, self.k_B, self.d_ion, self.e_charge, self.temp)
        elif self.ptype == 'slit':
            self.pore = slit(pwidth, plength, eshift, wall_atom_radius,
                              self.aion, self.k_B, self.d_ion, self.e_charge, self.temp)
        else:
            raise NameError(f"Unknown pore type '{self.ptype}': must be 'cyl' or 'slit'")

        self.mcparams = mc_settings(n_therm, n_sim)
        self.mcparams.set_prob(p_trans, p_widom)

        self.r_xy = self.pore.Rex - self.aion  # change pore functions to accept aion !!
        self.r_xy2 = self.r_xy**2
        self.d_xy = 2. * self.r_xy
        self.d_z = 0.2 * 2 * self.aion
        self.step = np.vstack([self.d_xy, self.d_xy, self.d_z])

        self.n_tot = self.pore.N_half_packing
        self.ion_type = np.random.randint(self.n_comp, size=self.n_tot)
        self.ion_coord = np.empty((self.n_tot, 3))

        for ion_num in np.arange(self.n_tot):
            r_position = self.pore.random_position(self.ion_type[ion_num], self.aion)
            while self.pore.overlap(r_position, self.ion_type[ion_num],
                                     self.ion_coord[:ion_num], self.ion_type[:ion_num],
                                     self.hc_dist2) == True:
                r_position = self.pore.random_position(self.ion_type[ion_num], self.aion)
            self.ion_coord[ion_num, :] = r_position

        self.output = output()

        # Only meaningfully driven by run_gcmc.py's --auto-resume loop.
        self.c_voltage = 0.0

    def kBT2eV(self, T):
        return T * self.k_B * self.factor_J_to_eV


class mcfunctions(object):

    def __init__(self, state):
        self.state = state

    def d_energy_displace(self, i_ion, r_pos):
        ''' the energy difference when moving n-th ion on a new position z '''
        if self.state.pore.overlap(r_pos, self.state.ion_type[i_ion],
                                    np.delete(self.state.ion_coord, i_ion, axis=0),
                                    np.delete(self.state.ion_type, i_ion, axis=0),
                                    self.state.hc_dist2) == True:
            return 0., True
        d_en = (self.state.q_comp[self.state.ion_type[i_ion]]**2) * \
            (self.state.pore.u1(r_pos) -
             self.state.pore.u1(self.state.ion_coord[i_ion]))
        d_en += np.sum(self.state.q_comp[np.delete(self.state.ion_type, i_ion, axis=0)] *
                       self.state.q_comp[self.state.ion_type[i_ion]] *
                       (self.state.pore.u2(np.delete(self.state.ion_coord, i_ion, axis=0), r_pos) -
                        self.state.pore.u2(np.delete(self.state.ion_coord, i_ion, axis=0), self.state.ion_coord[i_ion])))
        return d_en * self.state.U_conv_coeff, False

    def d_energy_remove(self, i_ion):
        d_en = (self.state.q_comp[self.state.ion_type[i_ion]]**2) * \
            self.state.pore.u1(self.state.ion_coord[i_ion])
        d_en += np.sum(self.state.q_comp[np.delete(self.state.ion_type, i_ion, axis=0)] *
                       self.state.q_comp[self.state.ion_type[i_ion]] *
                       self.state.pore.u2(np.delete(self.state.ion_coord, i_ion, axis=0), self.state.ion_coord[i_ion]))
        return -1.0 * d_en * self.state.U_conv_coeff

    def d_energy_insert(self, r_pos, i_comp):
        if (self.state.ion_type.size > 0):
            if (self.state.pore.overlap(r_pos, i_comp,
                                         self.state.ion_coord, self.state.ion_type,
                                         self.state.hc_dist2) == 1):
                return 0., True
            d_en = (self.state.q_comp[i_comp]**2) * self.state.pore.u1(r_pos)
            d_en += np.sum(self.state.q_comp[self.state.ion_type] *
                           self.state.q_comp[i_comp] *
                           self.state.pore.u2(self.state.ion_coord, r_pos))
        else:
            d_en = (self.state.q_comp[i_comp]**2) * self.state.pore.u1(r_pos)
        return d_en * self.state.U_conv_coeff, False

    def d_energy_swap(self, i_ion, new_comp):
        if self.state.pore.overlap(self.state.ion_coord[i_ion], new_comp,
                                    np.delete(self.state.ion_coord, i_ion, axis=0),
                                    np.delete(self.state.ion_type, i_ion, axis=0),
                                    self.state.hc_dist2) == 1:
            return 0., True
        d_en = (self.state.q_comp[new_comp]**2 -
                self.state.q_comp[self.state.ion_type[i_ion]]**2) * \
            self.state.pore.u1(self.state.ion_coord[i_ion])
        d_en += np.sum(self.state.q_comp[np.delete(self.state.ion_type, i_ion, axis=0)] *
                       (self.state.q_comp[new_comp] - self.state.q_comp[self.state.ion_type[i_ion]]) *
                       self.state.pore.u2(np.delete(self.state.ion_coord, i_ion, axis=0), self.state.ion_coord[i_ion]))
        return d_en * self.state.U_conv_coeff, False

    def displace_ion(self):
        if self.state.ion_type.size > 0:
            new_ion = np.random.randint(self.state.ion_type.size)
            new_coord = self.state.pore.random_displacement(self.state.ion_coord,
                                                              self.state.ion_type, self.state.step, new_ion, self.state.aion)
            d_en, overlap_bool = self.d_energy_displace(new_ion, new_coord)
            if overlap_bool == False:
                trans_prob = np.exp(-1.0 * d_en)
                if (trans_prob >= 1 or trans_prob > np.random.rand()):
                    self.state.ion_coord[new_ion, :] = new_coord
        return 0

    def remove_ion(self, mu_c_kBT):
        i_comp = np.random.randint(self.state.n_comp)
        if np.sum(self.state.ion_type == i_comp) > 0:
            i_ion = np.random.choice(np.where(self.state.ion_type == i_comp)[0])
            d_en = self.d_energy_remove(i_ion)
            t_w = mu_c_kBT[self.state.ion_type[i_ion]]
            t_prob = np.exp(-1.0 * (t_w + d_en)) * \
                np.sum(self.state.ion_type == i_comp) / \
                self.state.pore.volume_free[i_comp]
            if t_prob >= 1 or t_prob > np.random.rand():
                self.state.ion_coord[i_ion] = self.state.ion_coord[-1]
                self.state.ion_type[i_ion] = self.state.ion_type[-1]
                self.state.ion_coord = self.state.ion_coord[:-1]
                self.state.ion_type = self.state.ion_type[:-1]
                self.state.n_tot = self.state.ion_coord.shape[0]

    def insert_ion(self, mu_c_kBT):
        i_comp = np.random.randint(self.state.n_comp)
        t_w = mu_c_kBT[i_comp]
        new_coord = self.state.pore.random_position(i_comp, self.state.aion)
        d_en, overlap_bool = self.d_energy_insert(new_coord, i_comp)
        if overlap_bool == False:
            t_prob = np.exp(t_w - d_en) * self.state.pore.volume_free[i_comp] / \
                (np.sum(self.state.ion_type == i_comp) + 1)
            if t_prob >= 1 or t_prob > np.random.rand():
                self.state.ion_coord = np.append(self.state.ion_coord, new_coord).reshape(-1, 3)
                self.state.ion_type = np.append(self.state.ion_type, i_comp)

    def swap_ion(self, mu_c_kBT):
        if self.state.ion_type.size > 0:
            from_comp = np.random.choice(np.unique(self.state.ion_type))
            i_ion = np.random.choice(np.where(self.state.ion_type == from_comp)[0])
            to_comp = np.random.choice(np.delete(np.arange(self.state.n_comp), from_comp))
            if self.state.pore.box_bound(self.state.ion_coord[i_ion],
                                          self.state.ion_type[i_ion],
                                          self.state.aion) == False:
                return 0
            d_en, overlap_bool = self.d_energy_swap(i_ion, to_comp)
            if overlap_bool == False:
                t_prob = np.exp(mu_c_kBT[to_comp] - mu_c_kBT[from_comp] - d_en) * \
                    (np.sum(self.state.ion_type == from_comp) *
                     self.state.pore.volume_free[to_comp]) / \
                    ((np.sum(self.state.ion_type == to_comp) + 1) *
                     self.state.pore.volume_free[from_comp])
                if t_prob >= 1 or t_prob > np.random.rand():
                    self.state.ion_type[i_ion] = to_comp

    def mc_step(self, mu_c_kBT):
        rand_x = np.random.rand()
        if rand_x < self.state.mcparams.p_trans:
            self.displace_ion()
        elif rand_x < self.state.mcparams.p_widom:
            if np.random.rand() < 0.5:
                self.remove_ion(mu_c_kBT)
            else:
                self.insert_ion(mu_c_kBT)
        else:
            self.swap_ion(mu_c_kBT)

    def thermalization(self, mu_c_kBT, restart_freq=None, continuation=False):
        for i in np.arange(self.state.mcparams.c_therm, self.state.mcparams.n_therm):
            self.mc_step(mu_c_kBT)
            if continuation and restart_freq and i % restart_freq == 0 and i != 0 and i != self.state.mcparams.c_therm:
                self.state.mcparams.c_therm += restart_freq
                pickle_write(self.state, './cont.restart')
        if continuation:
            self.state.mcparams.c_therm += restart_freq

    def mc_simulate(self, mu_c_kBT, stat_freq, prefix, voltage, coord_dump=None,
                     restart_freq=None, continuation=False):
        for i in np.arange(self.state.mcparams.c_sim, self.state.mcparams.n_sim):
            self.mc_step(mu_c_kBT)

            if i % stat_freq == 0 and not (continuation and i == self.state.mcparams.c_sim):
                self.state.output.dump_count(self.state.ion_type,
                                              self.state.n_comp, prefix, voltage)
                self.state.output.print(i, self.state.ion_type, self.state.n_comp)

            if coord_dump and i % coord_dump == 0:
                self.state.output.dump_coords(self.state.ion_coord,
                                               self.state.ion_type, i, prefix, voltage)

            if continuation:
                if restart_freq and i % restart_freq == 0 and i != 0 and i != self.state.mcparams.c_sim:
                    self.state.mcparams.c_sim += restart_freq
                    pickle_write(self.state, './cont.restart')
            else:
                if restart_freq and i % restart_freq == 0:
                    pickle_write(self.state, f'{prefix}_{voltage}.restart')

        if continuation:
            self.state.mcparams.c_therm = 0
            self.state.mcparams.c_sim = 0
            self.state.mcparams.complete = 1
            pickle_write(self.state, './cont.restart')
            pickle_write(self.state, f'v_{voltage}.restart')
        else:
            pickle_write(self.state, f'{prefix}_{voltage}.restart')


def pickle_write(state_obj, filename):
    with open(filename, 'wb') as fdata:
        pickle.dump(state_obj, fdata, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as fdata:
        return pickle.load(fdata)
