import numpy as np 
import os 
from random import randint 

class PolymerSimulation:

    def __init__(self,system_parameters):
        self.system_parameters = system_parameters
        self.__initialize_system_parameters()
        self.__initialize_polymer_input()
        self.__initialize_lammps_script()

    def __initialize_system_parameters(self):
        self.sigma0 = self.system_parameters["sigma0"]
        self.mass0 = self.system_parameters["mass0"]
        self.eps0 = self.system_parameters["eps0"]
        self.rho_real = self.system_parameters["rho_real"]
        self.phi_hs = self.system_parameters["phi_hs"]
        self.nchain = self.system_parameters["nchain"]
        self.nmonomers = self.system_parameters["nmonomers"]
        self.filename = self.system_parameters["filename"]
        self.rho_real = self.rho_real / 1.66 
        self.rho_star = (self.rho_real * (self.sigma0**3.)) / self.mass0
        self.sigmaAA = self.sigma0/self.sigma0
        self.sigmaBB = 9.0 / self.sigma0
        self.sigmaAB = (self.sigmaAA + self.sigmaBB) / 2.0 
        self.epsAA = self.eps0 / self.eps0
        self.epsBB = 1.0 
        self.epsAB = 1.0 
        self.cut11 = (2.**(1./6.))*self.sigmaAA
        self.cut22 = (2.**(1./6.))*self.sigmaBB
        self.cut12 = (2.**(1./6.))*self.sigmaAB
        self.massA = self.mass0 / self.mass0 
        self.massB = 70000 / self.mass0
        self.volume = (self.nchain * self.nmonomers) / self.rho_star
        self.box_side = self.volume ** (1./3.)
        self.n_hs = int((6.*self.volume*self.phi_hs)/np.pi/(self.sigmaBB**3.))
        self.temperature = 1.0 
        self.gamma = 1.0 
        self.time_step = 0.01 
        self.bin = 1.0 
        self.number_of_steps = 10000000
        self.number_of_steps_equilibration = 2000000
        self.npart_tot = self.nchain*self.nmonomers + self.n_hs


    def __initialize_polymer_input(self):
        with open('def.chain2','w') as fdata:
            # First line is a comment line 
            fdata.write('Polymer chain definition\n\n')
            fdata.write('{}     rhostar\n'.format(self.rho_star))
            fdata.write('{}     random # seed (8 digits or less)\n'.format(randint(10000, 100000000)))
            fdata.write('{}     # of sets of chains (blank line + 6 values for each set)\n'.format(1))
            fdata.write('{}     molecule tag rule: 0 = by mol, 1 = from 1 end, 2 = from 2 ends\n\n'.format(0))
            fdata.write('{}     number of chains\n'.format(self.nchain))
            fdata.write('{}     monomers/chain\n'.format(self.nmonomers))
            fdata.write('{}     type of monomers (for output into LAMMPS file)\n'.format(1))
            fdata.write('{}     type of bonds (for output into LAMMPS file)\n'.format(1))
            fdata.write('{}     distance between monomers (in reduced units)\n'.format(0.97))
            fdata.write('{}     no distance less than this from site i-1 to i+1 (reduced unit)\n'.format(1.02))
        
        os.system("gfortran lib/chain.f -o chain")
        os.system("./chain < def.chain2 > "+str(self.filename)+"_poly_input.data")
        with open(str(self.filename)+"_poly_input.data", "r") as infile:
            lines = infile.readlines()

        with open(str(self.filename)+"_poly_input.data", "w") as outfile:
            for pos, line in enumerate(lines):
                if pos != 18 and pos != 20:
                    outfile.write(line)

    def __initialize_lammps_script(self):
        with open("lammps_"+str(self.filename)+".in", "w") as lmpScript:
            lmpScript.write("###################\n\n")
            lmpScript.write("# LAMMPS script generated from class PolymerSimulation\n\n")
            lmpScript.write("dimension 3 \n")
            lmpScript.write("atom_style  molecular \n")
            lmpScript.write("boundary   p p p \n\n")
            lmpScript.write("neighbor 1.0  multi\n")
            lmpScript.write("neigh_modify every 2 delay 10 check yes \n\n")
            lmpScript.write("read_data "+str(self.filename)+"_poly_input.data\n\n")
            lmpScript.write("region box block -{} {} -{} {} -{} {}\n".format(*(self.box_side/2.)*np.ones(6)))
            lmpScript.write("create_atoms 2 random {} {} box\n\n".format(self.n_hs,randint(10000, 100000000)))
            lmpScript.write("mass 1 {}\n".format(self.massA))
            lmpScript.write("mass 2 {}\n\n".format(self.massB))
            lmpScript.write("group polymer type 1\n")
            lmpScript.write("group HS type 2\n\n")
            lmpScript.write("pair_style hybrid/overlay lj/cut {}  lj/cut {}  lj/cut {}\n\n".format(self.cut11,self.cut12,self.cut22))
            lmpScript.write("pair_coeff      1 1 lj/cut 1 {} {} {}\n".format(self.epsAA,self.sigmaAA,self.cut11))
            lmpScript.write("pair_modify  shift yes\n\n")
            lmpScript.write("pair_coeff      1 2 lj/cut 2 {} {} {}\n".format(self.epsAB,self.sigmaAB,self.cut12))
            lmpScript.write("pair_modify  shift yes\n\n")
            lmpScript.write("pair_coeff      2 2 lj/cut 3 {} {} {}\n".format(self.epsBB,self.sigmaBB,self.cut22))
            lmpScript.write("pair_modify  shift yes\n\n")
            lmpScript.write("bond_style  fene \n")
            lmpScript.write("bond_coeff 1 30.0  2.25  1.0 1.0\n")
            lmpScript.write("special_bonds fene\n\n")
            lmpScript.write("minimize 1.0e-4 1.0e-6 1000 10000\n")
            lmpScript.write("reset_timestep 0\n")
            lmpScript.write("timestep {}\n".format(self.time_step))
            lmpScript.write("fix integrator all nve\n")
            lmpScript.write("fix dynamics all langevin {} {} {} {}\n".format(self.temperature,self.temperature,self.gamma,randint(10000, 100000000)))
            lmpScript.write("thermo_style  custom step temp pe ke etotal press \n")
            lmpScript.write("thermo 1000 \n")
            lmpScript.write("run {}\n".format(self.number_of_steps_equilibration))
            lmpScript.write("reset_timestep 0\n")
            lmpScript.write("dump img all custom 10000 simulation."+str(self.filename)+"_snap id type xs ys zs vx vy vz \n")
            lmpScript.write("compute  cmol all chunk/atom molecule \n")
            lmpScript.write("compute gyr all gyration/chunk cmol \n")
            lmpScript.write("variable  ave equal ave(c_gyr)*{} \n".format(self.sigma0))
            lmpScript.write("thermo_style  custom step temp etotal press v_ave\n")
            lmpScript.write("thermo 1000 \n")
            lmpScript.write("compute ficoll_msd HS msd \n")
            lmpScript.write("fix ficoll_msd HS ave/time 2 6 100 c_ficoll_msd[4] file "+str(self.filename)+"_msd_ficoll.dat\n")
            lmpScript.write("compute cc1 all chunk/atom molecule \n")
            lmpScript.write("compute myChunk all com/chunk cc1 \n")
            lmpScript.write("fix myCOM all ave/time 100 1 100 c_myChunk[*] file "+str(self.filename)+"_com_poly.dat mode vector\n")
            lmpScript.write("compute cc2 all chunk/atom molecule \n")
            lmpScript.write("compute myChunk2 all msd/chunk cc2\n")
            lmpScript.write("fix myPOLYmsd all ave/time 100 1 100 c_myChunk2[*] file "+str(self.filename)+"_msd_polymers.dat mode vector\n")
            lmpScript.write("run {}".format(self.number_of_steps))