import numpy as np 
import os 


class PolymerSimulation:
    """
    This is the class PolymerSimulation. After specifing the needed system parameters, 
    it will execute the chain.f, or chain_alone.f, code for the polymers input generation, 
    and then it will write the LAMMPS inputs for running. chain.f is for the case of mixing 
    polymers with obstacles, and the chain_alone.f is for when one wants only to simulate 
    polymers in empty space. 


    Parameters
    ----------
    system_parameters
        A dictionary with important system parameters.
    pathto 
        A path to the class directory location.
    """

    def __init__(self,system_parameters,pathto):
        self.system_parameters = system_parameters
        self.pathto = pathto
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
        self.num_files = self.system_parameters["num_files"]
        self.type_simulation = self.system_parameters["type_simulation"]
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
        self.time_step = 0.001 
        self.bin = 1.0 
        self.number_of_steps = self.system_parameters["number_of_steps"]
        self.number_of_steps_equilibration = self.system_parameters["number_of_steps_equilibration"]
        self.low_attraction = self.system_parameters["low_attraction"] 
        self.npart_tot = self.nchain*self.nmonomers + self.n_hs


    def __initialize_polymer_input(self):
        if self.type_simulation:
            os.system("gfortran "+str(self.pathto)+"/lib/chain.f -o chain")

            for i in range(self.num_files):
                with open('def.chain2','w') as fdata:
                    # First line is a comment line 
                    fdata.write('Polymer chain definition\n\n')
                    fdata.write('{}     rhostar\n'.format(self.rho_star))
                    fdata.write('{}     random # seed (8 digits or less)\n'.format(np.random.randint(10000, 100000000)))
                    fdata.write('{}     # of sets of chains (blank line + 6 values for each set)\n'.format(1))
                    fdata.write('{}     molecule tag rule: 0 = by mol, 1 = from 1 end, 2 = from 2 ends\n\n'.format(0))
                    fdata.write('{}     number of chains\n'.format(self.nchain))
                    fdata.write('{}     monomers/chain\n'.format(self.nmonomers))
                    fdata.write('{}     type of monomers (for output into LAMMPS file)\n'.format(1))
                    fdata.write('{}     type of bonds (for output into LAMMPS file)\n'.format(1))
                    fdata.write('{}     distance between monomers (in reduced units)\n'.format(0.97))
                    fdata.write('{}     no distance less than this from site i-1 to i+1 (reduced unit)\n'.format(1.02))
                
            
                os.system("./chain < def.chain2 > "+str(self.filename)+"_poly_input_"+str(i)+".data")
                with open(str(self.filename)+"_poly_input_"+str(i)+".data", "r") as infile:
                    lines = infile.readlines()

                with open(str(self.filename)+"_poly_input_"+str(i)+".data", "w") as outfile:
                    for pos, line in enumerate(lines):
                        if pos != 18 and pos != 20:
                            outfile.write(line)
        else:
            os.system("gfortran "+str(self.pathto)+"/lib/chain_alone.f -o chain")

            for i in range(self.num_files):
                with open('def.chain2','w') as fdata:
                    # First line is a comment line 
                    fdata.write('Polymer chain definition\n\n')
                    fdata.write('{}     rhostar\n'.format(self.rho_star))
                    fdata.write('{}     random # seed (8 digits or less)\n'.format(np.random.randint(10000, 100000000)))
                    fdata.write('{}     # of sets of chains (blank line + 6 values for each set)\n'.format(1))
                    fdata.write('{}     molecule tag rule: 0 = by mol, 1 = from 1 end, 2 = from 2 ends\n\n'.format(0))
                    fdata.write('{}     number of chains\n'.format(self.nchain))
                    fdata.write('{}     monomers/chain\n'.format(self.nmonomers))
                    fdata.write('{}     type of monomers (for output into LAMMPS file)\n'.format(1))
                    fdata.write('{}     type of bonds (for output into LAMMPS file)\n'.format(1))
                    fdata.write('{}     distance between monomers (in reduced units)\n'.format(0.97))
                    fdata.write('{}     no distance less than this from site i-1 to i+1 (reduced unit)\n'.format(1.02))
                
            
                os.system("./chain < def.chain2 > "+str(self.filename)+"_poly_input_"+str(i)+".data")
                with open(str(self.filename)+"_poly_input_"+str(i)+".data", "r") as infile:
                    lines = infile.readlines()

                with open(str(self.filename)+"_poly_input_"+str(i)+".data", "w") as outfile:
                    for pos, line in enumerate(lines):
                        if pos != 18 and pos != 20:
                            outfile.write(line)

    def __initialize_lammps_script(self):
        if self.type_simulation:
            for i in range(self.num_files):
                with open("lammps_"+str(self.filename)+"_"+str(i)+".in", "w") as lmpScript:
                    lmpScript.write("###################\n\n")
                    lmpScript.write("# LAMMPS script generated from class PolymerSimulation\n\n")
                    lmpScript.write("dimension 3 \n")
                    lmpScript.write("atom_style  molecular \n")
                    lmpScript.write("boundary   p p p \n\n")
                    lmpScript.write("neighbor 4.0  multi\n")
                    lmpScript.write("neigh_modify every 2 delay 10 check yes \n\n")
                    lmpScript.write("read_data "+str(self.filename)+"_poly_input_"+str(i)+".data\n\n")
                    lmpScript.write("region box block -{} {} -{} {} -{} {}\n".format(*(self.box_side/2.)*np.ones(6)))
                    lmpScript.write("create_atoms 2 random {} {} box\n\n".format(self.n_hs,np.random.randint(10000, 100000000)))
                    lmpScript.write("mass 1 {}\n".format(self.massA))
                    lmpScript.write("mass 2 {}\n\n".format(self.massB))
                    lmpScript.write("group polymer type 1\n")
                    lmpScript.write("group HS type 2\n\n")
                    lmpScript.write("#########################\n")
                    lmpScript.write("### POLYMER EQUILIBRATION \n")
                    lmpScript.write("pair_style soft 1.0 \n")
                    lmpScript.write("pair_coeff * *  0.0  1.0 \n")
                    lmpScript.write("variable prefactor equal ramp(0,60) \n")
                    lmpScript.write("fix        1   all adapt 1 pair    soft a * * v_prefactor \n")
                    lmpScript.write("bond_style     fene \n")
                    lmpScript.write("bond_coeff 1   30.0    1.5 1.0 1.0 \n ")
                    lmpScript.write("special_bonds fene \n")
                    lmpScript.write("reset_timestep 0 \n")
                    lmpScript.write("timestep   0.001 \n")
                    lmpScript.write("velocity all create 1.0 49589302   rot yes dist gaussian \n")
                    lmpScript.write("fix equilibrate1 all nve \n ")
                    lmpScript.write("fix equilibrate2 all langevin 1.0 1.0 1.0 87708 \n")
                    lmpScript.write("thermo_style  custom step temp pe ke etotal press\n")
                    lmpScript.write("thermo 1000\n")
                    lmpScript.write("run 100000\n")
                    lmpScript.write("unfix 1 \n")
                    lmpScript.write("unfix equilibrate1 \n")
                    lmpScript.write("unfix equilibrate2 \n")
                    lmpScript.write("# stop minimization \n")
                    lmpScript.write("#########################################\n")
                    lmpScript.write("pair_style hybrid/overlay lj/cut {}  lj/cut {}  lj/cut {}\n\n".format(self.cut11,self.cut12,self.cut22))
                    lmpScript.write("pair_coeff      1 1 lj/cut 1 {} {} {}\n".format(self.epsAA,self.sigmaAA,self.cut11))
                    lmpScript.write("pair_modify  shift yes\n\n")
                    lmpScript.write("pair_coeff      1 2 lj/cut 2 {} {} {}\n".format(self.epsAB,self.sigmaAB,self.cut12))
                    lmpScript.write("pair_modify  shift yes\n\n")
                    lmpScript.write("pair_coeff      2 2 lj/cut 3 {} {} {}\n".format(self.epsBB,self.sigmaBB,self.cut22))
                    lmpScript.write("pair_modify  shift yes\n\n")
                    lmpScript.write("bond_style  fene \n")
                    lmpScript.write("bond_coeff 1 30.0  1.5  1.0 1.0\n")
                    lmpScript.write("special_bonds fene\n\n")
                    lmpScript.write("minimize              0.00000001 0.000000001 10000 100000\n")
                    lmpScript.write("reset_timestep 0\n")
                    lmpScript.write("timestep {}\n".format(self.time_step))
                    lmpScript.write("fix integrator all nve\n")
                    lmpScript.write("fix dynamics all langevin {} {} {} {}\n".format(self.temperature,self.temperature,self.gamma,np.random.randint(10000, 100000000)))
                    lmpScript.write("thermo_style  custom step temp pe ke etotal press \n")
                    lmpScript.write("thermo 1000 \n")
                    lmpScript.write("run {}\n".format(self.number_of_steps_equilibration))
                    lmpScript.write("reset_timestep 0\n")
                    lmpScript.write("dump img all custom 10000 simulation."+str(self.filename)+"_"+str(i)+"_snap id type xs ys zs vx vy vz \n")
                    lmpScript.write("compute  cmol all chunk/atom molecule \n")
                    lmpScript.write("compute gyr all gyration/chunk cmol \n")
                    lmpScript.write("variable  ave equal ave(c_gyr)*{} \n".format(self.sigma0))
                    lmpScript.write("variable KineticEnergy equal ke  \n")
                    lmpScript.write("variable PotentialEnergy equal pe \n")
                    lmpScript.write("variable Temperature equal temp  \n")
                    lmpScript.write("fix             output all ave/time 20 50 1000 v_Temperature v_KineticEnergy v_PotentialEnergy v_ave file "+str(self.filename)+"_thermo_output_"+str(i)+".dat format %.15g \n\n")
                    lmpScript.write("thermo_style  custom step temp etotal press v_ave\n")
                    lmpScript.write("thermo 1000 \n")
                    lmpScript.write("compute ficoll_msd HS msd \n")
                    lmpScript.write("fix ficoll_msd HS ave/time 2 6 100 c_ficoll_msd[4] file "+str(self.filename)+"_msd_ficoll_"+str(i)+".dat\n")
                    lmpScript.write("compute cc1 all chunk/atom molecule \n")
                    lmpScript.write("compute myChunk all com/chunk cc1 \n")
                    lmpScript.write("fix myCOM all ave/time 100 1 100 c_myChunk[*] file "+str(self.filename)+"_com_poly_"+str(i)+".dat mode vector\n")
                    lmpScript.write("compute polymer_msd polymer chunk/atom molecule \n")
                    lmpScript.write("compute chainMSD polymer msd/chunk polymer_msd  \n")
                    lmpScript.write("fix polymer_MSDPrint polymer ave/time 1000 1 1000 c_chainMSD[4] file "+str(self.filename)+"_msd_polymer_"+str(i)+".dat mode vector \n")
                    lmpScript.write("run {}".format(self.number_of_steps))
        else:
            for i in range(self.num_files):
                with open("lammps_"+str(self.filename)+"_"+str(i)+".in", "w") as lmpScript:
                    lmpScript.write("###################\n\n")
                    lmpScript.write("# LAMMPS script generated from class PolymerSimulation\n\n")
                    lmpScript.write("dimension 3 \n")
                    lmpScript.write("atom_style  molecular \n")
                    lmpScript.write("boundary   p p p \n\n")
                    lmpScript.write("neighbor 4.0  multi\n")
                    lmpScript.write("neigh_modify every 2 delay 10 check yes \n\n")
                    lmpScript.write("read_data "+str(self.filename)+"_poly_input_"+str(i)+".data\n\n")
                    lmpScript.write("mass 1 {}\n".format(self.massA))
                    lmpScript.write("group polymer type 1\n")
                    lmpScript.write("#########################\n")
                    lmpScript.write("### POLYMER EQUILIBRATION \n")
                    lmpScript.write("pair_style soft 1.0 \n")
                    lmpScript.write("pair_coeff * *  0.0  1.0 \n")
                    lmpScript.write("variable prefactor equal ramp(0,60) \n")
                    lmpScript.write("fix        1   all adapt 1 pair    soft a * * v_prefactor \n")
                    lmpScript.write("bond_style     fene \n")
                    lmpScript.write("bond_coeff 1   30.0    1.5 1.0 1.0 \n ")
                    lmpScript.write("special_bonds fene \n")
                    lmpScript.write("reset_timestep 0 \n")
                    lmpScript.write("timestep   0.001 \n")
                    lmpScript.write("velocity all create 1.0 49589302   rot yes dist gaussian \n")
                    lmpScript.write("fix equilibrate1 all nve \n ")
                    lmpScript.write("fix equilibrate2 all langevin 1.0 1.0 1.0 87708 \n")
                    lmpScript.write("thermo_style  custom step temp pe ke etotal press\n")
                    lmpScript.write("thermo 1000\n")
                    lmpScript.write("run 100000\n")
                    lmpScript.write("unfix 1 \n")
                    lmpScript.write("unfix equilibrate1 \n")
                    lmpScript.write("unfix equilibrate2 \n")
                    lmpScript.write("# stop minimization \n")
                    lmpScript.write("#########################################\n")
                    lmpScript.write("pair_style hybrid/overlay lj/cut {} lj/cut {} \n\n".format(self.cut11, 2.5))
                    lmpScript.write("pair_coeff      1 1 lj/cut 1 {} {} {}\n".format(self.epsAA,self.sigmaAA,self.cut11))
                    lmpScript.write("pair_modify  shift yes\n\n")
                    lmpScript.write("pair_coeff      1 1 lj/cut 2 {} 1.0 2.5 \n\n".format(self.low_attraction))
                    lmpScript.write("bond_style  fene \n")
                    lmpScript.write("bond_coeff 1 30.0  1.5  1.0 1.0\n")
                    lmpScript.write("special_bonds fene\n\n")
                    lmpScript.write("minimize              0.00000001 0.000000001 10000 100000\n")
                    lmpScript.write("reset_timestep 0\n")
                    lmpScript.write("timestep {}\n".format(self.time_step))
                    lmpScript.write("fix integrator all nve\n")
                    lmpScript.write("fix dynamics all langevin {} {} {} {}\n".format(self.temperature,self.temperature,self.gamma,np.random.randint(10000, 100000000)))
                    lmpScript.write("thermo_style  custom step temp pe ke etotal press \n")
                    lmpScript.write("thermo 1000 \n")
                    lmpScript.write("run {}\n".format(self.number_of_steps_equilibration))
                    lmpScript.write("reset_timestep 0\n")
                    lmpScript.write("dump img all custom 10000 simulation."+str(self.filename)+"_"+str(i)+"_snap id type xs ys zs vx vy vz \n")
                    lmpScript.write("compute  cmol all chunk/atom molecule \n")
                    lmpScript.write("compute gyr all gyration/chunk cmol \n")
                    lmpScript.write("variable  ave equal ave(c_gyr)*{} \n".format(self.sigma0))
                    lmpScript.write("variable KineticEnergy equal ke  \n")
                    lmpScript.write("variable PotentialEnergy equal pe \n")
                    lmpScript.write("variable Temperature equal temp  \n")
                    lmpScript.write("fix             output all ave/time 20 50 1000 v_Temperature v_KineticEnergy v_PotentialEnergy v_ave file "+str(self.filename)+"_thermo_output_"+str(i)+".dat format %.15g \n\n")
                    lmpScript.write("thermo_style  custom step temp etotal press v_ave\n")
                    lmpScript.write("thermo 1000 \n")
                    lmpScript.write("compute cc1 all chunk/atom molecule \n")
                    lmpScript.write("compute myChunk all com/chunk cc1 \n")
                    lmpScript.write("fix myCOM all ave/time 100 1 100 c_myChunk[*] file "+str(self.filename)+"_com_poly_"+str(i)+".dat mode vector\n")
                    lmpScript.write("compute polymer_msd polymer chunk/atom molecule \n")
                    lmpScript.write("compute chainMSD polymer msd/chunk polymer_msd  \n")
                    lmpScript.write("fix polymer_MSDPrint polymer ave/time 1000 1 1000 c_chainMSD[4] file "+str(self.filename)+"_msd_polymer_"+str(i)+".dat mode vector \n")
                    lmpScript.write("run {}".format(self.number_of_steps))