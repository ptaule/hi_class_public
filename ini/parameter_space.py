import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

plt.rc("text", usetex=True)

# Make sure to use correct version of class
import sys
sys.path.insert(1, "/home/ptaule/repos/hi_class_public/python/build/lib.linux-x86_64-3.8/")
from classy import Class
from classy import CosmoSevereError, CosmoComputationError

params = {
    "Omega_b"                : 0.05 ,
    "Omega_cdm"              : 0.25 ,
    "output"                 : "mPk" ,
    "input_verbose"          : "1" ,
    "background_verbose"     : "3" ,
    "thermodynamics_verbose" : "1" ,
    "perturbations_verbose"  : "1" ,
    "transfer_verbose"       : "1" ,
    "primordial_verbose"     : "1" ,
    "spectra_verbose"        : "1" ,
    "nonlinear_verbose"      : "1" ,
    "lensing_verbose"        : "1" ,
    "output_verbose"         : "1" ,
}

cosmo = Class()
cosmo.set(params)

mg_params = {
    "Omega_Lambda"           : 0 ,
    "Omega_fld"              : 0 ,
    "Omega_smg"              : -1 ,
    "gravity_model"          : "quasi_static_alphas_power_law",
    "expansion_model"        : "lcdm",
}

gammaM = 0
eta = 3

gammaB_arr = np.linspace(-1,1,31)
gammaT_arr = np.linspace(-1,0,21)

unstable = np.zeros((len(gammaB_arr), len(gammaT_arr)))

for i in range(len(gammaB_arr)):
    for j in range(len(gammaT_arr)):
        print(gammaB_arr[i])
        print(gammaT_arr[j])

        cosmo.set({
            "parameters_smg"  : "%f, %f, %f, %f" % (gammaB_arr[i], gammaT_arr[j], gammaM, eta)
        })

        try:
            cosmo.compute()
        except CosmoComputationError:
            unstable[i][j] = 1
        except CosmoSevereError:
            raise("Something went wrong running class.")
        except KeyboardInterrupt:
            raise("Keyboard interruption.")

        cosmo.struct_cleanup()


x_unstable = gammaB_arr[np.where(unstable==1)[0]]
y_unstable = gammaT_arr[np.where(unstable==1)[1]]

x_stable = gammaB_arr[np.where(unstable==0)[0]]
y_stable = gammaT_arr[np.where(unstable==0)[1]]

fig = plt.figure(constrained_layout=True)

plt.scatter(x_unstable, y_unstable, c="yellow", s=20)
plt.scatter(x_stable, y_stable, c="blue", s=20)

plt.xlabel("$\\gamma_{B}$")
plt.ylabel("$\\gamma_{T}$")

plt.savefig("plots/out.pdf")
plt.clf()
