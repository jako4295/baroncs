# Table of Contents
- [Table of Contents](#table-of-contents)
- [baroncs ](#baroncs-)
- [Installation ](#installation-)
- [Usage ](#usage-)
  - [Charge state distribution ](#charge-state-distribution-)
  - [Dataframe ](#dataframe-)
  - [Mean charge state and standard deviation ](#mean-charge-state-and-standard-deviation-)
- [Sources ](#sources-)

# baroncs <a name="baroncs"></a>
This package is based on Baron's corrected formula [^1]. It calculates the charge state distribution of ions in a gas target. 

# Installation <a name="installation"></a>
To install the package, run the following command in the terminal:
```bash
pip install baroncs
```
For developer installation, run the following command in the terminal:
```bash
git clone https://github.com/jako4295/baroncs.git
pip install -e baroncs
```

# Usage <a name="usage"></a>

## Charge state distribution <a name="charge_state_distribution"></a>
By default the charge state distribution is using `e0=931.5` and `energy=4.2` to calculate the relativistic beta factor. However, these energies can also be specified. The beta factor can also be specified directly. These examples are shown below:

To obtain the charge state distribution for a given atom and energy, run the following python code:
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
distx, disty = charge_obj.charge_state_distribution(
    atomic_nr=82,  # atomic number of the projectile
    energy=4.2,  # energy of the projectile in MeV/u
    e0=931.5,  # rest energy in MeV
    dist_onesided_len=5,  # length of the on each side of the mean charge state
    plot=True,  # plot the distribution
)
```
This returns the following plot:
![Charge state distribution](plots/charge_state_distribution5.png)

To calculate the charge state distribution given a beta factor, run the following python code:
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
distx, disty = charge_obj.charge_state_distribution(
    atomic_nr=82,  # atomic number of the projectile
    beta=0.09,  # beta factor
    dist_onesided_len=5,  # length of the on each side of the mean charge state
    plot=True,  # plot the distribution
)
```
This returns the following plot:
![Charge state distribution specified beta](plots/charge_state_distribution.png)

If the dist `dist_onesided_len` is changed then we get the following:
```python
distx, disty = charge_obj.charge_state_distribution(
    atomic_nr=82,  # atomic number of the projectile
    energy=4.2,  # energy of the projectile in MeV/u
    e0=931.5,  # rest energy in MeV
    dist_onesided_len=10,  # length of the on each side of the mean charge state
    plot=True,  # plot the distribution
)
```
![Charge state distribution varied dist onesided len](plots/charge_state_distribution10.png)

If it is desired to create the plot in a different way then `distx` and `disty` corresponds to the x and y values of the `plt.bar` function. 

## Dataframe <a name="dataframe"></a>
To create a dataframe with information about each atomic number, run the following python code:
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
df = charge_obj.dataframe()  # this uses default e0=931.5, energy=4.2
```
The dataframe has the following columns:
- `Name`: atomic number of the projectile (this is the index of the dataframe)
- `Atomic Number`: atomic number of the projectile
- `Commonest Isotope`: the most common isotope of the projectile
- `Abundance (%)`: the abundance of the most common isotope
- `Melting Point`: the melting point of the projectile
- `Nearest Q (A/Q=8)`: the nearest charge state to A/Q=8
- `Actual Q/A`: the actual charge state of the projectile
- `A/Q`: the A/Q of the projectile
- `Mean Charge`: the mean charge state of the projectile
- `Most Probable`: the most probable charge state of the projectile
- `Standard Deviation`: the standard deviation of the charge state distribution
- `Probable %`: the probability of the most probable charge state
- `Magnetic Rigidity`: the magnetic rigidity of the projectile

If it is not desired to get all atoms, then we can specify which atoms to get by using the `atoms` argument:
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
df = charge_obj.dataframe(atoms=[54, 82])  # this uses default e0=931.5, energy=4.2
```
This will only return the dataframe for Xenon and Lead (atomic number 54 and 82 respectively).

## Mean charge state and standard deviation <a name="mean_charge_state"></a>
Both the mean charge state and the standard deviation uses the relativistic beta factor. If `charge_obj.charge_state_distribution` is called then the relativistic beta function is calculated and stored in `charge_obj.beta` (this is the one that will be used if `charge_obj.mean_charge_state` or `charge_obj.std_charge_state` is called). 

If it is desired to calculate the mean charge state and standard deviation with a different beta factor, then it can be specified either using `e0` and `energy`, or specifying the beta factor directly. These examples are shown below:
*Using `e0` and `energy`*
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
charge_obj.beta = (4.2, 931.5)  # set the beta factor (energy, e0)
mean_charge_state = charge_obj.mean_charge_state(
    atomic_nr=82,  # atomic number of the projectile
)
std_charge_state = charge_obj.std_charge_state(
    atomic_nr=82,  # atomic number of the projectile
)
```
*Specifying `beta` directly*
```python
import baroncs as bcs

charge_obj = bcs.ChargeState()
charge_obj.beta = 0.09  # set the beta factor directly
mean_charge_state = charge_obj.mean_charge_state(
    atomic_nr=82,  # atomic number of the projectile
)
std_charge_state = charge_obj.std_charge_state(
    atomic_nr=82,  # atomic number of the projectile
)
```


# Sources <a name="sources"></a>
[^1]: [Charge exchange of very heavy ions in carbon foils and in the residual gas of GANIL cyclotrons](https://www.sciencedirect.com/science/article/pii/016890029390622O)