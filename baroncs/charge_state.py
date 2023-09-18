import os
from typing import Union
import numpy as np
import pandas as pd


class ChargeState:
    energy: float = None
    e0: float = None
    _beta: float = None

    def __init__(self):
        pass

    def charge_state_distribution(
        self,
        atomic_nr: float,
        beta: float = None,
        energy: float = 4.2,
        e0: float = 931.5,
        c_factor: float = 1,
        dist_onesided_len: int = 5,
        plot: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate the charge state distribution for a given atom.

        The formula for the mean charge state and the standard deviation depends on the atomic number.
        for Z < 54, the formula is given by equation 1 and 2 in the paper given as:
        mean_charge_p = Z * (c_factor - exp(-83.28 * (beta / (Z ** (0.447)))))
        std_charge = 0.5 * sqrt(mean_charge * (1 - (mean_charge / Z) ** (1.67)))

        For Z >= 54, the formula is given by equation 3 and 4 in the paper given as:
        mean_charge = mean_charge * (1 - exp(-12.905 + 0.2124 * Z - 0.00122 * Z ** 2))
        std_charge = sqrt(mean_charge_p * (0.07535 + 0.19 * (mean_charge_p / Z) - 0.2654 * (mean_charge_p / Z) ** 2))
        Note that the mean_charge_p is the mean charge state calculated using equation 1. The article always uses
        the mean_charge_p in equation 4.

        :param atomic_nr: The atomic number of the atom
        :param beta: The relativistic beta factor
        :param energy: The kinetic energy
        :param e0: The rest mass energy
        :param c_factor: The c factor
        :param dist_onesided_len: The length of the distribution on one side of the mean.
                                  This will affect the amount of charge states included in the distribution.
                                  If set to 5, the distribution will include 5 charge states on each side of the
                                  mean charge state.
        :param plot: Whether to plot the distribution or not. The plot is made as a bar plot:
                     >>> import matplotlib.pyplot as plt
                     >>> fig, ax = plt.subplots()
                     >>> ax.bar(distx, disty)
                     >>> fig.show()
        :return: The charge state distribution
        """
        self.energy = energy
        self.e0 = e0
        if beta is None:
            self.beta = (energy, e0)
        else:
            self.beta = beta

        mean_charge = self.mean_charge_state(atomic_nr, c_factor)
        std_charge = self.std_charge_state(atomic_nr, c_factor)
        mean_int = round(mean_charge)
        charge_x = (
            np.arange(0, 2 * mean_int + 1)
            if dist_onesided_len > mean_int
            else np.arange(
                mean_int - dist_onesided_len, mean_int + dist_onesided_len + 1
            )
        )
        charge_y = (
            1
            / (std_charge * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((charge_x - mean_charge) / std_charge) ** 2)
        )

        if plot:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
            ax.bar(charge_x, charge_y * 100)
            ax.set_title(f"Charge State Distribution for atom number {atomic_nr}")
            ax.set_ylabel("Percentage")
            ax.set_xlabel("Charge State")
            ax.text(
                s=f"Mean Charge: {mean_charge:.2f},\nStandard Deviation: {std_charge:.2f}",
                x=0.05,
                y=0.95,
                transform=ax.transAxes,
                va="top",
                ha="left",
            )
            fig.savefig("charge_state_distribution.png")
            fig.show()

        return charge_x, charge_y

    def dataframe(
        self,
        beta: int = None,
        e0: float = 931.5,
        energy: float = 4.2,
        atoms: list[int] = "all",
        c_factor: float = 1,
    ) -> pd.DataFrame:
        """
        Generate a dataframe for each atom, which contains the following:
        - Name of the element
        - Atomic Number
        - Commonest Isotope
        - Abundance (%)
        - Melting Point
        - Nearest Q (A/Q=8)
        - Actual Q/A
        - A/Q
        - Mean Charge
        - Most Probable
        - Standard Deviation
        - Probable %
        - Magnetic Rigidity

        If atoms is "all", all atoms will be included in the dataframe. Otherwise provide it as a list of atomic
        numbers. e.g. [54, 82] to only use Xe and Pb.

        If beta is provided it will be used to create the dataframe.
        Otherwise it will be calculated from the energy and e0.

        :param beta: The relativistic factor
        :param e0: The rest mass energy
        :param energy: The kinetic energy
        :param atoms: The atoms to include in the dataframe. If "all", all atoms will be included.
                      Otherwise provide it as a list of atomic numbers. e.g. [54, 82]
        :param c_factor: The c factor
        :return: A dataframe with the charge state distribution for each atom
        """
        df = self.__generate(
            atoms=atoms, beta=beta, e0=e0, energy=energy, c_factor=c_factor
        )
        return df

    def mean_charge_state(
        self, atomic_nr: Union[pd.Series, float], c_factor: float = 1
    ):
        """
        Calculates the mean charge state.

        It is based on either Equation 3 or 1 in the paper. This depends on the atomic number.
        If the atomic number is less than 54, the equation 1 is used. Otherwise equation 3 is used.

        :param atomic_nr: The atomic number
        :param c_factor: The c factor
        :return: The mean charge state
        """
        mean_charge = self.__mean_charge_state_p(atomic_nr, c_factor)  # Equation 1

        if isinstance(atomic_nr, pd.Series):
            # This is the equation 3 in the paper
            mean_charge[atomic_nr >= 54] = self.__mean_charge_state(
                mean_charge[atomic_nr >= 54], atomic_nr[atomic_nr >= 54]
            )
            return mean_charge

        if atomic_nr >= 54:
            # This is the equation 3 in the paper
            mean_charge = self.__mean_charge_state(mean_charge, atomic_nr)
        return mean_charge

    def std_charge_state(self, atomic_nr: Union[pd.Series, float], c_factor: float = 1):
        """
        Calculates the standard deviation of the charge state distribution.

        It is based on either Equation 4 or 2 in the paper. This depends on the atomic number.
        If the atomic number is less than 54, the equation 2 is used. Otherwise equation 4 is used.

        :param atomic_nr: The atomic number
        :param c_factor: The c factor
        :return: The standard deviation of the charge state distribution
        """
        mean_charge = self.__mean_charge_state_p(atomic_nr, c_factor)
        y = mean_charge / atomic_nr
        if isinstance(atomic_nr, pd.Series):
            std_charge = 0.5 * np.sqrt(mean_charge * (1 - y ** (1.67)))
            std_charge[atomic_nr >= 54] = np.sqrt(
                mean_charge[atomic_nr >= 54]
                * (
                    0.07535
                    + 0.19 * y[atomic_nr >= 54]
                    - 0.2654 * y[atomic_nr >= 54] ** 2
                )
            )
            return std_charge

        if atomic_nr >= 54:
            # This is the equation 4 in the paper
            return np.sqrt(mean_charge * (0.07535 + 0.19 * y - 0.2654 * y**2))
        else:
            # This is the equation 2 in the paper
            return 0.5 * np.sqrt(mean_charge * (1 - y ** (1.67)))

    @property
    def beta(self):
        """
        The relativistic beta factor.

        If beta is set as a float, it will be used as the beta factor.
        If beta is set as a tuple of (energy, e0), the beta factor will be calculated from the energy and e0.

        If used together with std_charge_state() and mean_charge_state(), the beta factor can be set as either:
        >>> import baroncs as bcs
        >>> tmp = bcs.ChargeState()
        >>> tmp.beta = 0.09  # Set beta directly
        >>> tmp.beta = (4, 931.5)  # Set beta from energy and e0
        and then the mean charge state and standard deviation can be calculated as:
        >>> tmp.mean_charge_state(82)
        >>> tmp.std_charge_state(82)

        :return: The relativistic beta factor
        """
        if isinstance(self._beta, float):
            return self._beta
        else:
            raise ValueError("Use charge_state_distribution() first or set beta")

    @beta.setter
    def beta(self, value):
        if isinstance(value, float):
            self._beta = value
        elif isinstance(value, tuple):
            self.energy, self.e0 = value
            self._beta = np.sqrt(1 - (1 + self.energy / self.e0) ** (-2))
        else:
            raise ValueError(
                "For the beta setter, provide a tuple of (energy, e0) or a float for beta"
            )

    def __mean_charge_state(
        self, mean_charge: Union[pd.Series, float], atomic_nr: Union[pd.Series, float]
    ):
        return mean_charge * (
            1 - np.exp(-12.905 + 0.2124 * atomic_nr - 0.00122 * atomic_nr**2)
        )

    def __mean_charge_state_p(
        self, atomic_nr: Union[pd.Series, float], c_factor: float = 1
    ):
        """
        Equation 1 in the paper
        """
        return atomic_nr * (
            c_factor - np.exp(-83.28 * (self.beta / (atomic_nr ** (0.447))))
        )

    def __generate(
        self,
        atoms: list[int] = "all",
        beta: int = None,
        e0: float = 931.5,
        energy: float = 4.2,
        c_factor: float = 1,
    ) -> pd.DataFrame:
        data_dir = os.path.join(os.path.dirname(__file__), "table")
        data_path = os.path.join(data_dir, "elements.csv")
        with open(data_path, "r") as f:
            df = pd.read_csv(f, index_col=0)

        df["Nearest Q (A/Q=8)"] = df["Commonest Isotope"].apply(
            lambda x: int(x / 8) if int(x / 8) * 8 == x else int(x / 8) + 1
        )
        df["Actual Q/A"] = df["Nearest Q (A/Q=8)"] / df["Commonest Isotope"]
        df["A/Q"] = df["Commonest Isotope"] / df["Nearest Q (A/Q=8)"]

        # If beta is provided, use it. Otherwise calculate it from the energy and e0
        if beta is None:
            self.beta = (energy, e0)
        else:
            self.beta = beta

        df["Mean Charge"] = self.mean_charge_state(df["Atomic Number"], c_factor)
        df["Most Probable"] = (df["Mean Charge"] + 0.5).astype(int)
        df["Standard Deviation"] = self.std_charge_state(df["Atomic Number"], c_factor)
        df["Probable %"] = (
            100
            / (df["Standard Deviation"] * np.sqrt(2 * np.pi))
            * np.exp(
                -0.5
                * ((df["Most Probable"] - df["Mean Charge"]) / df["Standard Deviation"])
                ** 2
            )
        )
        df["Magnetic Rigidity"] = (
            df["Commonest Isotope"]
            * self.beta
            * 3
            * 1.838
            / (1.759 * df["Most Probable"])
        )

        if atoms == "all":
            return df
        elif isinstance(atoms, list):
            return df[df["Atomic Number"].isin(atoms)]
        else:
            raise ValueError("Provide atoms as a list of atomic numbers")
