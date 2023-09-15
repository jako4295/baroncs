import os
import numpy as np
import pandas as pd


class ChargeState:
    atomic_nr: float = None
    energy: float = None
    c_factor: float = None
    e0: float = None
    _beta: float = None

    def __init__(self):
        pass

    def charge_state_distribution(
        self,
        atomic_nr: float = None,
        energy: float = 4.2,
        c_factor: float = 1,
        e0: float = 931.5,
        dist_onesided_len: int = 5,
    ) -> tuple[np.ndarray, np.ndarray]:
        self.atomic_nr = atomic_nr
        self.energy = energy
        self.c_factor = c_factor
        self.e0 = e0

        mean_charge = self.mean_charge_state
        std_charge = self.std_charge_state
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

        return charge_x, charge_y

    def dataframe(
        self,
        beta: int = None,
        e0: float = 931.5,
        energy: float = 4.2,
        atoms: int = "all",
    ) -> pd.DataFrame:
        """
        If beta is provided it will be used to create the dataframe.

        If beta is not provided, it will be calculated from the energy and e0. If self.beta is not None, it will be used instead.
        Here is a way to instantiate self.beta:
        >>> tmp = ChargeState()
        >>> tmp.beta = (4.2, 931.5)  # (energy, e0)

        :param beta: The relativistic factor
        :param e0: The rest mass energy
        :param energy: The kinetic energy
        :param atoms: The atoms to include in the dataframe. If "all", all atoms will be included.
                      Otherwise provide it as a list of atomic numbers. e.g. [54, 82]
        :return: A dataframe with the charge state distribution for each atom
        """
        df = self.__generate(atoms, beta, e0, energy)
        return df

    @property
    def mean_charge_state(self):
        """
        Equation 3 or 1 in the paper. This depends on the atomic number.
        """
        if self.atomic_nr is None:
            raise ValueError("Use charge_state_distribution() first")

        mean_charge = self.__mean_charge_state_p()  # Equation 1
        if self.atomic_nr >= 54:
            # This is the equation 3 in the paper
            mean_charge = mean_charge * (
                1
                - np.exp(
                    -12.905 + 0.2124 * self.atomic_nr - 0.00122 * self.atomic_nr**2
                )
            )
        return mean_charge

    @property
    def std_charge_state(self):
        """
        Equation 4 or 2 in the paper. This depends on the atomic number.
        """
        if self.atomic_nr is None:
            raise ValueError("Use charge_state_distribution() first")

        mean_charge = self.__mean_charge_state_p()
        if self.atomic_nr >= 54:
            # This is the equation 4 in the paper
            y = mean_charge / self.atomic_nr
            return np.sqrt(mean_charge * (0.07535 + 0.19 * y - 0.2654 * y**2))
        else:
            # This is the equation 2 in the paper
            return 0.5 * np.sqrt(
                mean_charge * (1 - (mean_charge / self.atomic_nr) ** (1.67))
            )

    @property
    def beta(self):
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

    def __mean_charge_state_p(
        self, atomic_nr: pd.Series = None, c_factor: float = None
    ):
        """
        Equation 1 in the paper
        """
        if atomic_nr is None and self.atomic_nr is None:
            return self.atomic_nr * (
                self.c_factor
                - np.exp(-83.28 * (self.beta / (self.atomic_nr ** (0.447))))
            )
        elif (atomic_nr is not None) and (c_factor is not None):
            return atomic_nr * (
                c_factor - np.exp(-83.28 * (self.beta / (atomic_nr ** (0.447))))
            )
        else:
            raise ValueError(
                "Provide both atomic_nr and c_factor or use charge_state_distribution() first"
            )

    def __generate(
        self,
        atoms: str = "all",
        beta: int = None,
        e0: float = 931.5,
        energy: float = 4.2,
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

        if beta is None and self._beta is None:
            self.beta = (energy, e0)
            beta = self.beta
        elif beta is None:
            beta = self.beta
        else:
            self.beta = beta

        # Change mean charge to use self.mean_charge_state - probably remove property
        # from mean_charge_state and std_charge_state.
        df["Mean Charge"] = self.__mean_charge_state_p(df["Atomic Number"], 1)
        df["Monst Probable"] = (df["Mean Charge"] + 0.5).astype(int)
        df["Standard Deviation"] = 


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    tmp = ChargeState()
    tmp.dataframe()
    hist = tmp.charge_state_distribution(atomic_nr=82)
    plt.plot(hist)
    plt.show()
