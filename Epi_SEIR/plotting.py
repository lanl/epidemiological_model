"""Plotting function for disease model output.

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def load_output():
    df = pd.read_csv("human_model_output/dengue_model_output.csv")
    Iha = df['Asymptomatic Infected Humans'].tolist()
    Ihs = df['Symptomatic Infected Humans'].tolist()

    return Iha, Ihs


def graph_model():
    """Plots output of WNV model"""
    # Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T
    Iha, Ihs = load_output()

    fig = plt.figure(facecolor='w', figsize=[1.5 * 6.4, 4.8])
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Cases')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')

    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    t = list(range(len(Iha)))
    ax.plot(t, Iha, 'r', alpha=0.5, lw=2, label='Asymptomatic Infected Humans')
    ax.plot(t, Ihs, 'b', alpha=0.5, lw=2, label='Symptomatic Infected Humans')
    print(np.array(Ihs) - np.array(Iha))
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    plt.title("Dengue Incidence")

    plt.show()


def main():
    graph_model()


if __name__ == "__main__":
    main()
