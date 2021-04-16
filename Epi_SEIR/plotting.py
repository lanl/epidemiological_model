"""Plotting function for disease model output.

"""

import matplotlib.pyplot as plt


def graph_model(self):
    """Plots output of WNV model"""
    Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T

    fig = plt.figure(facecolor='w', figsize=[2*6.4, 2*4.8])
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Cases')
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')

    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)

    ax.plot(self.t, Sh, 'b', alpha=0.5, lw=2, label='Susceptible Humans')
    ax.plot(self.t, Iha+Ihs, 'r', alpha=0.5, lw=2, label='Infected Humans')
    ax.plot(self.t, Iv, 'k', alpha=0.5, lw=2, label='Infected Vectors')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    plt.title("WNV Incidence")

    plt.show()
