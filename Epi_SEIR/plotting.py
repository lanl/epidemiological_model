"""Plotting function for disease model output.

"""
import sys
from utils import create_arg_parser_plot
import matplotlib.pyplot as plt
import pandas as pd


def load_output(disease_name, output_file):
    """Read model output into the program."""
    df = pd.read_csv(output_file)
    human_vec = [x for x in df.columns if "Human" in x or "Time" in x]
    vector_vec = [x for x in df.columns if "Vector" in x or "Time" in x]
    
    if disease_name.lower() == "wnv":
        bird_vec = [x for x in df.columns if "Bird" in x or "Time" in x]
        
    
    human = df[human_vec]
    vector = df[vector_vec]
    
    if disease_name.lower() == "wnv":
        bird = df[bird_vec]
        data = {'human': human, 'vector': vector, 'bird': bird}
    elif disease_name.lower() == "dengue":
        data = {'human': human, 'vector': vector}

    return data


def graph_model():
    """Plots output of human disease model"""
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser_plot()
        args, unknown = parser.parse_known_args()
    
    disease_name = args.disease_name.lower()
    output_file = args.output_file  
    plot_data = load_output(disease_name, output_file)
    
    name1 = output_file.replace('model_output.csv', '')
    name2 = name1.replace('human_model_output/', '')
    
    for i in range(0, len(plot_data.keys())):
        k = plot_data[list(plot_data.keys())[i]]
        n_plot = len(k.columns) -1 
        k.plot(x='Time',subplots=True, figsize=(7.5,n_plot*2.5))
        plt.savefig(f'plots/{name2}{list(plot_data.keys())[i]}.png')
        
    
#     # Sh, Eh, Iha, Ihs, Rh, Sv, Ev, Iv = self.model_output.T
#     Ih, Sh, Rh = load_output()

#     fig = plt.figure(facecolor='w', figsize=[1.5 * 6.4, 4.8])
#     ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#     ax.set_xlabel('Time (days)')
#     ax.set_ylabel('Cases')
#     ax.yaxis.set_tick_params(length=0)
#     ax.xaxis.set_tick_params(length=0)
#     ax.grid(b=True, which='major', c='w', lw=2, ls='-')

#     for spine in ('top', 'right', 'bottom', 'left'):
#         ax.spines[spine].set_visible(False)

#     t = list(range(len(Ih)))
#     ax.plot(t, Sh, alpha=0.5, lw=2, label='Susceptible Humans')
#     ax.plot(t, Ih, alpha=0.5, lw=2, label='Infected Humans')
#     ax.plot(t, Rh, alpha=0.5, lw=2, label='Recovered Humans')
#     legend = ax.legend()
#     legend.get_frame().set_alpha(0.5)
#     plt.title("Dengue Incidence")

#     plt.show()


def main():
    graph_model()


if __name__ == "__main__":
    main()
