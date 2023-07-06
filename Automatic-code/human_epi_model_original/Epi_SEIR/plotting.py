"""Plotting function for disease model output. 
Typical usage example:

    python plotting.py -of <absolute path model output csv> -d <disease name
    [dengue][wnv]> -sf <save model output figures>
    
    Will automatically show figure with -of and -d, if you -sf it will not automatically show the figure
    Note: naming of saved figures will only work if output data in /human_model_output and 
    files created from models_main.py or models_params.py

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
        if args.save_figure == False:
            plt.show()
            plt.close()
        if args.save_figure == True:
            plt.savefig(f'plots/{name2}{list(plot_data.keys())[i]}.png')
            plt.close()


def main():
    graph_model()


if __name__ == "__main__":
    main()
