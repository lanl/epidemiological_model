import logging
import VBDM

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')

file_handler = logging.FileHandler('logs/model.log', mode='w')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

logger.info('testing')

def main():
    config_dir = 'config/'
    config_file = config_dir + 'config.yaml'

    den = VBDM.DengueSEIRModel(config_file, days=14)
    den.run_model()
    #den.graph_model()

if __name__ == "__main__":
    main()
