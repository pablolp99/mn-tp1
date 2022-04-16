import numpy as np

from pathlib import Path
import yaml

def create_instances(config, fixed_temp=None):
    print(config["INTERNAL_RADIUS"],config["EXTERNAL_RADIUS"],config["M_1_RADIUS"],config["N_DEGREES"],config["ISOTERMA"],config["INSTANCES"])
    
    for _ in range(config["INSTANCES"]):
        external_temp = [config["EXT_TEMP"] for _ in range(config["N_DEGREES"])]
        if fixed_temp is not None:
            internal_temps = [fixed_temp for _ in range(config["N_DEGREES"])]
        else:
            internal_temps = np.random.uniform(low=config["INT_LOW"], high=config["INT_HIGH"], size=config["N_DEGREES"])
        
        print(" ".join(list(map(str, external_temp))), " ".join(list(map(str, internal_temps))))

        
if __name__ == "__main__":
    # Read config
    conf = yaml.safe_load(Path('config.yml').read_text())
    for i in [4,6,8,10,12,14,16,18,20,24,28,32,36,42,48,56,64]:
        conf["M_1_RADIUS"] = i
        conf["N_DEGREES"] = i
        conf["INSTANCES"] = 30
        create_instances(conf, fixed_temp=None)
        print("\n")