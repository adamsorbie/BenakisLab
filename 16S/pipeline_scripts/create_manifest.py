import pandas as pd 
import numpy as np
import sys

manifest = pd.read_csv(sys.argv[1], names=["sample-id","forward-absolute-filepath","reverse-absolute-filepath"], sep="\t")

manifest.to_csv(sys.argv[2],index=False, sep="\t")

