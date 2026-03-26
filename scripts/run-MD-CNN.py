import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras.models import load_model
import os
import sys

# load model
path_to_model = "/Users/amos/Documents/repos/amr_phenotype_cnn/models/paper_md_models/MDCNN_saved_model/"

try:
    # Lower-level API - sometimes works when Keras fails
    # model = tf.saved_model.load(
    #     path_to_model, tags=[], export_dir=path_to_model)
    model = tf.keras.models.load_model(path_to_model, compile=False)
    print("✓ Loaded successfully")
    print(type(model))
    print(model.summary())
except Exception as e:
    print(f"✗ Failed: {e}")
