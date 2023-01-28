# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

import os
import datetime

# DATA DIR prefix
DATA_DIR_PREFIX = "data"
CIF_DIR_PREFIX = "cifs"
DATA_DIR = os.path.join(__file__, DATA_DIR_PREFIX)
CIF_DIR = os.path.join(__file__, CIF_DIR_PREFIX)

# Default input arguments
INT_THRESHOLD = 1.0e-5
DEFAULT_SYMPREC = 1.0e-2
DEFAULT_ANGLE_TOLERANCE = 5.0

# Global state configuration
DISABLE_PROGRESSBAR = False

# Display formatting
LINEWIDTH = 88
HLINE = "-" * LINEWIDTH
TQDM_CONF = {"ncols": LINEWIDTH}

# Date and time, also used for (default) random seed
NOW = datetime.datetime.now()
DEFAULT_SEED = int(NOW.timestamp())
