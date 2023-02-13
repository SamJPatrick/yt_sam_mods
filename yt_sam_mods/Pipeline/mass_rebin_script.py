import os
from yt.extensions.sam_mods.model_profiles import create_profile_cube



PROFS_DIR = "Profiles/Normal_profiles"
OUTDIR = '.'



if __name__ == "__main__":
    create_profile_cube(None, output_dir= OUTDIR, data_dir= os.path.join(os.getcwd(), PROFS_DIR))
