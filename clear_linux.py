import os
import shutil

script_dir = os.path.dirname(os.path.abspath(__file__))
build_dir = os.path.join(script_dir, "build")

if os.path.exists(build_dir) and os.path.isdir(build_dir):
    for item in os.listdir(build_dir):
        path = os.path.join(build_dir, item)
        try:
            if os.path.isfile(path) or os.path.islink(path):
                os.unlink(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)
        except Exception as e:
            pass

print("Build directory cleared.")
