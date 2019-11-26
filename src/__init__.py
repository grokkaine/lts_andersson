import os

# We set the directory of the __init__.py file as the current directory
source_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(source_path)
print("Current directory set to:", source_path)

