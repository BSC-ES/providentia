import os

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

# function that joins paths making sure they are all in the right direction
def join(*args):
    return os.path.join(*args).replace("\\","/")
