import sys

# run from app/ to keep root dir same for tests
from app import view, controller

# If False is used as a command line argument then the remote 
# MTI service is called from controller.py 
if __name__ == "__main__":
	view.run()