if __name__ == '__main__' and __package__ is None:
    #__package__ = get_package()
    print('ikiki')

import os, sys
sys.path.insert(0, os.path.dirname(__file__))
print(sys.path)