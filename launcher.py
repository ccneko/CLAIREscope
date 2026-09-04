#!/usr/bin/env python3
"""CLAIREscope Desktop Launcher Entry Point"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from clairescope.gui import main

if __name__ == "__main__":
    main()
