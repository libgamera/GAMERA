import sys
import os
sys.path.append(os.path.abspath('../breathe'))
extensions = ['sphinx.ext.pngmath', 'sphinx.ext.todo', 'breathe' ]
breathe_projects = { "GAMERA": os.path.abspath('../xml') }
breathe_default_project = "GAMERA"
