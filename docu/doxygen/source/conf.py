import sys
import os
extensions = ['sphinx.ext.pngmath', 'sphinx.ext.todo', 'breathe' ]
breathe_projects = { "GAMERA": os.path.abspath('../xml') }
breathe_default_project = "GAMERA"
