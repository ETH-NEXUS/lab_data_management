# To initialize the jupyter notebook
# you simply need to put this script to
# /root/.ipython/profile_default/startup/00-init.py

import django
from os import chdir, environ
from sys import path

APP_DIR = "/app"
path.insert(0, APP_DIR)
environ["DJANGO_ALLOW_ASYNC_UNSAFE"] = "true"
chdir(APP_DIR)
django.setup()
