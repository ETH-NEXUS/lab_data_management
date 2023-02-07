from django.core.management.base import BaseCommand
from django.core.management import call_command
import traceback
from friendlylog import colored_logger as log
from django.conf import settings
from django.db import connection, connections


class Command(BaseCommand):
  def add_arguments(self, parser):
    parser.add_argument('action', type=str, help='The action to execute')

  def init(self):
    call_command('loaddata', 'plate_dimensions')
    call_command('loaddata', 'well_types')

  def clean(self):
    call_command('flush')
    call_command('makemigrations', 'core')
    call_command('migrate')
    call_command('initadmin')
    self.init()

  def reset(self):
    """
    Resets the whole database.
    The django_extensions package must be installed!
    """
    dbname = settings.DATABASES['default']['NAME']
    answer = input(
        f"Do you really want to drop the database {dbname}? !!! YOU WILL LOSE ALL DATA !!! [y/N] ")
    if answer and len(answer) > 0 and answer[0].lower() == 'y':
      with connection.cursor() as cursor:
        # First disconnect all connection except our own
        cursor.execute(
            f"SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE pid <> pg_backend_pid() AND datname='{dbname}'")
        # disconnect our own connection
        cursor.execute(
            "SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname='{dbname}'")
      # make sure everything is disconnected
      connections.close_all()
      # reset the database DROP / CREATE
      call_command('reset_db', '--no-input')
      self.stdout.write(self.style.SUCCESS(
          f"Successfully reset database '{dbname}'"))

  def handle(self, *args, **options):
    try:
      if options.get('action') == 'clean':
        self.clean()
      elif options.get('action') == 'init':
        self.init()
      elif options.get('action') == 'reset':
        self.reset()
    except Exception as ex:
      log.error(ex)
      traceback.print_exc()
