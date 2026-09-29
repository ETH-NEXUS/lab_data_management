"""
The result zips of an experiment are deleted with the experiment, so nothing is
left in the media folder that no page shows any more.
"""

import os
import shutil

from django.db.models.signals import post_delete
from django.dispatch import receiver

from analysis import tasks
from core.models import Experiment


@receiver(post_delete, sender=Experiment)
def delete_analysis_results(sender, instance, **kwargs) -> None:
    """
    Called for every deleted experiment, also when its project is deleted.
    Removes <MEDIA_ROOT>/analysis/<experiment id>/ with all its zips.
    """
    experiment_folder = os.path.join(tasks.ANALYSIS_FOLDER, str(instance.id))
    # An experiment that was never analysed has no folder
    shutil.rmtree(experiment_folder, ignore_errors=True)
