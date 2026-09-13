"""
The threshold that decides which library wells count as running low.
"""

from rest_framework import viewsets
from rest_framework.permissions import IsAuthenticated


class ThresholdViewSet(viewsets.ModelViewSet):
    """
    Lets every logged in user read and change the one threshold, for example
    from the form on the messages page. Creating a second threshold or deleting
    the only one is not possible through the API: without it no well would be
    judged against the values the users set. The admin can still do both.
    The endpoint adds the queryset and the serializer.
    """

    permission_classes = [IsAuthenticated]
    # PUT is left out as well: the page only ever changes single values (PATCH).
    http_method_names = ["get", "patch", "head", "options"]
